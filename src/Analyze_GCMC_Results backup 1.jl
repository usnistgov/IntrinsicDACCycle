"""Function to convert temperature to thermodynamic β."""
function T_to_β(x)
    #input x in Kelvin
    #return β in mol/kJ
    β = (1. / x) * (1. / kB) * (1. / Na) *1000.
    return β
end

"""Function to convert thermodynamic β to temperature."""
function β_to_T(x)
    #input x in mol/kJ
    #return T in Kelvin
    T = (1. / x) * (1. / kB) * (1. / Na) *1000.
    return T
end

"""This function reads in the JSON files for the GCMC results as well as the material discription file."""
function read_jsons(directory, name)
    #strip off the "_clean" suffix and copy that name
    name_parsed = replace(name, "_clean" => "")

    # material_string = directory*name*"_clean.json"
    material_string = directory*"/CSD_FEASST_Materials/Materials/"*name*".json"
    material = JSON.parsefile(material_string)

    # Kh_N₂string = directory*"results_300_N2_"*name*".json"
    Kh_N₂string = directory*"/Results/"*"results_300_N2_"*name_parsed*".json"
    Kh_N₂ = JSON.parsefile(Kh_N₂string)

    # Kh_CO₂string = directory*name*"_clean.results.json"
    Kh_CO₂string = directory*"/CSD_FEASST_Materials/Results/"*name*".results.json"
	Kh_CO₂ = JSON.parsefile(Kh_CO₂string)

	# One_atm_string = directory*"results_300K_1atmN2_"*name*".json"
	One_atm_string = directory*"/Atmosphere_check/"*"results_300K_1atmN2_"*name_parsed*".json"
	One_atm_N₂ = JSON.parsefile(One_atm_string)

    return material, Kh_N₂, Kh_CO₂, One_atm_N₂
end

"""This function takes the predicted Henry constant from the GMCC results
 and extrapolates it to new thermodyanmic β."""
function Kh_extrapolate(β, Kh_results, material)
    #Extrapolate the Kh to a new beta using the Kh* coeffs.

    #Define the factor for the t-distribution for standard deivation
    tvalue = 1.959
    
    #get the material density
    ρₛ = material["cell_density"] #Densitiy in amu/angstrom^3
    ρₛ *= amu * 1e30 #convert density to kg/m^3
    
    #Find the beta where Kh was calculated
    β₀ = Kh_results["beta"] #[mol/kJ]
    T = β_to_T(β₀) # [K]

    Δβ = abs.(β .- β₀) #[mol/kJ]
    if any(Δβ .> 0.3)
        print("Warning: Check Extrapolation Range")
        print("GCMC at ", T, " K")
        print("Extrapolating to", β_to_T.(β))
    end

    δᵦ = β .- β₀ #Distance in beta-space [mol/kJ]

    Kcoeff = Kh_results["Kcoeff"] #K* coefficients
    power = length(Kcoeff)
    powers = collect(range(0, power-1, power)) #the polynomial degree for each coeff
    reps = length(δᵦ)
    power_reps = ones((reps,1)) .* powers' #Matrix of each degree for each δᵦ
    δᵦ_m =  δᵦ .* ones(power)' 
    # δᵦᴾ = δᵦ' .^ power_reps #each δᵦ to each power.
    δᵦᴾ = δᵦ_m .^ power_reps #each δᵦ to each power.
    kᵦ = δᵦᴾ .* Kcoeff'
    x = sum(kᵦ, dims = 2) #Sum the terms of the polynomial for each δᵦ

    Kh = (β./ρₛ).*x #Kh [mmol/(kg Pa)] Mulitply by the pre-factor

    Covariance = hcat(Kh_results["Kcoeff_covar"]...)
    Jacobian = (reshape(β, (:,1))./ ρₛ) .* δᵦᴾ
    
    poly_var = Jacobian*Covariance*Jacobian' #the whole covariance matrix for the results
    trials = Kh_results["trials"]
    Kh_std = reshape(sqrt.(diag(poly_var)./trials)*tvalue, (:,1)) #get the std from the covariance matrix

    return Kh, Kh_std # [mmol/(kg Pa)]
end

"""This function takes the predictions of the heat of adsorption at infinate dilution,
and extrapolates to new thermodyanmic β."""
function qₐ∞(β, Kh_results)
    #Calculate the heat of adsorption at infinate dilution
    #Define the factor for the t-distribution for standard deivation
    tvalue = 1.959
    
    β₀ = Kh_results["beta"] #[mol/kJ]

    δᵦ = β .- β₀ #Distance in beta-space

    Kcoeff = Kh_results["Kcoeff"] #K* coefficients  
    power = length(Kcoeff)
    powers = collect(range(0, power-1, power)) #the polynomial degree for each coeff
    reps = length(δᵦ)
    power_reps = ones((reps,1)) .* powers' #Matrix of each degree for each δᵦ
    δᵦ_m =  δᵦ .* ones(power)' #matrix of each  δᵦ
    # δᵦᴾ = δᵦ' .^ power_reps #each δᵦ to each power.
    δᵦᴾ = δᵦ_m .^ power_reps #each δᵦ to each power.

    kᵦ = δᵦᴾ .* Kcoeff'

    denominator = sum(kᵦ, dims = 2) #Sum the terms of the polynomial for each δᵦ
    
    Kcoeff_1 = Kcoeff[2:end]
    powers_1 = powers[2:end]
    power_reps_1 = power_reps[:, 2:end]
    δᵦ_m_1 =  δᵦ .* ones(power-1)'
    δᵦᴾ_1 = δᵦ_m_1 .^(power_reps_1 .- 1) 
    kᵦ_1 = δᵦᴾ_1 .* Kcoeff_1' .* powers_1'

    numerator = sum(kᵦ_1, dims =2)

    q_ads_∞ = 1 ./β .+ (numerator ./ denominator)

    #Uncertainty Propagation:
    #Add a leading zero to the k_beta for df/dKh0
    
    kᵦ_1_0 = hcat(reshape(0 .* β', (length(β), 1)), kᵦ_1)
    A = sum(kᵦ_1_0, dims = 2)
    B = sum(kᵦ, dims = 2)
    
    first_term = power_reps .* B .* δᵦ_m .^(power_reps .-1)
    # @show first_term
    second_term = A .* δᵦ_m .^(power_reps)
    # @show size(second_term)
    Jacobian = (first_term .- second_term) ./(B.^2)
    #when δᵦ = 0 use q = 1/β + Kh2/Kh1
    #^ in which case the Jacobian becomes:
    Jacobian_0 = zeros(power)
    Jacobian_0[1] = -Kcoeff[2]/(Kcoeff[1]^2)
    Jacobian_0[2] = 1/Kcoeff[1]
    #find indexes where δ_beta == 0
    mask_δᵦ = findall(x -> x == 0, reshape(δᵦ, (:)))
    Jacobian[mask_δᵦ, CartesianIndex.(1:power)] .= reshape(Jacobian_0, (1,:)) 
     

    Covariance = hcat(Kh_results["Kcoeff_covar"]...)
    q_ads_∞_var = Jacobian*Covariance*Jacobian' #The whole covariance matrix for the results
    trials = Kh_results["trials"]
    
    q_ads_∞_std = reshape(sqrt.(diag(q_ads_∞_var)./trials)*tvalue, (:,1)) #get the std from the covariance matrix
    
    return q_ads_∞, q_ads_∞_std #kJ/mol of gas
end

"""This function compares the adsorption as predicted by the extrapolated Henry constant
to the adsorption as predicted by the direct GCMC simulation."""
function Uptake_Kh_v_MC(material, Kh_results, MC_results)
    #Function to compare the adsorption from the extrapolated Kh
    # to the directly simulated GCMC adsorption
    
    #Define the factor for the t-distribution for standard deivation
    tvalue = 1.959
    
    molar_mass_supercell = material["cell_mass"] #[g/mol] for sorbent

    #find the β of the GCMC
    # β_MC = convert(Array{Float64,1},MC_results["beta"]) #[mol/kJ]
    β_MC = Array([parse(Float64,MC_results["beta"])]) #[mol/kJ]
    
    Kh_β, σ_Kh_β = Kh_extrapolate(β_MC, Kh_results, material) #[mmol/(kg Pa)]

    target_P = MC_results["target_p"] #[kJ/(mol Å³)]
    target_P *= 1000/(Na*1e-30) #[Pa]

    #Calculate the uptake from the extrapolated Kh
    Uptake_Kh = Kh_β[1] * target_P #[mmol/kg sorbent]
    σ_Uptake_Kh = σ_Kh_β[1] * target_P #[mmol/kg sorbent]

    #Calculate the uptake from the GCMC
    adsorb_MC = MC_results["molecules_adsorptives0_average"] #[#N₂/Supercell]
    Uptake_MC = adsorb_MC*1e6/molar_mass_supercell #[mmol/kg sorbent]
    σ_adsorb_MC = MC_results["molecules_adsorptives0_block_std"]
    σ_Uptake_MC = σ_adsorb_MC*1e6*tvalue/molar_mass_supercell #[mmol/kg sorbent]

    return (Uptake_Kh, σ_Uptake_Kh), (Uptake_MC, σ_Uptake_MC)
    
end

"""This function determines if the extrapolated Henry constant 
is close enough to the direct GCMC simulation to safely use."""
function Close_enough(material, Kh_results, MC_results)
    #Calculate the uptake from Kh extrapolation and GCMC
    Uptake_Kh, Uptake_MC = Uptake_Kh_v_MC(material, Kh_results, MC_results)

    μ_Uptake_Kh, σ_Uptake_Kh = Uptake_Kh
    μ_Uptake_MC, σ_Uptake_MC = Uptake_MC
    

    #Factor for agreement
    factor = 0.1
    #Limits for the Kh result
    lower_Kh = μ_Uptake_Kh*(1-factor)
    upper_Kh = μ_Uptake_Kh*(1+factor)
    #Limits for the MC result
    lower_MC = μ_Uptake_MC - 2*σ_Uptake_MC
    upper_MC = μ_Uptake_MC + 2*σ_Uptake_MC
    
    #Is the lower bound of the MC result within a factor of the bounds of the Kh result?
    test1 = (lower_MC > lower_Kh) & (lower_MC < upper_Kh)

    #Is the mean of the MC result within a faftor of the bounds of the Kh result?
    test2 = (μ_Uptake_MC > lower_Kh) & (μ_Uptake_MC < upper_Kh)

    #Is the upper bound of the MC result within a factor of the bounds of the Kh result?
    test3 = (upper_MC > lower_Kh) & (upper_MC < upper_Kh)

    tests = any([test1, test2, test3])
    return tests 
end
