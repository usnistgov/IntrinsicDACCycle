### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# ╔═╡ 6c576153-4bca-411a-b745-3a3cb88d6ad9
using Pkg

# ╔═╡ a7a44fe8-29c9-4b65-93f0-ae2672eb1960
Pkg.activate("/users/asm6/Julia_scripts/IntrinsicDACCycle")

# ╔═╡ d177558c-ff0e-4089-a79b-11c27007c759
begin
	using Revise
	using IntrinsicDACCycle
end

# ╔═╡ 7927d0de-3ca7-408d-bdb2-93d252af9a80
using Plots

# ╔═╡ 5ee89c88-5efe-40e2-bf24-af8facfa2199
using JSON

# ╔═╡ d31a012a-515c-4521-9e72-600eb36b6083
using LinearAlgebra

# ╔═╡ 2012b243-b921-4c28-a31c-74969c46276f
using Metaheuristics

# ╔═╡ 316b64c9-029d-456b-b1b9-ec82a14a7f1b
begin
	using NaNStatistics
	using Distributed
	using Distributions
end

# ╔═╡ c6ccef8c-4db6-11ef-0f54-0f3bae2f93b4
cd("/users/asm6/Julia_scripts/IntrinsicDACCycle")

# ╔═╡ 3c8fa7bc-b4c2-41cd-b050-a98e3407b43e
begin
	directory = "/users/asm6/DAC_data"
	name = "OKILEA_clean"
	# name = "ABAVIJ_clean"
	# name = "ROLCEC13_clean"

	material, Kh_N₂, Kh_CO₂, One_atm_N₂ = IntrinsicDACCycle.read_jsons(directory,name)
	

end

# ╔═╡ ce51f824-77fa-402a-a16e-372d5f61838b
Kh_CO₂

# ╔═╡ f6e3c4b8-b4dd-463e-8246-c49305237a31
begin
	Ts = collect(LinRange(200, 400, 5))
	Ps = collect(LinRange(1.1*101325, 1.0*101325, 5))

	βs = IntrinsicDACCycle.T_to_β.(Ts)

	# Kh_CO2s = IntrinsicDACCycle.Kh_extrapolate(βs, Kh_CO₂, material)
	# q = IntrinsicDACCycle.qₐ∞(β, Kh_CO₂)
end

# ╔═╡ edea02bf-4658-40b4-b11f-2faabf1a7200
begin
	α = 400/1000000
	# thing = IntrinsicDACCycle.Optimize_Intrinsic_Refresh_path_distributions(directory, name, α)
end

# ╔═╡ 1420ac5b-c8ba-4ccd-a047-9a21ba6b4b97
begin
	
	saturation_string = directory*"/Saturation/Saturation/CO2_sat_200K_"*name*".json"
	saturation_dict = JSON.parsefile(saturation_string)

	step = saturation_dict["step_number"]
	loading = saturation_dict["isotherm"]["$step"]["loading"] #[molecules CO2]
    loading_stdev = saturation_dict["isotherm"]["$step"]["stdev"] #[molecules CO2]

	material_string = directory*"/CSD_FEASST_Materials/Materials/"*name*".json"
	material_dict = JSON.parsefile(material_string)
	mass = material_dict["cell_mass"] #[amu] or [g/mol]

	sat_loading = loading*1e6/mass # [mmol/kg] millimoles of CO2 per kg of sorbent

	sat_Kh = sat_loading/(α * Ps[1])
end 

# ╔═╡ ea1ccc4d-4a86-4d53-b865-1c7986d3598a
function saturation_adsorb_co2(directory, name)
	#read in the saturation file
	saturation_string = directory*"/Saturation/Saturation/CO2_sat_200K_"*name*".json"
	saturation_dict = JSON.parsefile(saturation_string)

	#Read the loading at saturation
	step = saturation_dict["step_number"]
	loading = saturation_dict["isotherm"]["$step"]["loading"] #[molecules CO2]
    loading_stdev = saturation_dict["isotherm"]["$step"]["stdev"] #[molecules CO2]

	#Test if the loading is larger than 0.9 molecules
	test = loading > 0.9 
	return test
end
	

# ╔═╡ 80cd44e7-af7c-4751-9e78-77d94a82f32b
saturation_dict

# ╔═╡ 1d569789-5ddb-4e14-a063-4103cfb75de1
sat_Kh

# ╔═╡ 62606add-5fc3-4518-a88f-4bd0b253ce39
sat_loading

# ╔═╡ cbbe86a4-06cf-4344-a2db-a90980ccf481
begin 
	loadings = []
	fugacities = []
	for steps in range(1, saturation_dict["step_number"])
		push!(loadings, saturation_dict["isotherm"][string(steps)]["loading"])
		push!(fugacities, saturation_dict["isotherm"][string(steps)]["fugacity"])
	end
end

# ╔═╡ 7b17395c-7cc2-4dc0-a361-17114e08bbec
begin
	plot(fugacities, loadings )
	scatter!(fugacities, loadings)
end

# ╔═╡ 1e37e81f-6d2f-4b9a-815f-6f8699ab7083
begin
	plot(fugacities, loadings ./ fugacities)
	scatter!(fugacities, loadings ./ fugacities)
end

# ╔═╡ 00cd9198-483e-4d7f-9add-11596183583b


# ╔═╡ 63b7a1fe-b67a-4d15-a5b5-486723daaaf5
begin
	#Define useful constants
	const kB = 1.380649e-23  #J/K
	const Na = 6.02214076e23   #1/mol
	const amu = 1.66053906660e-27 #kg
	const R = Na*kB #[J/(K mol)]
end 

# ╔═╡ 63d7bd63-181a-4a99-84b6-2931a468726b
function Kh_extrapolate(β, Kh_results, material)
    #Extrapolate the Kh to a new beta using the Kh* coeffs.

    #Define the factor for the t-distribution for standard deivation
    tvalue = 1.959
    
    #get the material density
    ρₛ = material["cell_density"] #Densitiy in amu/angstrom^3
    ρₛ *= amu * 1e30 #convert density to kg/m^3
    
    #Find the beta where Kh was calculated
    β₀ = Kh_results["beta"] #[mol/kJ]
    T = IntrinsicDACCycle.β_to_T(β₀) # [K]

    Δβ = abs.(β .- β₀) #[mol/kJ]
    if any(Δβ .> 0.3)
        print("Warning: Check Extrapolation Range")
        print("GCMC at ", T, " K")
        print("Extrapolating to", IntrinsicDACCycle.β_to_T.(β))
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

# ╔═╡ 159eee30-a3c4-4a9d-848e-3f7c2f75f7cf
Kh_CO2s, Kh_std = Kh_extrapolate(βs, Kh_CO₂, material)

# ╔═╡ 69c4ed3f-dd56-4a3a-9740-704a10ce07f5
Kh_N2s, Kh_N2std = Kh_extrapolate(βs, Kh_CO₂, material)

# ╔═╡ fe759187-2198-48d4-be35-309608b16451
test = IntrinsicDACCycle.Kh_extrapolate(βs, Kh_CO₂, material)

# ╔═╡ 78a2b7b2-3299-4ab5-b99b-97d37ac23d9c
Kh_CO2s .* α .* Ps 

# ╔═╡ fe85ba0d-68de-46ba-a6de-434bf81af8b4
Kh_CO2s .* α .* Ps .< sat_loading

# ╔═╡ f3c54cae-0ce5-48e3-883a-b64da8818050
IntrinsicDACCycle.Close_enough(material, Kh_N₂, One_atm_N₂)

# ╔═╡ c7004efc-4d62-4ab8-a787-e25f9e84c706
Ps

# ╔═╡ d3e9d752-dcff-4486-93dd-82d68a6be81e
function truncate_to_saturation(directory, name, α,
                                Ts, Ps, βs,
                                Henry_CO2, Henry_CO2_err,
                                Henry_N2, Henry_N2_err)

    #Read in the Saturatuion JSON
    saturation_string = directory*"/Saturation/Saturation/CO2_sat_200K_"*name*".json"
    saturation_dict = JSON.parsefile(saturation_string)

    #Find the saturation CO2 molecules per unit cell
    step = saturation_dict["step_number"]
    loading = saturation_dict["isotherm"]["$step"]["loading"] #[molecules CO2]
    loading_stdev = saturation_dict["isotherm"]["$step"]["stdev"] #[molecules CO2]

    #Read in the unit cell info
    material_string = directory*"/CSD_FEASST_Materials/Materials/"*name*".json"
    material = JSON.parsefile(material_string)
    mass = material["cell_mass"] #[amu] or [g/mol]

    #calculate the loading at saturation
    sat_loading = loading*1e6/mass # [mmol/kg] millimoles of CO2 per kg of sorbent

    #calculate the loadings from Henry Isotherms
    loadings_Kh = Henry_CO2 .* α .* Ps

    #Mask where the loadings from Henry Isotherms would be bigger than saturation
    mask = loadings_Kh .< sat_loading
    # #calculate the straight line isotherm to saturation from starting partial pressure
    # sat_Kh = sat_loading/(α * Ps[1])

    # #Mask off anywhere the Henery Constant is bigger than the straight line isotherm to saturation
    # mask = Henry_CO2 .< sat_Kh

    #Use the mask to truncate the vectors
	mask = reshape(mask, :)
    Ts = Ts[mask]
    Ps = Ps[mask]
    βs = βs[mask]

    Henry_CO2 = Henry_CO2[mask]
    Henry_CO2_err = Henry_CO2_err[mask]

    Henry_N2 = Henry_N2[mask]
    Henry_N2_err = Henry_N2_err[mask]

    return Ts, Ps, βs, Henry_CO2, Henry_CO2_err, Henry_N2, Henry_N2_err
end

# ╔═╡ 49169e8a-868d-42d7-9fc7-5c7b9eaea90c
sat_test = IntrinsicDACCycle.truncate_to_saturation(directory, name, α, Ts, Ps, βs, Kh_CO2s, Kh_std, Kh_N2s, Kh_N2s)

# ╔═╡ 438c46d9-5872-4c58-9902-de159e90f16a
begin
	thingamabob = IntrinsicDACCycle.Intrinisic_refresh_objectives_posterior_dist(directory, name,
																					 Ts, Ps, α, 100)
end

# ╔═╡ b6f58e4f-2072-4d55-9be1-be4f357f74fa
Ts

# ╔═╡ 43853edf-20bd-4bb6-bdf0-b60a0aee225e
anothertest = IntrinsicDACCycle.Optimize_Intrinsic_Refresh(directory, name, α)

# ╔═╡ 59b733d3-e361-423a-a010-966fe16c8176
begin
	function ScorePath(parameters, dT_max, dP_max, 
	                   Base_directory, name, α,  
	                   return_flag = "objectives", i = 0)
	
	    T_start = parameters[1]
	    ΔT = parameters[2]
	    P_start = parameters[3]
	    ΔP = parameters[4]
	    
	    outputs = try
	        
	        #Calculate the number of steps from start to start+Δ with the max T step or P step size
	        T_steps = length(T_start:dT_max:T_start+ΔT)
	        P_steps = length(P_start:dP_max:P_start+ΔP)
	
	        #Choose the larger number of steps
	        steps = maximum([T_steps, P_steps])
	        #Create the T and P path with those steps
	        Ts = collect(LinRange(T_start, T_start+ΔT, steps))
	        Ps = collect(LinRange(P_start, P_start+ΔP, steps))

	        #Perform the Intrinsic Refresh Analysis with those steps
	        ξ, α_end =  IntrinsicDACCycle.Intrinisic_refresh_objectives(Base_directory, name,
	                                                                Ts, Ps, α)
			# if ismissing(ξ)
			# 	ξ = NaN
			# end
			# if ismissing(α_end)
			# 	α_end = NaN
			# end
	        (ξ, α_end, dT_max, dP_max, i)
	    catch e
	        if isa(e, DomainError) 
	            #If there is a DomainError it's likely from a numerical instabilty with IAST 
	            # So take a smaller step size and try again. Do this up to 10 times
	            if i <= 10
	                dT_max *= 0.75
	                dP_max *= 0.75
	                i += 1
	                println("Re-trying with smaller steps!")
	                ScorePath(parameters, dT_max, dP_max, Base_directory, name, α, "outputs", i)
	            else
	                println("Too many iterations, likely something other than stepsize is wrong")
	                throw(e)
	            end
	        end
	    end

	    ξ = outputs[1]
	    α_end = outputs[2]
		if ismissing(ξ)
			ξ = NaN
		end
		if ismissing(α_end)
			α_end = NaN
		end
		
	    objectives = [1/ξ, 1-α_end]
	    gx = [0.0] # inequality constraints
	    hx = [0.0] # equality constraints
	
	    if return_flag == "outputs"
	        return outputs
	    elseif return_flag == "objectives"
	        return objectives, gx, hx
	    elseif return_flag == "steps"
	        dT_max = outputs[3]
	        dP_max = outputs[4]
	        return objectives, gx, hx, dT_max, dP_max
	    else 
	        return "Need to define return_flag"
	    end
	
	end
end 

# ╔═╡ a8700a52-51a2-49d6-8c22-cf331a54731e
begin
	t1 = false
	t2 = false

	t3 = t1 & t2
end

# ╔═╡ b0405caf-1431-4c04-bb7c-56fe9cc468c8
begin
	
"""Function to optimize the start and end temperatures and pressures 
for a given material and inlet CO2 concentration.
Assumes a linear path,"""
function Optimize_Intrinsic_Refresh(Base_directory::String, name::String,  
                                    α::Real)
	#Initialize a dictionary to store all the results
	Results_Dict = sort(Dict{String, Any}("Name" => name))
    
	#do the close_enough to linear adsorption test
	#and the saturation test before attempting the calculation
	material, Kh_N₂, Kh_CO₂, One_atm_N₂ = IntrinsicDACCycle.read_jsons(directory,name)
	close_enough = IntrinsicDACCycle.Close_enough(material, Kh_N₂, One_atm_N₂)
	saturation_any_co2 = saturation_adsorb_co2(directory, name)
	continue_test = close_enough & saturation_any_co2
	Results_Dict["Close Enough to Linear Adsorption"] = close_enough
	Results_Dict["Saturation Adsorb Any CO2"] = saturation_any_co2
	if continue_test == false
		return Results_Dict
	end
	
	
	#Define the limits of the parameters	
    #specify Start and total Step
	T_start = 273.0 #[K] 
	ΔT = 200.0 #[K]
    dT_max = 0.25

	T_start_lower = 200.0 #[K]
	T_start_upper = 400.0 #[K]

	ΔT_lower = 0.0 #[K]
	ΔT_upper = 200.0 #[K]

	P_start = 101325.0 #[Pa]
	#limit the ΔP to only reach "rough vaccuum" 100 Pa.
	ΔP = 100.0 - P_start #[Pa] the offset of 1 Pa ensures that the path of P never gets to 0 Pa. 
    dP_max = -125.0

	P_start_lower = 101325.0 #[Pa] 
	P_start_upper = 1.1 .* [101325.0] #[Pa]

	#Lower limit of ΔP is to reach "rough vaccum" 100 Pa
	ΔP_lower = 100.0 - P_start #[Pa]
	# ΔP_lower = (500 .- P_start[1])./steps .+ zero(ΔP) #[Pa]
	ΔP_upper = 0.0 #[Pa]

    #Combine the bounds into one object for pasing to optimizer
    lower_bound = cat(T_start_lower, ΔT_lower, P_start_lower, ΔP_lower, dims = 1)
	upper_bound = cat(T_start_upper, ΔT_upper, P_start_upper, ΔP_upper, dims = 1)
	bounds = cat(lower_bound, upper_bound, dims = 2)'



    #Perform the multi objective optimization
    N= 800
	method = Metaheuristics.NSGA3(N= N)

    function f(x)
        return ScorePath(x, dT_max, dP_max, Base_directory, name, α, "objectives", 0)
    end

	Metaheuristics.optimize!(f, bounds, method)

    #Unpack the results (Pareto front of perfrormance metrics and their parameters)
    results_state = Metaheuristics.get_result(method)
	results = pareto_front(results_state)
    
	non_dom_solutions = Metaheuristics.get_non_dominated_solutions(method.status.population)
	num_solutions = length(non_dom_solutions)

	results_parameters = Vector{Vector}(undef, num_solutions)

	for (i,solution) in enumerate(non_dom_solutions)
		
		results_parameters[i] = solution.x
	end
	results_parameters = cat(results_parameters..., dims = 2)'
    
    path_Ts = Vector{Vector}(undef, length(results_parameters[:,1]))
	path_Ps = Vector{Vector}(undef, length(results_parameters[:,1]))
	for (i,param) in enumerate(eachrow(results_parameters))
        
        #To recover the stepsize run ScorePath again
        objectives = ScorePath(param, dT_max, dP_max, Base_directory, name, α, "steps")
        dT_max = objectives[4]
        dP_max = objectives[5]

		path_T_start = param[1]
		path_ΔT = param[2]
		path_P_start = param[3]
		path_ΔP = param[4]
	
		path_T_steps = length(path_T_start:dT_max:path_T_start+path_ΔT)
		path_P_steps = length(path_P_start:-dP_max:path_P_start+path_ΔP)
	
		path_steps = maximum([path_T_steps, path_P_steps])
	
		path_Ts[i] = collect(LinRange(path_T_start, path_T_start+path_ΔT, path_steps))
		path_Ps[i] = collect(LinRange(path_P_start, path_P_start+path_ΔP, path_steps))
	end

    #Initialize a dictionary to store all the intermediates
    Path_Dict = sort(Dict{String, Any}("Refresh_Path" => "Definition of refresh path and material properties along that path (in per kg of sorbent basis)."))
    E_Balance_Dict = sort(Dict{String, Any}("E_Balance" => "Energy balance along path"))
    Step_1_Dict = sort(Dict{String, Any}("Step_1" => "Adsorption"))
    Step_2_Dict = sort(Dict{String, Any}("Step_2" => "Desorption"))
    Step_3_Dict = sort(Dict{String, Any}("Step_3" => "Waste energy recovery"))

	Path_Dict["Temperatures"] = Vector{Vector}(undef,num_solutions)
    Path_Dict["Temperature_units"] = "K"
    Path_Dict["Pressures"] = Vector{Vector}(undef,num_solutions)
    Path_Dict["Pressure_units"] = "Pa"
    Path_Dict["Betas"] = Vector{Vector}(undef,num_solutions)
    Path_Dict["Beta_units"] = "mol/kJ"

	Path_Dict["Henry_CO2"] = Vector{Vector}(undef,num_solutions)
    Path_Dict["Henry_CO2_err"] = Vector{Vector}(undef,num_solutions)
    Path_Dict["Henry_N2"] = Vector{Vector}(undef,num_solutions)
    Path_Dict["Henry_N2_err"] = Vector{Vector}(undef,num_solutions)
    Path_Dict["Henry_units"] = "mmol/(kg Pa)"

    Path_Dict["Moles_CO2"] = Vector{Vector}(undef,num_solutions)
    Path_Dict["Moles_N2"] = Vector{Vector}(undef,num_solutions)
    Path_Dict["Moles_units"] = "mol/kg"

    Path_Dict["Heat_of_adsorb_CO2"] = Vector{Matrix}(undef,num_solutions)
    Path_Dict["Heat_of_adsorb_CO2_err"] = Vector{Matrix}(undef,num_solutions)
    Path_Dict["Heat_of_adsorb_N2"] = Vector{Matrix}(undef,num_solutions)
    Path_Dict["Heat_of_adsorb_N2_err"] = Vector{Matrix}(undef,num_solutions)
    Path_Dict["Heat_of_adsorb_units"] = "J/mol"

    Path_Dict["Specific_heat_sorbent"] = Vector{Matrix}(undef,num_solutions)
    Path_Dict["Specific_heat_sorbent_err"] = Vector{Matrix}(undef,num_solutions)
    Path_Dict["Specific_heat_sorbent_units"] = "J/(kg K)"

	Path_Dict["Flag"] = Vector{String}(undef,num_solutions)

    Step_1_Dict["Heat_to_adsorb_CO2"] = Vector{Real}(undef,num_solutions)
    Step_1_Dict["Heat_to_adsorb_N2"] = Vector{Real}(undef,num_solutions)
    Step_1_Dict["Work_to_adsorb_CO2"] = Vector{Real}(undef,num_solutions)
    Step_1_Dict["Work_to_adsorb_N2"] = Vector{Real}(undef,num_solutions)
    Step_1_Dict["E_units"] = "J/kg_sorb"

    E_Balance_Dict["E1"] = Vector{Real}(undef,num_solutions)

    Step_2_Dict["Heat_to_desorb_CO2"] = Vector{Matrix}(undef,num_solutions)
    Step_2_Dict["Heat_to_desorb_N2"] = Vector{Matrix}(undef,num_solutions)
    Step_2_Dict["Work_to_desorb_CO2"] = Vector{Vector}(undef,num_solutions)
    Step_2_Dict["Work_to_desorb_N2"] = Vector{Vector}(undef,num_solutions)
    Step_2_Dict["E_to_heat_adsorbed_CO2"] = Vector{Vector}(undef,num_solutions)
    Step_2_Dict["E_to_heat_adsorbed_N2"] = Vector{Vector}(undef,num_solutions)
    Step_2_Dict["E_to_heat_sorbent"] = Vector{Matrix}(undef,num_solutions)
    Step_2_Dict["E_to_change_pressure"] = Vector{Vector}(undef,num_solutions)
    Step_2_Dict["E_units"] = "J/kg_sorb"

    E_Balance_Dict["E2"] = Vector{Real}(undef,num_solutions)

    Step_3_Dict["E_recovered"] = Vector{Real}(undef,num_solutions)
    E_Balance_Dict["E3"] = Vector{Real}(undef,num_solutions)

    E_Balance_Dict["Total_E_of_cycle"] = Vector{Real}(undef,num_solutions)
    E_Balance_Dict["E_units"] = "J/kg_sorb"

    Results_Dict["Captured_CO2"] = Vector{Real}(undef,num_solutions)
    Results_Dict["Captured_N2"] = Vector{Real}(undef,num_solutions)
    Results_Dict["Captured_gas_units"] = "mol/kg_sorb"

    Results_Dict["Intrinsic_capture_efficiency"] = Vector{Real}(undef,num_solutions)
    Results_Dict["Intrinsic_capture_efficiency_units"] = "mol/J"
    Results_Dict["Purity_captured_CO2"] = Vector{Real}(undef,num_solutions)

    #Re-evalutate the path at just the pareto front parameters
    for i in 1:num_solutions
        #Re-evalutate the path
        ith_results = IntrinsicDACCycle.Intrinisic_refresh_path(Base_directory, name,
                                                path_Ts[i], path_Ps[i], α)

        #Re-package resutls
		if "Refresh_Path" in keys(ith_results)
	        ith_path_dict = ith_results["Refresh_Path"]
	        ith_e_balance_dict = ith_results["E_Balance"]
	        ith_step_1_dict = ith_e_balance_dict["Step_1"]
	        ith_step_2_dict = ith_e_balance_dict["Step_2"]
	        ith_step_3_dict = ith_e_balance_dict["Step_3"]
	
	        Path_Dict["Temperatures"][i] = ith_path_dict["Temperatures"]
	        Path_Dict["Pressures"][i] = ith_path_dict["Pressures"]
	        Path_Dict["Betas"][i] = ith_path_dict["Betas"]
	        
	        Path_Dict["Henry_CO2"][i] = ith_path_dict["Henry_CO2"]
	        Path_Dict["Henry_CO2_err"][i] = ith_path_dict["Henry_CO2_err"]
	        Path_Dict["Henry_N2"][i] = ith_path_dict["Henry_N2"]
	        Path_Dict["Henry_N2_err"][i] = ith_path_dict["Henry_N2_err"]
	        
	        Path_Dict["Moles_CO2"][i] = ith_path_dict["Moles_CO2"]
	        Path_Dict["Moles_N2"][i] = ith_path_dict["Moles_N2"]
	        
	        Path_Dict["Heat_of_adsorb_CO2"][i] = ith_path_dict["Heat_of_adsorb_CO2"]
	        Path_Dict["Heat_of_adsorb_CO2_err"][i] = ith_path_dict["Heat_of_adsorb_CO2_err"]
	        Path_Dict["Heat_of_adsorb_N2"][i] = ith_path_dict["Heat_of_adsorb_N2"]
	        Path_Dict["Heat_of_adsorb_N2_err"][i] = ith_path_dict["Heat_of_adsorb_N2_err"]
	        
	        Path_Dict["Specific_heat_sorbent"][i] = ith_path_dict["Specific_heat_sorbent"]
	        Path_Dict["Specific_heat_sorbent_err"][i] = ith_path_dict["Specific_heat_sorbent_err"]
	        
	        Step_1_Dict["Heat_to_adsorb_CO2"][i] = ith_step_1_dict["Heat_to_adsorb_CO2"]
	        Step_1_Dict["Heat_to_adsorb_N2"][i] = ith_step_1_dict["Heat_to_adsorb_N2"]
	        Step_1_Dict["Work_to_adsorb_CO2"][i] = ith_step_1_dict["Work_to_adsorb_CO2"]
	        Step_1_Dict["Work_to_adsorb_N2"][i] = ith_step_1_dict["Work_to_adsorb_N2"]
	
	        E_Balance_Dict["E1"][i] = ith_e_balance_dict["E1"]
	        
	        Step_2_Dict["Heat_to_desorb_CO2"][i] = ith_step_2_dict["Heat_to_desorb_CO2"]
	        Step_2_Dict["Heat_to_desorb_N2"][i] = ith_step_2_dict["Heat_to_desorb_N2"]
	        Step_2_Dict["Work_to_desorb_CO2"][i] = ith_step_2_dict["Work_to_desorb_CO2"]
	        Step_2_Dict["Work_to_desorb_N2"][i] = ith_step_2_dict["Work_to_desorb_N2"]
	        Step_2_Dict["E_to_heat_adsorbed_CO2"][i] = ith_step_2_dict["E_to_heat_adsorbed_CO2"]
	        Step_2_Dict["E_to_heat_adsorbed_N2"][i] = ith_step_2_dict["E_to_heat_adsorbed_N2"]
	        Step_2_Dict["E_to_heat_sorbent"][i] = ith_step_2_dict["E_to_heat_sorbent"]
	        Step_2_Dict["E_to_change_pressure"][i] = ith_step_2_dict["E_to_change_pressure"]
	        
	        E_Balance_Dict["E2"][i] = ith_e_balance_dict["E2"]
	    
	        Step_3_Dict["E_recovered"][i] = ith_step_3_dict["E_recovered"]
	        E_Balance_Dict["E3"][i] = ith_e_balance_dict["E3"]
	    
	        E_Balance_Dict["Total_E_of_cycle"][i] = ith_e_balance_dict["Total_E_of_cycle"]
	        
	        Results_Dict["Captured_CO2"][i] = ith_results["Captured_CO2"]
	        Results_Dict["Captured_N2"][i] = ith_results["Captured_N2"]
	        
	        Results_Dict["Intrinsic_capture_efficiency"][i] = ith_results["Intrinsic_capture_efficiency"]
	        Results_Dict["Purity_captured_CO2"][i] = ith_results["Purity_captured_CO2"]
		else
			Path_Dict["Flag"][i] = "No path, check Henry constant or saturation calculations"
		end
    end

    #Package all the dictionaries together
    E_Balance_Dict["Step_1"] = Step_1_Dict
    E_Balance_Dict["Step_2"] = Step_2_Dict
    E_Balance_Dict["Step_3"] = Step_3_Dict

    Results_Dict["Refresh_Path"] = Path_Dict
    Results_Dict["E_Balance"] = E_Balance_Dict

    return Results_Dict

end








"""Function to optimize the start and end temperatures and pressures 
for a given material and inlet CO2 concentration.
Assumes a linear path,
With the uncertainties in performance metrics"""
function Optimize_Intrinsic_Refresh_w_err(Base_directory::String, name::String,  
                                    α::Real)

	#Initialize a dictionary to store all the results
	Results_Dict = sort(Dict{String, Any}("Name" => name))
    
	#do the close_enough to linear adsorption test
	#and the saturation test before attempting the calculation
	material, Kh_N₂, Kh_CO₂, One_atm_N₂ = IntrinsicDACCycle.read_jsons(directory,name)
	close_enough = IntrinsicDACCycle.Close_enough(material, Kh_N₂, One_atm_N₂)
	saturation_any_co2 = saturation_adsorb_co2(directory, name)
	continue_test = close_enough & saturation_any_co2
	Results_Dict["Close Enough to Linear Adsorption"] = close_enough
	Results_Dict["Saturation Adsorb Any CO2"] = saturation_any_co2
	if continue_test == false
		return Results_Dict
	end
	
    #Define the limits of the parameters	
    #specify Start and total Step
	T_start = 273.0 #[K] 
	ΔT = 200.0 #[K]
    dT_max = 0.25

	T_start_lower = 200.0 #[K]
	T_start_upper = 400.0 #[K]

	ΔT_lower = 0.0 #[K]
	ΔT_upper = 200.0 #[K]

	P_start = 101325.0 #[Pa]
	#limit the ΔP to only reach "rough vaccuum" 100 Pa.
	ΔP = 100.0 - P_start #[Pa] the offset of 1 Pa ensures that the path of P never gets to 0 Pa. 
    dP_max = -125.0

	P_start_lower = 101325.0 #[Pa] 
	P_start_upper = 1.1 .* [101325.0] #[Pa]

	#Lower limit of ΔP is to reach "rough vaccum" 100 Pa
	ΔP_lower = 100.0 - P_start #[Pa]
	# ΔP_lower = (500 .- P_start[1])./steps .+ zero(ΔP) #[Pa]
	ΔP_upper = 0.0 #[Pa]

    #Combine the bounds into one object for pasing to optimizer
    lower_bound = cat(T_start_lower, ΔT_lower, P_start_lower, ΔP_lower, dims = 1)
	upper_bound = cat(T_start_upper, ΔT_upper, P_start_upper, ΔP_upper, dims = 1)
	bounds = cat(lower_bound, upper_bound, dims = 2)'

    #Perform the multi objective optimization
    N= 800
	method = Metaheuristics.NSGA3(N= N)
    function f(x)
        return ScorePath(x, dT_max, dP_max, Base_directory, name, α, "objectives", 0)
    end

	Metaheuristics.optimize!(f, bounds, method)

    #Unpack the results (Pareto front of perfrormance metrics and their parameters)
    results_state = Metaheuristics.get_result(method)
	results = pareto_front(results_state)
    
	non_dom_solutions = Metaheuristics.get_non_dominated_solutions(method.status.population)
	num_solutions = length(non_dom_solutions)

	results_parameters = Vector{Vector}(undef, num_solutions)

	for (i,solution) in enumerate(non_dom_solutions)
		
		results_parameters[i] = solution.x
	end
	results_parameters = cat(results_parameters..., dims = 2)'
    
    path_Ts = Vector{Vector}(undef, length(results_parameters[:,1]))
	path_Ps = Vector{Vector}(undef, length(results_parameters[:,1]))
	for (i,param) in enumerate(eachrow(results_parameters))

        #To recover the stepsize run ScorePath again
        objectives = ScorePath(param, dT_max, dP_max, Base_directory, name, α, "steps")
        dT_max = objectives[4]
        dP_max = objectives[5]
	
		path_T_start = param[1]
		path_ΔT = param[2]
		path_P_start = param[3]
		path_ΔP = param[4]
	
		path_T_steps = length(path_T_start:dT_max:path_T_start+path_ΔT)
		path_P_steps = length(path_P_start:-dP_max:path_P_start+path_ΔP)
	
		path_steps = maximum([path_T_steps, path_P_steps])
	
		path_Ts[i] = collect(LinRange(path_T_start, path_T_start+path_ΔT, path_steps))
		path_Ps[i] = collect(LinRange(path_P_start, path_P_start+path_ΔP, path_steps))
	end

    #Initialize a dictionary to store all the intermediates
    Path_Dict = sort(Dict{String, Any}("Refresh_Path" => "Definition of refresh path and material properties along that path (in per kg of sorbent basis)."))
    E_Balance_Dict = sort(Dict{String, Any}("E_Balance" => "Energy balance along path"))
    Step_1_Dict = sort(Dict{String, Any}("Step_1" => "Adsorption"))
    Step_2_Dict = sort(Dict{String, Any}("Step_2" => "Desorption"))
    Step_3_Dict = sort(Dict{String, Any}("Step_3" => "Waste energy recovery"))

	Path_Dict["Temperatures"] = Vector{Vector}(undef,num_solutions)
    Path_Dict["Temperature_units"] = "K"
    Path_Dict["Pressures"] = Vector{Vector}(undef,num_solutions)
    Path_Dict["Pressure_units"] = "Pa"
    Path_Dict["Betas"] = Vector{Vector}(undef,num_solutions)
    Path_Dict["Beta_units"] = "mol/kJ"

	Path_Dict["Henry_CO2"] = Vector{Vector}(undef,num_solutions)
    Path_Dict["Henry_CO2_err"] = Vector{Vector}(undef,num_solutions)
    Path_Dict["Henry_N2"] = Vector{Vector}(undef,num_solutions)
    Path_Dict["Henry_N2_err"] = Vector{Vector}(undef,num_solutions)
    Path_Dict["Henry_units"] = "mmol/(kg Pa)"

    Path_Dict["Moles_CO2"] = Vector{Vector}(undef,num_solutions)
    Path_Dict["Moles_N2"] = Vector{Vector}(undef,num_solutions)
    Path_Dict["Moles_units"] = "mol/kg"

    Path_Dict["Heat_of_adsorb_CO2"] = Vector{Matrix}(undef,num_solutions)
    Path_Dict["Heat_of_adsorb_CO2_err"] = Vector{Matrix}(undef,num_solutions)
    Path_Dict["Heat_of_adsorb_N2"] = Vector{Matrix}(undef,num_solutions)
    Path_Dict["Heat_of_adsorb_N2_err"] = Vector{Matrix}(undef,num_solutions)
    Path_Dict["Heat_of_adsorb_units"] = "J/mol"

    Path_Dict["Specific_heat_sorbent"] = Vector{Matrix}(undef,num_solutions)
    Path_Dict["Specific_heat_sorbent_err"] = Vector{Matrix}(undef,num_solutions)
    Path_Dict["Specific_heat_sorbent_units"] = "J/(kg K)"

	Path_Dict["Flag"] = Vector{String}(undef,num_solutions)

    Step_1_Dict["Heat_to_adsorb_CO2"] = Vector{Real}(undef,num_solutions)
    Step_1_Dict["Heat_to_adsorb_N2"] = Vector{Real}(undef,num_solutions)
    Step_1_Dict["Work_to_adsorb_CO2"] = Vector{Real}(undef,num_solutions)
    Step_1_Dict["Work_to_adsorb_N2"] = Vector{Real}(undef,num_solutions)
    Step_1_Dict["E_units"] = "J/kg_sorb"

    E_Balance_Dict["E1"] = Vector{Real}(undef,num_solutions)

    Step_2_Dict["Heat_to_desorb_CO2"] = Vector{Matrix}(undef,num_solutions)
    Step_2_Dict["Heat_to_desorb_N2"] = Vector{Matrix}(undef,num_solutions)
    Step_2_Dict["Work_to_desorb_CO2"] = Vector{Vector}(undef,num_solutions)
    Step_2_Dict["Work_to_desorb_N2"] = Vector{Vector}(undef,num_solutions)
    Step_2_Dict["E_to_heat_adsorbed_CO2"] = Vector{Vector}(undef,num_solutions)
    Step_2_Dict["E_to_heat_adsorbed_N2"] = Vector{Vector}(undef,num_solutions)
    Step_2_Dict["E_to_heat_sorbent"] = Vector{Matrix}(undef,num_solutions)
    Step_2_Dict["E_to_change_pressure"] = Vector{Vector}(undef,num_solutions)
    Step_2_Dict["E_units"] = "J/kg_sorb"

    E_Balance_Dict["E2"] = Vector{Real}(undef,num_solutions)

    Step_3_Dict["E_recovered"] = Vector{Real}(undef,num_solutions)
    E_Balance_Dict["E3"] = Vector{Real}(undef,num_solutions)

    E_Balance_Dict["Total_E_of_cycle"] = Vector{Real}(undef,num_solutions)
    E_Balance_Dict["E_units"] = "J/kg_sorb"

    Results_Dict["Captured_CO2"] = Vector{Real}(undef,num_solutions)
    Results_Dict["Captured_CO2_std"] = Vector{Real}(undef,num_solutions)
    Results_Dict["Captured_N2"] = Vector{Real}(undef,num_solutions)
    Results_Dict["Captured_N2_std"] = Vector{Real}(undef,num_solutions)
    Results_Dict["Captured_gas_units"] = "mol/kg_sorb"

    Results_Dict["Intrinsic_capture_efficiency"] = Vector{Real}(undef,num_solutions)
    Results_Dict["Intrinsic_capture_efficiency_std"] = Vector{Real}(undef,num_solutions)
    Results_Dict["Intrinsic_capture_efficiency_units"] = "mol/J"
    Results_Dict["Purity_captured_CO2"] = Vector{Real}(undef,num_solutions)
    Results_Dict["Purity_captured_CO2_std"] = Vector{Real}(undef,num_solutions)

    #Re-evalutate the path at just the pareto front parameters
    for i in 1:num_solutions
        #Re-evalutate the path
        # if there is a sampling that give a domain error, re-try with smaller step sizes a few times

        trial_Ts = path_Ts[i]
        trial_Ps = path_Ps[i]

        #Initialize containers for the results
        ith_results = nothing
        objectives_dist = nothing

        for attempts in 1:10
            ith_results = IntrinsicDACCycle.Intrinisic_refresh_path(Base_directory, name,
                                                    trial_Ts, trial_Ps, α)

            objectives_dist = IntrinsicDACCycle.Intrinisic_refresh_objectives_posterior_dist(Base_directory, name,
                                                                                    trial_Ts, trial_Ps, α,
                                                                                    100)
			if "Refresh_Path" in keys(ith_results)
	            #Test if Ts_s are not NaNs (i.e. failed Close_enough test earlier)
	            test_1 = ~any(isnan.(path_Ts[i]))
	            #Test for domain error
	            test_2 = any(isnan.(objectives_dist[1]))
	            #if it passed the Close_enough test but still created a domain error 
	            if test_1 & test_2
	                new_Ts_start = path_Ts[i][1]
	                new_Ts_end = path_Ts[i][end]
	                #Create new T steps at finer resultion
	                new_dT = (path_Ts[i][2] - path_Ts[i][1])*0.75^attempts
	
	                new_Ps_start = path_Ps[i][1]
	                new_Ps_end = path_Ps[i][end]
	                #Create new P steps at finer resoluion
	                new_dP = (path_Ps[i][2] - path_Ps[i][1]) *0.75^attempts
	
	                new_path_T_steps = length(new_Ts_start:new_dT:new_Ts_end)
	                new_path_P_steps = length(new_Ps_start:new_dP:new_Ps_end)
	
	                new_path_steps = maximum([new_path_T_steps, new_path_P_steps])
	
	            
	                #Re-define the path Ts and Ps
	                trial_Ts = collect(LinRange(new_Ts_start, new_Ts_end, new_path_steps))
	                trial_Ps = collect(LinRange(new_Ps_start, new_Ps_end, new_path_steps))
	                print("Uncertainty calculation lead to numerical instabilty, re-trying at finer step size.")
	            else
	                break
	            end
			end
        end

        if "Refresh_Path" in keys(ith_results)
	        mean_ξ = mean(objectives_dist[1])
	        std_ξ = std(objectives_dist[1])
	        mean_α = mean(objectives_dist[2])
	        std_α = std(objectives_dist[2])
	        
	        mean_Δn_CO2 = mean(objectives_dist[3])
	        std_Δn_CO2 = std(objectives_dist[3])
	        mean_Δn_N2 = mean(objectives_dist[4])
	        std_Δn_N2 = std(objectives_dist[4])
	        

        #Re-package resutls
		
	        ith_path_dict = ith_results["Refresh_Path"]
	        ith_e_balance_dict = ith_results["E_Balance"]
	        ith_step_1_dict = ith_e_balance_dict["Step_1"]
	        ith_step_2_dict = ith_e_balance_dict["Step_2"]
	        ith_step_3_dict = ith_e_balance_dict["Step_3"]
	
	        Path_Dict["Temperatures"][i] = ith_path_dict["Temperatures"]
	        Path_Dict["Pressures"][i] = ith_path_dict["Pressures"]
	        Path_Dict["Betas"][i] = ith_path_dict["Betas"]
	        
	        Path_Dict["Henry_CO2"][i] = ith_path_dict["Henry_CO2"]
	        Path_Dict["Henry_CO2_err"][i] = ith_path_dict["Henry_CO2_err"]
	        Path_Dict["Henry_N2"][i] = ith_path_dict["Henry_N2"]
	        Path_Dict["Henry_N2_err"][i] = ith_path_dict["Henry_N2_err"]
	        
	        Path_Dict["Moles_CO2"][i] = ith_path_dict["Moles_CO2"]
	        Path_Dict["Moles_N2"][i] = ith_path_dict["Moles_N2"]
	        
	        Path_Dict["Heat_of_adsorb_CO2"][i] = ith_path_dict["Heat_of_adsorb_CO2"]
	        Path_Dict["Heat_of_adsorb_CO2_err"][i] = ith_path_dict["Heat_of_adsorb_CO2_err"]
	        Path_Dict["Heat_of_adsorb_N2"][i] = ith_path_dict["Heat_of_adsorb_N2"]
	        Path_Dict["Heat_of_adsorb_N2_err"][i] = ith_path_dict["Heat_of_adsorb_N2_err"]
	        
	        Path_Dict["Specific_heat_sorbent"][i] = ith_path_dict["Specific_heat_sorbent"]
	        Path_Dict["Specific_heat_sorbent_err"][i] = ith_path_dict["Specific_heat_sorbent_err"]
	        
	        Step_1_Dict["Heat_to_adsorb_CO2"][i] = ith_step_1_dict["Heat_to_adsorb_CO2"]
	        Step_1_Dict["Heat_to_adsorb_N2"][i] = ith_step_1_dict["Heat_to_adsorb_N2"]
	        Step_1_Dict["Work_to_adsorb_CO2"][i] = ith_step_1_dict["Work_to_adsorb_CO2"]
	        Step_1_Dict["Work_to_adsorb_N2"][i] = ith_step_1_dict["Work_to_adsorb_N2"]
	
	        E_Balance_Dict["E1"][i] = ith_e_balance_dict["E1"]
	        
	        Step_2_Dict["Heat_to_desorb_CO2"][i] = ith_step_2_dict["Heat_to_desorb_CO2"]
	        Step_2_Dict["Heat_to_desorb_N2"][i] = ith_step_2_dict["Heat_to_desorb_N2"]
	        Step_2_Dict["Work_to_desorb_CO2"][i] = ith_step_2_dict["Work_to_desorb_CO2"]
	        Step_2_Dict["Work_to_desorb_N2"][i] = ith_step_2_dict["Work_to_desorb_N2"]
	        Step_2_Dict["E_to_heat_adsorbed_CO2"][i] = ith_step_2_dict["E_to_heat_adsorbed_CO2"]
	        Step_2_Dict["E_to_heat_adsorbed_N2"][i] = ith_step_2_dict["E_to_heat_adsorbed_N2"]
	        Step_2_Dict["E_to_heat_sorbent"][i] = ith_step_2_dict["E_to_heat_sorbent"]
	        Step_2_Dict["E_to_change_pressure"][i] = ith_step_2_dict["E_to_change_pressure"]
	        
	        E_Balance_Dict["E2"][i] = ith_e_balance_dict["E2"]
	    
	        Step_3_Dict["E_recovered"][i] = ith_step_3_dict["E_recovered"]
	        E_Balance_Dict["E3"][i] = ith_e_balance_dict["E3"]
	    
	        E_Balance_Dict["Total_E_of_cycle"][i] = ith_e_balance_dict["Total_E_of_cycle"]
	        
	        Results_Dict["Captured_CO2"][i] = mean_Δn_CO2
	        Results_Dict["Captured_CO2_std"][i] = std_Δn_CO2
	        Results_Dict["Captured_N2"][i] = mean_Δn_N2
	        Results_Dict["Captured_N2_std"][i] = std_Δn_N2
	        
	        Results_Dict["Intrinsic_capture_efficiency"][i] = mean_ξ
	        Results_Dict["Intrinsic_capture_efficiency_std"][i] = std_ξ
	        Results_Dict["Purity_captured_CO2"][i] = mean_α
	        Results_Dict["Purity_captured_CO2_std"][i] = std_α
		else
			Path_Dict["Flag"][i] = "No path, check Henry constant or saturation calculations"
		end
				
    end

    #Package all the dictionaries together
    E_Balance_Dict["Step_1"] = Step_1_Dict
    E_Balance_Dict["Step_2"] = Step_2_Dict
    E_Balance_Dict["Step_3"] = Step_3_Dict

    Results_Dict["Refresh_Path"] = Path_Dict
    Results_Dict["E_Balance"] = E_Balance_Dict

    return Results_Dict

end




"""Function to optimize the start and end temperatures and pressures 
for a given material and inlet CO2 concentration.
Assumes a linear path,
Reports the path parameters sampled and thier performances"""
function Optimize_Intrinsic_Refresh_path_distributions(Base_directory::String, name::String,  
                                    α::Real)

    #Define the limits of the parameters	
    #specify Start and total Step
	T_start = 273.0 #[K] 
	ΔT = 200.0 #[K]
    dT_max = 0.25

	T_start_lower = 200.0 #[K]
	T_start_upper = 400.0 #[K]

	ΔT_lower = 0.0 #[K]
	ΔT_upper = 200.0 #[K]

	P_start = 101325.0 #[Pa]
	#limit the ΔP to only reach "rough vaccuum" 100 Pa.
	ΔP = 100.0 - P_start #[Pa] the offset of 1 Pa ensures that the path of P never gets to 0 Pa. 
    dP_max = -125.0

	P_start_lower = 101325.0 #[Pa] 
	P_start_upper = 1.1 .* [101325.0] #[Pa]

	#Lower limit of ΔP is to reach "rough vaccum" 100 Pa
	ΔP_lower = 100.0 - P_start #[Pa]
	# ΔP_lower = (500 .- P_start[1])./steps .+ zero(ΔP) #[Pa]
	ΔP_upper = 0.0 #[Pa]

    #Combine the bounds into one object for pasing to optimizer
    lower_bound = cat(T_start_lower, ΔT_lower, P_start_lower, ΔP_lower, dims = 1)
	upper_bound = cat(T_start_upper, ΔT_upper, P_start_upper, ΔP_upper, dims = 1)
	bounds = cat(lower_bound, upper_bound, dims = 2)'


    #Perform the multi objective optimization
    N= 800
	method = Metaheuristics.NSGA3(N= N)
    function f(x)
        return ScorePath(x, dT_max, dP_max, Base_directory, name, α, "objectives", 0)
    end
	Metaheuristics.optimize!(f, bounds, method)

    #Unpack the results (Pareto front of perfrormance metrics and their parameters)
    results_state = Metaheuristics.get_result(method)
	results = pareto_front(results_state)
    

    return method.status.population
end




end

# ╔═╡ 68539c35-deef-4187-b2f6-00d47ca65b29
moreofatest = Optimize_Intrinsic_Refresh(directory, name, α)

# ╔═╡ 4e02108e-bc90-473f-ac0d-faea6014c5f7
begin
	dT_max = 0.25
	dP_max = -125.0
	parameters = [199, 20, 101325.0, 0.0]
	
	result = ScorePath(parameters, dT_max, dP_max, directory, name, α)
end 

# ╔═╡ fd131f2e-b3c8-4581-b436-1fb45b867193
begin
	a = missing

	if ismissing(a)
		@show "Hi"
	end
	b = 10 
end

# ╔═╡ 9d80a406-f462-4ce5-99ea-c10f0bde60b9
typeof(a + b)

# ╔═╡ 9ac5e77d-e413-4e9b-8626-75b51f977b56
moreofatest2 = Optimize_Intrinsic_Refresh(directory, name, α)

# ╔═╡ a19b4eaa-73a5-4580-8d74-aa55afaa8237


# ╔═╡ c7924bae-b66a-4e0a-9953-b9bff9e019f5


# ╔═╡ 33d813a7-2b5d-4cab-be26-2a61cbd71ea4
size(sat_test[1])[1] == 0 

# ╔═╡ 7ba46d2f-2a76-410d-af67-6524512ab0ed
size(sat_test[1])

# ╔═╡ 5e12128a-6b5c-4185-8cb2-7b30c7b65fb9
"terminated_early" in keys(saturation_dict["isotherm"])

# ╔═╡ 2202ad42-8619-409e-bdd5-804ceae70200
function foobar(x)
	if x > 5
		return "a thing"

	end
	@show x
	return x + 2
end

# ╔═╡ b0ae97ea-c78e-454a-a755-183cc280c2b1
foobar(6)

# ╔═╡ Cell order:
# ╠═c6ccef8c-4db6-11ef-0f54-0f3bae2f93b4
# ╠═6c576153-4bca-411a-b745-3a3cb88d6ad9
# ╠═a7a44fe8-29c9-4b65-93f0-ae2672eb1960
# ╠═7927d0de-3ca7-408d-bdb2-93d252af9a80
# ╠═5ee89c88-5efe-40e2-bf24-af8facfa2199
# ╠═d31a012a-515c-4521-9e72-600eb36b6083
# ╠═2012b243-b921-4c28-a31c-74969c46276f
# ╠═316b64c9-029d-456b-b1b9-ec82a14a7f1b
# ╠═d177558c-ff0e-4089-a79b-11c27007c759
# ╠═3c8fa7bc-b4c2-41cd-b050-a98e3407b43e
# ╠═ce51f824-77fa-402a-a16e-372d5f61838b
# ╠═f6e3c4b8-b4dd-463e-8246-c49305237a31
# ╠═edea02bf-4658-40b4-b11f-2faabf1a7200
# ╠═1420ac5b-c8ba-4ccd-a047-9a21ba6b4b97
# ╠═ea1ccc4d-4a86-4d53-b865-1c7986d3598a
# ╠═80cd44e7-af7c-4751-9e78-77d94a82f32b
# ╠═1d569789-5ddb-4e14-a063-4103cfb75de1
# ╠═62606add-5fc3-4518-a88f-4bd0b253ce39
# ╠═cbbe86a4-06cf-4344-a2db-a90980ccf481
# ╠═7b17395c-7cc2-4dc0-a361-17114e08bbec
# ╠═1e37e81f-6d2f-4b9a-815f-6f8699ab7083
# ╠═00cd9198-483e-4d7f-9add-11596183583b
# ╠═63b7a1fe-b67a-4d15-a5b5-486723daaaf5
# ╠═63d7bd63-181a-4a99-84b6-2931a468726b
# ╠═159eee30-a3c4-4a9d-848e-3f7c2f75f7cf
# ╠═69c4ed3f-dd56-4a3a-9740-704a10ce07f5
# ╠═fe759187-2198-48d4-be35-309608b16451
# ╠═78a2b7b2-3299-4ab5-b99b-97d37ac23d9c
# ╠═fe85ba0d-68de-46ba-a6de-434bf81af8b4
# ╠═f3c54cae-0ce5-48e3-883a-b64da8818050
# ╠═c7004efc-4d62-4ab8-a787-e25f9e84c706
# ╠═d3e9d752-dcff-4486-93dd-82d68a6be81e
# ╠═49169e8a-868d-42d7-9fc7-5c7b9eaea90c
# ╠═438c46d9-5872-4c58-9902-de159e90f16a
# ╠═b6f58e4f-2072-4d55-9be1-be4f357f74fa
# ╠═68539c35-deef-4187-b2f6-00d47ca65b29
# ╠═43853edf-20bd-4bb6-bdf0-b60a0aee225e
# ╠═59b733d3-e361-423a-a010-966fe16c8176
# ╠═a8700a52-51a2-49d6-8c22-cf331a54731e
# ╠═b0405caf-1431-4c04-bb7c-56fe9cc468c8
# ╠═4e02108e-bc90-473f-ac0d-faea6014c5f7
# ╠═fd131f2e-b3c8-4581-b436-1fb45b867193
# ╠═9d80a406-f462-4ce5-99ea-c10f0bde60b9
# ╠═9ac5e77d-e413-4e9b-8626-75b51f977b56
# ╠═a19b4eaa-73a5-4580-8d74-aa55afaa8237
# ╠═c7924bae-b66a-4e0a-9953-b9bff9e019f5
# ╠═33d813a7-2b5d-4cab-be26-2a61cbd71ea4
# ╠═7ba46d2f-2a76-410d-af67-6524512ab0ed
# ╠═5e12128a-6b5c-4185-8cb2-7b30c7b65fb9
# ╠═2202ad42-8619-409e-bdd5-804ceae70200
# ╠═b0ae97ea-c78e-454a-a755-183cc280c2b1
