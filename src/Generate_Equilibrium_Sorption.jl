

"""Set up the Python libraries."""
# pyiast = PyCall.pyimport("pyiast")
# pd = PyCall.pyimport("pandas")

"""Function to solve for the concentration of CO2 in the outlet.
This assumes that the change in moles of gas in the absolute adsorption is what is sent to the outlet. 
In this case the CO2 concentration in the outlet should be the same as the CO2 concentration
 in the resulting absolute adsroption at the new conditions."""
function solve_for_desorb(α_guess, P,  mol_co2, mol_n2, isotherm_CO2, isotherm_N2)
    #Find the partial pressures desorbing gas for the isotherms at Temperature (beta), 
     # and constant total presure and constant moles

    #Find the absolute adsorbed concentration of CO2 at the new (T, P) conditions
    #such that the gas given off to the outlet has the same concentration.
    function fun(y)
        #calculate loadings as a function of partial pressures:
        x = logistic(y) 
        partial_pressures = [x*P,(1-x)*P]
        
        iast_component_loadings = pyiast.iast(partial_pressures, 
                                                      [isotherm_CO2, isotherm_N2], 
                                                      verboseflag=false)
        n_co2 = iast_component_loadings[1] #absolute CO2 adsorption at proposed concentration
        n_n2 = iast_component_loadings[2] #absolute N2 adsorption at proposed concentration
        
        #Find the difference absolute loadings at previous desorption step  
        # and those at the new (T,P) conditions at the proposed concentration 
        # Then what concentration the gas would yield to the outlet.
        
        #Moles of gas released to outlet
        d_co2 = mol_co2 - n_co2
        d_n2 = mol_n2 - n_n2
        #Concentration of CO2 in outlet gas stream
        α_free_gas = d_co2/(d_co2+d_n2)
        
        #Difference between the concentration in the outlet
        # and the absolute adsorption concentration. 
        # In equlibrium, this difference will be zero.
        d_alpha = α_free_gas-x
        
        return d_alpha
    end
    
    #Find the new concentration of the released gas! 
    y = logit(α_guess)
    solution = find_zero(fun, y)
    alpha = logistic(solution)
    
    #calculate loadings as a function of partial pressures:
    partial_pressures = [alpha*P,(1-alpha)*P]
        
    iast_component_loadings = pyiast.iast(partial_pressures, 
                                                      [isotherm_CO2, isotherm_N2], 
                                                      verboseflag=false)
    n_co2 = iast_component_loadings[1]
    n_n2 = iast_component_loadings[2]
        
    #Find the difference and what concentration the gas would yield
    d_co2 = mol_co2 - n_co2
    d_n2 = mol_n2 - n_n2    

    return alpha, n_co2, n_n2, d_co2, d_n2
end

"""Function to calculate the equlibrium absolute adsorption, 
while desorbing from the initial conditions along the temperature and total pressure path."""
function Generate_sorption_path(Ts, Ps, α, Kh_CO2, Kh_N2, material)
    #Where T and P are the temperature [K] and total Pressure [Pa] steps 
    #alpha is the inlet CO2 concentration
    #Kh_N2, Kh_CO2 are the reslults of the Henry constant calculations {dict}
    #material is the dictionary of the material properties

    # Pressures = 10 .^ range(log10(0.01), log10(111458), 100) #[Pa] for iast calculations
    Pressures = 10 .^ range(log10(0.01), log10(111458*5), 150) #[Pa] for iast calculations
    
    βs = T_to_β.(Ts)

    #Adsorbtion:   Equilibrium with inlet T, P, and alpha
    #Generate an isotherm
    #Extrapolate the CO2 isotherm to the current β
    #convert to [μmol/(kg Pa)] for numerical stability
    Henry_CO2, Henry_CO2_err = Kh_extrapolate(βs[1], Kh_CO2, material) .* 1e6 #[nmol/(kg Pa)]
    PHenry_CO2 =  reshape(Henry_CO2 .* Pressures,:) #[nmol/kg]
    CO2_df = pd.DataFrame( Dict( "P"=> Pressures,
                                 "loading" => PHenry_CO2))
    isotherm_CO2 = pyiast.ModelIsotherm(CO2_df,
                                        loading_key = "loading",
                                        pressure_key = "P",
                                        model = "Henry")
        
    #Extrapolate the N2 isotherm to the current β
    Henry_N2, Henry_N2_err = Kh_extrapolate(βs[1], Kh_N2, material) .* 1e6 #[nmol/(kg Pa)]
    PHenry_N2 =  reshape(Henry_N2 .* Pressures,:) #[nmol/Pa]
    N2_df = pd.DataFrame( Dict( "P"=> Pressures,
                                 "loading" => PHenry_N2))
    isotherm_N2 = pyiast.ModelIsotherm(N2_df,
                                       loading_key = "loading",
                                       pressure_key = "P",
                                       model = "Henry")
                                
    partial_pressures = [α*Ps[1], (1-α)*Ps[1]]
    iast_component_loadings = pyiast.iast(partial_pressures, 
                                          [isotherm_CO2, isotherm_N2], 
                                          verboseflag=false)
    n_CO2 = [iast_component_loadings[1]] #[mmol/kg]
    d_CO2 = [NaN] #[mmol/kg]
    n_N2 = [iast_component_loadings[2]] #[mmol/kg]
    d_N2 = [NaN] #[mmol/kg]
    αs = [α]
    
    #Desorption: In equilibrium with the outlet
    for (i, (β, P)) in enumerate(zip(βs[2:end], Ps[2:end]))
        #Extrapolate the CO2 isotherm to the current β
        Henry_CO2, Henry_CO2_err = Kh_extrapolate(β, Kh_CO2, material) .* 1e6 #[nmol/(kg Pa)]
        PHenry_CO2 =  reshape(Henry_CO2 .* Pressures,:) #[nmol/kg]
        CO2_df = pd.DataFrame( Dict( "P"=> Pressures,
                                     "loading" => PHenry_CO2))
        isotherm_CO2 = pyiast.ModelIsotherm(CO2_df,
                                            loading_key = "loading",
                                            pressure_key = "P",
                                            model = "Henry")
        
        #Extrapolate the N2 isotherm to the current β
        Henry_N2, Henry_N2_err = Kh_extrapolate(β, Kh_N2, material) .* 1e6 #[nmol/(kg Pa)]
        PHenry_N2 =  reshape(Henry_N2 .* Pressures,:) #[nmol/kg]
        N2_df = pd.DataFrame( Dict( "P"=> Pressures,
                                     "loading" => PHenry_N2))
        isotherm_N2 = pyiast.ModelIsotherm(N2_df,
                                            loading_key = "loading",
                                            pressure_key = "P",
                                            model = "Henry")

        #Find conentration of released gas, and new adsorbed amounts of CO2 and N2
        mol_co2 = n_CO2[end] #[mmol/kg]
        mol_n2 = n_N2[end] #[mmol/kg]
        alpha_guess = αs[end]
        alpha_i, n_CO2_i, n_N2_i, d_CO2_i, d_N2_i = solve_for_desorb(alpha_guess, P, 
                                                                     mol_co2, mol_n2, 
                                                                     isotherm_CO2, isotherm_N2)
                                
        append!(n_CO2, [n_CO2_i]) #[mmol/kg]
        append!(d_CO2, [d_CO2_i]) #[mmol/kg]
        append!(n_N2, [n_N2_i]) #[mmol/kg]
        append!(d_N2, [d_N2_i]) #[mmol/kg]
        append!(αs, [alpha_i])
    end
    n_CO2 = n_CO2 .* 1e-6 #[mmol/kg]
    n_N2 = n_N2 .* 1e-6 #[mmol/kg]
    d_CO2 = d_CO2 .* 1e-6 #[mmol/kg]
    d_N2 = d_N2 .* 1e-6 #[mmol/kg]

    return n_CO2, n_N2, d_CO2, d_N2, αs
end

"""Function to calculate the equlibrium absolute adsorption, 
while desorbing from the initial conditions along the temperature and total pressure path."""
function Analytical_Henry_Generate_sorption_path(βs, Ps, α, Henry_CO2, Henry_N2)
    #Where β and P are the inverse temperature [mol/kJ] and total Pressure [Pa] steps 
    #alpha is the inlet CO2 concentration
    #Henry_CO2, Henry_CO2 are Henry constants along the path
    #material is the dictionary of the material properties

    #Adsorbtion:   Equilibrium with inlet T, P, and alpha
    n_CO2 = [Henry_CO2[1] * Ps[1] * α] #[mmol/kg]
    n_N2 = [Henry_N2[1] * Ps[1] * (1-α)] #[mmol/kg]
    d_CO2 = [NaN] #[mmol/kg]
    d_N2 = [NaN] #[mmol/kg]
    αs = [α]

    #Desorption: In equilibrium with the outlet
    for (β, P, henry_co2, henry_n2) in zip(βs[2:end], Ps[2:end], Henry_CO2[2:end], Henry_N2[2:end])

        A = henry_n2 * P - henry_co2 * P

        B = n_CO2[end] + n_N2[end] + henry_co2 * P - henry_n2 * P

        C = -1 * n_CO2[end]

        x1 = (-1 .* B + sqrt(B.^2 .- 4 .* A .* C))./(2 .* A)
        x2 = (-1 .* B - sqrt(B.^2 .- 4 .* A .* C))./(2 .* A)
        xs = [x1, x2]
 
        x = xs[argmin(abs.(xs .- αs[end]))]
    
        n_CO2_i = henry_co2 * P * x
        n_N2_i = henry_n2 * P * (1-x)
        d_CO2_i = n_CO2[end] - n_CO2_i
        d_N2_i = n_N2[end] - n_N2_i

        append!(n_CO2, [n_CO2_i]) #[mmol/kg]
        append!(d_CO2, [d_CO2_i]) #[mmol/kg]
        append!(n_N2, [n_N2_i]) #[mmol/kg]
        append!(d_N2, [d_N2_i]) #[mmol/kg]
        append!(αs, [x])

    end


    return n_CO2, n_N2, d_CO2, d_N2, αs

    
    
end

#test comment

