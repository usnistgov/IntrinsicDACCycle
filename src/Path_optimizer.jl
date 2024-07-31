"""
Define a function to try the parameters of the path
    and return the objectives and constraint costs (gx, hx)
""" 
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
    # Metaheuristics doesn't handle `missing`s but does handle `Nan`s
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




"""Function to optimize the start and end temperatures and pressures 
for a given material and inlet CO2 concentration.
Assumes a linear path,"""
function Optimize_Intrinsic_Refresh(Base_directory::String, name::String,  
                                    α::Real)

	#Initialize a dictionary to store all the results
	Results_Dict = sort(Dict{String, Any}("Name" => name))
    
	#do the close_enough to linear adsorption test
	#and the saturation test before attempting the calculation
	material, Kh_N₂, Kh_CO₂, One_atm_N₂ = read_jsons(Base_directory,name)
	close_enough = Close_enough(Base_material, Kh_N₂, One_atm_N₂)
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
        ith_results = Intrinisic_refresh_path(Base_directory, name,
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
	material, Kh_N₂, Kh_CO₂, One_atm_N₂ = read_jsons(Base_directory,name)
	close_enough = Close_enough(material, Kh_N₂, One_atm_N₂)
	saturation_any_co2 = saturation_adsorb_co2(Base_directory, name)
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
            ith_results = Intrinisic_refresh_path(Base_directory, name,
                                                    trial_Ts, trial_Ps, α)

            objectives_dist = Intrinisic_refresh_objectives_posterior_dist(Base_directory, name,
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


