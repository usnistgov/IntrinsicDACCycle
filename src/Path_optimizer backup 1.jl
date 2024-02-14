

"""Function to optimize the start and end temperatures and pressures 
for a givin material and inlet CO2 concentration.
Assumes a linear path,"""
function Optimize_Intrinsic_Refresh(Base_directory::String, name::String,  
                                    α::Real)

    #Define the limits of the parameters	
    #specify Start and total Step
	T_start = 273.0 #[K] 
	ΔT = 200.0 #[K]

	T_start_lower = 200.0 #[K]
	T_start_upper = 400.0 #[K]

	ΔT_lower = 0.0 #[K]
	ΔT_upper = 200.0 #[K]

	P_start = 101325.0 #[Pa]
	#limit the ΔP to only reach "rough vaccuum" 100 Pa.
	ΔP = 100.0 - P_start #[Pa] the offset of 1 Pa ensures that the path of P never gets to 0 Pa. 

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

    """
    Define a function to try the parameters of the path
        and return the objectives and constraint costs (gx, hx)
    """ 
    function ScorePath(parameters)
        T_start = parameters[1]
        ΔT = parameters[2]
        P_start = parameters[3]
        ΔP = parameters[4]
        
        T_steps = length(T_start:0.5:T_start+ΔT)
        P_steps = length(P_start:-250:P_start+ΔP)

        steps = maximum([T_steps, P_steps])

        Ts = collect(LinRange(T_start, T_start+ΔT, steps))
        Ps = collect(LinRange(P_start, P_start+ΔP, steps))

        ξ, α_end =  IntrinsicDACCycle.Intrinisic_refresh_objectives(Base_directory, name,
                                                                Ts, Ps, α)
        
        objectives = [1/ξ, 1-α_end]
        gx = [0.0] # inequality constraints
        hx = [0.0] # equality constraints
        return objectives, gx, hx
    end

    #Perform the multi objective optimization
    N= 400
	method = Metaheuristics.NSGA3(N= N)
	Metaheuristics.optimize!(ScorePath, bounds, method)

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
	
		path_T_start = param[1]
		path_ΔT = param[2]
		path_P_start = param[3]
		path_ΔP = param[4]
	
		path_T_steps = length(path_T_start:0.5:path_T_start+path_ΔT)
		path_P_steps = length(path_P_start:-250:path_P_start+path_ΔP)
	
		path_steps = maximum([path_T_steps, path_P_steps])
	
		path_Ts[i] = collect(LinRange(path_T_start, path_T_start+path_ΔT, path_steps))
		path_Ps[i] = collect(LinRange(path_P_start, path_P_start+path_ΔP, path_steps))
	end

    #Initialize a dictionary to store all the results and intermediates
    Results_Dict = sort(Dict{String, Any}("Name" => name))
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

	Path_Dict["Henry_CO2"] = Vector{Matrix}(undef,num_solutions)
    Path_Dict["Henry_CO2_err"] = Vector{Matrix}(undef,num_solutions)
    Path_Dict["Henry_N2"] = Vector{Matrix}(undef,num_solutions)
    Path_Dict["Henry_N2_err"] = Vector{Matrix}(undef,num_solutions)
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
    end

    #Package all the dictionaries together
    E_Balance_Dict["Step_1"] = Step_1_Dict
    E_Balance_Dict["Step_2"] = Step_2_Dict
    E_Balance_Dict["Step_3"] = Step_3_Dict

    Results_Dict["Refresh_Path"] = Path_Dict
    Results_Dict["E_Balance"] = E_Balance_Dict

    return Results_Dict

end












# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""Function to find the set of pareto optimal paths
for the intrinsic refresh cycle of each material.
"""
# function find_pareto_path(directory::String, name::String, α)
#     #Specify the prior distributions of the path
#     steps = 202
#     T_start = Uniform(250, 350) #[K] Start anywhere between 250 K and 350 K. 
#     dT_s = Uniform(0,2, steps) #[K] Each of the steps can be between 0 K and 2 K.
#     P_start = Uniform(101325 *0.9, 101325 *1.1) #[Pa] Start anywhere between +/- 10% of 1 atm
#     dP_s = Uniform(0, -2000, steps) #[Pa] Each of the steps can be between 0 Pa and -2000 Pa.

#     #Make the path and return the objective funcitons
#     function path_evaluate(T_start, dT_s, P_start, dP_s)
#         #Make the path
#         Ts = T_start .+ cumsum(dT_s)
#         Ps = P_start .+ cumsum(dP_s)

#         #Perform the Intrinsic_refresh along that path
#         Results = Intrinisic_refresh_path(directory, name,
#                                           Ts, Ps, α)
#         #Objective to maximize:
#         objective1 = Results["Purity_captured_CO2"]
#         #Objective to maximize:
#         objective2 = Results["Intrinsic_capture_efficiency"]

#         return objective1, objective2, Results
#     end

#     #Find the pareto optimal path
#     pareto_paths = ParetoOptimize((T_start, dT_s, P_start, dP_s) -> path_evaluate(T_start, dT_s, P_start, dP_s), (1,2))

#     #Save the results on the set pareto paths
      



# end