### A Pluto.jl notebook ###
# v0.19.47

using Markdown
using InteractiveUtils

# ╔═╡ fda35b50-71d7-11ef-23da-850f3a6532ea
using Pkg

# ╔═╡ 63ff5602-89e2-45ee-925c-543db2bc8e3e
Pkg.activate("/users/asm6/Julia_scripts/IntrinsicDACCycle")


# ╔═╡ b8f5f181-367d-41ca-b463-cc3dc8f3b3b6
begin
	using Plots
	using Distributed
	using JSON
	using DataFrames

end


# ╔═╡ bca614e9-ef9d-4eef-9e34-7b7f26c7b6cf
using LinearAlgebra

# ╔═╡ d866129f-b8d0-4996-ad83-c55569749fc8
using NaNStatistics

# ╔═╡ d139dcb6-66a4-402b-b46c-8e519909fe0d
using Metaheuristics 

# ╔═╡ 7ce56fb6-a859-4678-b1ae-57a179e648cd
using IntrinsicDACCycle

# ╔═╡ bb9814ac-5a96-4022-9011-22afee662fdb
Base_directory = "/users/asm6/DAC_data"

# ╔═╡ c8d4c5ec-5eb1-4856-9590-01efaeb97f98
begin
	#get all the material files
	list_of_material_files = filter(x -> occursin.(".json",x), readdir(Base_directory*"/CSD_FEASST_Materials/Materials/"))
	#strip off the .json tag
    list_of_materials = replace.(list_of_material_files, ".json" => "")
	#filter for _clean matierals
	# list_of_clean_materials = filter(x -> occursin.("_clean", x), list_of_materials)
end


# ╔═╡ 10515cb5-4c36-455c-9f5a-e5db0ed82e63
length(list_of_materials)

# ╔═╡ 4fa763cd-df93-4e58-845c-0a98d45ac858
begin
	#get all the Optimized Intrinsic Cycle files
	list_of_OpIDACC_files = filter(x -> occursin.(".json",x), readdir(Base_directory*"/Optimized_Intrinsic_Cycle/"))
	#strip off the .json tag
    list_of_OpIDACC = replace.(list_of_OpIDACC_files, ".json" => "")
	list_of_OpIDACC = replace.(list_of_OpIDACC, "OptIDC_" => "")
end

# ╔═╡ db3c9b0c-da6c-4fdb-9813-ac949c1b096c
length(list_of_OpIDACC)

# ╔═╡ 01cf8b49-387e-46aa-bf6c-9a6d828f33a5
list_of_OpIDACC

# ╔═╡ 70230652-6811-4911-ba9b-0d0113144f21
begin
	remaining_files = setdiff(list_of_materials, list_of_OpIDACC)
end

# ╔═╡ 3517897b-dc70-4db4-8e39-4e1d471fda74
length(remaining_files)

# ╔═╡ dfe4e04f-7b76-4691-abc6-bf9abc0f096c
function Run_an_Opitimize_IntrinsicDACCycle(name)
    directory = "/users/asm6/DAC_data/"
    α = 400/1000000
    results = IntrinsicDACCycle.Optimize_Intrinsic_Refresh_w_err(directory, name, α)

    filename = directory*"/Optimized_Intrinsic_Cycle/OptIDC_"*name*".json"
    open(filename,"w") do f
        JSON.print(f, results, 4)
    end
end


# ╔═╡ 525abae2-dd07-48be-8aaf-ebb37389a628
for file in remaining_files
	@show "running", file
	try
		Run_an_Opitimize_IntrinsicDACCycle(file)
		@show "done!"
	catch
		@show "Skipping"
	end
	
end

# ╔═╡ 1084bedb-a0f2-4c6e-a1ec-330bdb339613
remaining_files

# ╔═╡ 92718c52-6d3d-49a9-b0a2-94c18655200b


# ╔═╡ 279e087b-f1ca-4bdd-b1aa-4ee7121ac751
remaining_files[[1,3]]

# ╔═╡ 07d997f5-0758-48c7-b19e-40c57591bd21
begin
	directory = "/users/asm6/DAC_data/"
    α = 400/1000000


	
	#For each remaining file
	for name in remaining_files[[1,3]]
		complete_flag = 0
		@show name
		#try 10 times 
		for i in 1:10 
			@show complete_flag
			
			if complete_flag <= 10
				#calculate the optimized cycle
				global results = IntrinsicDACCycle.Optimize_Intrinsic_Refresh_w_err(directory, name, α)
				#iterate the flag
				complete_flag += 1
				#flag if a good calc
				if isassigned(results["Captured_CO2"])
					complete_flag = 1000
				end
			end
		end
		
		#Clean up the undefined results	
		if ~isassigned(results["Captured_CO2"])
			results["Captured_CO2"] = []
		end
	
		if ~isassigned(results["Captured_CO2_std"])
			results["Captured_CO2_std"] = []
		end
	
		if ~isassigned(results["Captured_N2"])
			results["Captured_N2"] = []
		end
	
		if ~isassigned(results["Captured_N2_std"])
			results["Captured_N2_std"] = []
		end
	
		if ~isassigned(results["Intrinsic_capture_efficiency"])
			results["Intrinsic_capture_efficiency"] = []
		end
	
		if ~isassigned(results["Intrinsic_capture_efficiency_std"])
			results["Intrinsic_capture_efficiency_std"] = []
		end
		
		if ~isassigned(results["Purity_captured_CO2"])
			results["Purity_captured_CO2"] = []
		end
	
		if ~isassigned(results["Purity_captured_CO2"])
			results["Purity_captured_CO2"] = []
		end
		
		if ~isassigned(results["Purity_captured_CO2_std"])
			results["Purity_captured_CO2_std"] = []
		end
	
		if ~isassigned(results["Captured_CO2"])
			results["Refresh_Path"] = []
		end
	
		if ~isassigned(results["Captured_CO2"])
			results["E_Balance"] = []
		end
	
		#Write the file to disc
		filename = directory*"/Optimized_Intrinsic_Cycle/OptIDC_"*name*".json"
	    open(filename,"w") do f
	        JSON.print(f, results, 4)
	    end
	end
end
		

# ╔═╡ ea14a359-fd6d-4401-9b68-4634fec3582b
length(results["Captured_CO2"])

# ╔═╡ 9bed4fa9-df5b-4c18-bfd4-d027c4c46a64
# begin
#     directory = "/users/asm6/DAC_data/"
#     α = 400/1000000
#     results = IntrinsicDACCycle.Optimize_Intrinsic_Refresh_w_err(directory, "DIBJIJ_clean", α)
# end

# ╔═╡ 9ea586d8-9e96-4f00-a4ed-005b898e6063
# begin
# 	if ~isassigned(results["Captured_CO2"])
# 		results["Captured_CO2"] = []
# 	end

# 	if ~isassigned(results["Captured_CO2_std"])
# 		results["Captured_CO2_std"] = []
# 	end

# 	if ~isassigned(results["Captured_N2"])
# 		results["Captured_N2"] = []
# 	end

# 	if ~isassigned(results["Captured_N2_std"])
# 		results["Captured_N2_std"] = []
# 	end

# 	if ~isassigned(results["Intrinsic_capture_efficiency"])
# 		results["Intrinsic_capture_efficiency"] = []
# 	end

# 	if ~isassigned(results["Intrinsic_capture_efficiency_std"])
# 		results["Intrinsic_capture_efficiency_std"] = []
# 	end
	
# 	if ~isassigned(results["Purity_captured_CO2"])
# 		results["Purity_captured_CO2"] = []
# 	end

# 	if ~isassigned(results["Purity_captured_CO2"])
# 		results["Purity_captured_CO2"] = []
# 	end

# 	if ~isassigned(results["Captured_CO2"])
# 		results["Refresh_Path"] = []
# 	end

# 	if ~isassigned(results["Captured_CO2"])
# 		results["E_balance"] = []
# 	end
# end

# ╔═╡ 418fa3b5-58ce-42d8-aa18-e4e009b1af21
# results

# ╔═╡ 93b65e57-c80a-49e0-ab4e-56a29010d1d4
# isassigned(results["Name"])

# ╔═╡ 6f931e86-767c-4113-a30a-47dd81af0159
begin
	thing = sort(Dict{String, Any}("some_propery" => []))
end

# ╔═╡ 546d3951-ce1a-434a-908a-a6930cb39482
~isassigned(thing["some_propery"])

# ╔═╡ 96fe7099-2607-45b4-8cd0-b84cfd33a38d
length(thing["some_propery"])

# ╔═╡ 46bdd140-0839-4133-b87d-28fecf144ec4
# begin
# 	name = "DIBJIJ_clean"
# 	material, Kh_N₂, Kh_CO₂, One_atm_N₂ = IntrinsicDACCycle.read_jsons(directory, name)

#     # #Test if the Henry constant at 1 atm is close enough to the direct GCMC at 1 atm
#     # close_enough_test = Close_enough(material, Kh_N₂, One_atm_N₂)
#     # #If not close enough
#     # if close_enough_test == false
#     #     return nothing
#     # end
    
#     T_start = 300
#     ΔT = 50
# 	dT_max = 0.25

#     P_start = 1.1 .* 101325.0
#     ΔP = 100.0 - P_start
#     dP_max = -125.0
        
# 	#Calculate the number of steps from start to start+Δ with the max T step or P step size
# 	T_steps = length(T_start:dT_max:T_start+ΔT)
# 	P_steps = length(P_start:dP_max:P_start+ΔP)

# 	#Choose the larger number of steps
# 	steps = maximum([T_steps, P_steps])
# 	#Create the T and P path with those steps
# 	Ts = collect(LinRange(T_start, T_start+ΔT, steps))
# 	Ps = collect(LinRange(P_start, P_start+ΔP, steps))
# 	# α = 400/1000000 #400 ppm is the concentration of CO2 in ambient air
#     βs = IntrinsicDACCycle.T_to_β.(Ts) #[mol/kJ]

#     #Test if the Henry constant at 1 atm is close enough to the direct GCMC at 1 atm
#     close_enough_test = IntrinsicDACCycle.Close_enough(material, Kh_N₂, One_atm_N₂)
#     #If not close enough
#     if close_enough_test == false
#         #set Ts and Ps and βs to NaN
#         Ts = Ts .* NaN
#         Ps = Ps .* NaN
#         βs = βs .* NaN
#     end

#     #Extrapolate Henry constants along the path
#     #Extrapolate the CO2 isotherm to the βs
#     Henry_CO2, Henry_CO2_err = IntrinsicDACCycle.Kh_extrapolate(βs, Kh_CO₂, material) #[mmol/(kg Pa)]

#     #Extrapolate the N2 isotherm to the βs
#     Henry_N2, Henry_N2_err = IntrinsicDACCycle.Kh_extrapolate(βs, Kh_N₂, material)  #[mmol/(kg Pa)]

#     #If extrapolating outside a reasonable range (for that material),
#     # the Henry Constants will un-physically increase with increasing temperature
#     #Keep only the monotonically decreasing parts of the Henry constants
#     mono_Henry_CO2, CO2_indices = IntrinsicDACCycle.keep_monotonic_decreasing(Henry_CO2)
#     mono_Henry_N2, N2_indices = IntrinsicDACCycle.keep_monotonic_decreasing(Henry_N2)
#     #Choose whichever set of indices is smaller
#     index_lenghts = [length(CO2_indices), length(N2_indices)]
# 	choice_of_indices = [CO2_indices, N2_indices]
# 	indices = choice_of_indices[argmin(index_lenghts)]

#     #Apply the index truncation to all relevant variables
#     Henry_CO2 = Henry_CO2[indices]
#     Henry_CO2_err = Henry_CO2_err[indices]

#     Henry_N2 = Henry_N2[indices]
#     Henry_N2_err = Henry_N2_err[indices]

#     Ts = Ts[indices]
#     Ps = Ps[indices]
#     βs = βs[indices]
# 	@show Ts
#     #Compare the CO2 Henry constant to Saturation uptake of CO2
#     #Truncate to a sensible range
#     Ts, Ps, βs, Henry_CO2, Henry_CO2_err, Henry_N2, Henry_N2_err = IntrinsicDACCycle.truncate_to_saturation(directory, name, α,
#                                                                                           Ts, Ps, βs, 
#                                                                                           Henry_CO2, Henry_CO2_err, 
#                                                                                           Henry_N2, Henry_N2_err) 
# end

# ╔═╡ 08f327a4-1f42-4a1a-b71e-0597142f1382


# ╔═╡ 3a767149-49e6-4c1c-af7d-34c230563ce8
# results["Captured_CO2"]*2

# ╔═╡ 7f5ec415-73d4-49db-a22e-cf3121134767


# """Function to optimize the start and end temperatures and pressures 
# for a given material and inlet CO2 concentration.
# Assumes a linear path,
# With the uncertainties in performance metrics"""
# function Optimize_Intrinsic_Refresh_w_err(Base_directory::String, name::String,  
#                                     α::Real)

    
#     #Initialize a dictionary to store all the results
# 	Results_Dict = sort(Dict{String, Any}("Name" => name))
    
# 	#do the close_enough to linear adsorption test
# 	#and the saturation test before attempting the calculation
# 	material, Kh_N₂, Kh_CO₂, One_atm_N₂ = IntrinsicDACCycle.read_jsons(Base_directory,name)
# 	close_enough = IntrinsicDACCycle.Close_enough(material, Kh_N₂, One_atm_N₂)
# 	saturation_any_co2 = IntrinsicDACCycle.saturation_adsorb_co2(Base_directory, name)
# 	continue_test = close_enough & saturation_any_co2
# 	Results_Dict["Close Enough to Linear Adsorption"] = close_enough
# 	Results_Dict["Saturation Adsorb Any CO2"] = saturation_any_co2
# 	if continue_test == false
# 		return Results_Dict
# 	end
# 	@show close_enough
# 	@show saturation_any_co2
    
#     #Define the limits of the parameters	
#     #specify Start and total Step
# 	T_start = 273.0 #[K] 
# 	ΔT = 200.0 #[K]
#     dT_max = 0.25

# 	T_start_lower = 200.0 #[K]
# 	T_start_upper = 400.0 #[K]

# 	ΔT_lower = 0.0 #[K]
# 	ΔT_upper = 200.0 #[K]

# 	P_start = 101325.0 #[Pa]
# 	#limit the ΔP to only reach "rough vaccuum" 100 Pa.
# 	ΔP = 100.0 - P_start #[Pa] the offset of 1 Pa ensures that the path of P never gets to 0 Pa. 
#     dP_max = -125.0

# 	P_start_lower = 101325.0 #[Pa] 
# 	P_start_upper = 1.1 .* [101325.0] #[Pa]

# 	#Lower limit of ΔP is to reach "rough vaccum" 100 Pa
# 	ΔP_lower = 100.0 - P_start #[Pa]
# 	# ΔP_lower = (500 .- P_start[1])./steps .+ zero(ΔP) #[Pa]
# 	ΔP_upper = 0.0 #[Pa]

#     #Combine the bounds into one object for pasing to optimizer
#     lower_bound = cat(T_start_lower, ΔT_lower, P_start_lower, ΔP_lower, dims = 1)
# 	upper_bound = cat(T_start_upper, ΔT_upper, P_start_upper, ΔP_upper, dims = 1)
# 	bounds = cat(lower_bound, upper_bound, dims = 2)'

#     #Perform the multi objective optimization
# 	@show "trying optimization"
#     N= 800
# 	method = Metaheuristics.NSGA3(N= N)
#     function f(x)
#         return ScorePath(x, dT_max, dP_max, Base_directory, name, α, "objectives", 0)
#     end

# 	Metaheuristics.optimize!(f, bounds, method)
# 	@show "unpacking optimization"
#     #Unpack the results (Pareto front of perfrormance metrics and their parameters)
#     results_state = Metaheuristics.get_result(method)
# 	results = pareto_front(results_state)
    
# 	non_dom_solutions = Metaheuristics.get_non_dominated_solutions(method.status.population)
# 	num_solutions = length(non_dom_solutions)

# 	results_parameters = Vector{Vector}(undef, num_solutions)

# 	for (i,solution) in enumerate(non_dom_solutions)
		
# 		results_parameters[i] = solution.x
# 	end
# 	results_parameters = cat(results_parameters..., dims = 2)'
    
#     path_Ts = Vector{Vector}(undef, length(results_parameters[:,1]))
# 	path_Ps = Vector{Vector}(undef, length(results_parameters[:,1]))
# 	for (i,param) in enumerate(eachrow(results_parameters))

#         #To recover the stepsize run ScorePath again
#         objectives = ScorePath(param, dT_max, dP_max, Base_directory, name, α, "steps")
#         dT_max = objectives[4]
#         dP_max = objectives[5]
	
# 		path_T_start = param[1]
# 		path_ΔT = param[2]
# 		path_P_start = param[3]
# 		path_ΔP = param[4]
	
# 		path_T_steps = length(path_T_start:dT_max:path_T_start+path_ΔT)
# 		path_P_steps = length(path_P_start:-dP_max:path_P_start+path_ΔP)
	
# 		path_steps = maximum([path_T_steps, path_P_steps])
	
# 		path_Ts[i] = collect(LinRange(path_T_start, path_T_start+path_ΔT, path_steps))
# 		path_Ps[i] = collect(LinRange(path_P_start, path_P_start+path_ΔP, path_steps))
# 	end

#     #Initialize a dictionary to store all the intermediates
#     Path_Dict = sort(Dict{String, Any}("Refresh_Path" => "Definition of refresh path and material properties along that path (in per kg of sorbent basis)."))
#     E_Balance_Dict = sort(Dict{String, Any}("E_Balance" => "Energy balance along path"))
#     Step_1_Dict = sort(Dict{String, Any}("Step_1" => "Adsorption"))
#     Step_2_Dict = sort(Dict{String, Any}("Step_2" => "Desorption"))
#     Step_3_Dict = sort(Dict{String, Any}("Step_3" => "Waste energy recovery"))

# 	Path_Dict["Temperatures"] = Vector{Vector}(undef,num_solutions)
#     Path_Dict["Temperature_units"] = "K"
#     Path_Dict["Pressures"] = Vector{Vector}(undef,num_solutions)
#     Path_Dict["Pressure_units"] = "Pa"
#     Path_Dict["Betas"] = Vector{Vector}(undef,num_solutions)
#     Path_Dict["Beta_units"] = "mol/kJ"

# 	Path_Dict["Henry_CO2"] = Vector{Vector}(undef,num_solutions)
#     Path_Dict["Henry_CO2_err"] = Vector{Vector}(undef,num_solutions)
#     Path_Dict["Henry_N2"] = Vector{Vector}(undef,num_solutions)
#     Path_Dict["Henry_N2_err"] = Vector{Vector}(undef,num_solutions)
#     Path_Dict["Henry_units"] = "mmol/(kg Pa)"

#     Path_Dict["Moles_CO2"] = Vector{Vector}(undef,num_solutions)
#     Path_Dict["Moles_N2"] = Vector{Vector}(undef,num_solutions)
#     Path_Dict["Moles_units"] = "mol/kg"

#     Path_Dict["Heat_of_adsorb_CO2"] = Vector{Matrix}(undef,num_solutions)
#     Path_Dict["Heat_of_adsorb_CO2_err"] = Vector{Matrix}(undef,num_solutions)
#     Path_Dict["Heat_of_adsorb_N2"] = Vector{Matrix}(undef,num_solutions)
#     Path_Dict["Heat_of_adsorb_N2_err"] = Vector{Matrix}(undef,num_solutions)
#     Path_Dict["Heat_of_adsorb_units"] = "J/mol"

#     Path_Dict["Specific_heat_sorbent"] = Vector{Matrix}(undef,num_solutions)
#     Path_Dict["Specific_heat_sorbent_err"] = Vector{Matrix}(undef,num_solutions)
#     Path_Dict["Specific_heat_sorbent_units"] = "J/(kg K)"

#     Path_Dict["Flag"] = Vector{String}(undef,num_solutions)

#     Step_1_Dict["Heat_to_adsorb_CO2"] = Vector{Real}(undef,num_solutions)
#     Step_1_Dict["Heat_to_adsorb_N2"] = Vector{Real}(undef,num_solutions)
#     Step_1_Dict["Work_to_adsorb_CO2"] = Vector{Real}(undef,num_solutions)
#     Step_1_Dict["Work_to_adsorb_N2"] = Vector{Real}(undef,num_solutions)
#     Step_1_Dict["E_units"] = "J/kg_sorb"

#     E_Balance_Dict["E1"] = Vector{Real}(undef,num_solutions)

#     Step_2_Dict["Heat_to_desorb_CO2"] = Vector{Matrix}(undef,num_solutions)
#     Step_2_Dict["Heat_to_desorb_N2"] = Vector{Matrix}(undef,num_solutions)
#     Step_2_Dict["Work_to_desorb_CO2"] = Vector{Vector}(undef,num_solutions)
#     Step_2_Dict["Work_to_desorb_N2"] = Vector{Vector}(undef,num_solutions)
#     Step_2_Dict["E_to_heat_adsorbed_CO2"] = Vector{Vector}(undef,num_solutions)
#     Step_2_Dict["E_to_heat_adsorbed_N2"] = Vector{Vector}(undef,num_solutions)
#     Step_2_Dict["E_to_heat_sorbent"] = Vector{Matrix}(undef,num_solutions)
#     Step_2_Dict["E_to_change_pressure"] = Vector{Vector}(undef,num_solutions)
#     Step_2_Dict["E_units"] = "J/kg_sorb"

#     E_Balance_Dict["E2"] = Vector{Real}(undef,num_solutions)

#     Step_3_Dict["E_recovered"] = Vector{Real}(undef,num_solutions)
#     E_Balance_Dict["E3"] = Vector{Real}(undef,num_solutions)

#     E_Balance_Dict["Total_E_of_cycle"] = Vector{Real}(undef,num_solutions)
#     E_Balance_Dict["E_units"] = "J/kg_sorb"

#     Results_Dict["Captured_CO2"] = Vector{Real}(undef,num_solutions)
#     Results_Dict["Captured_CO2_std"] = Vector{Real}(undef,num_solutions)
#     Results_Dict["Captured_N2"] = Vector{Real}(undef,num_solutions)
#     Results_Dict["Captured_N2_std"] = Vector{Real}(undef,num_solutions)
#     Results_Dict["Captured_gas_units"] = "mol/kg_sorb"

#     Results_Dict["Intrinsic_capture_efficiency"] = Vector{Real}(undef,num_solutions)
#     Results_Dict["Intrinsic_capture_efficiency_std"] = Vector{Real}(undef,num_solutions)
#     Results_Dict["Intrinsic_capture_efficiency_units"] = "mol/J"
#     Results_Dict["Purity_captured_CO2"] = Vector{Real}(undef,num_solutions)
#     Results_Dict["Purity_captured_CO2_std"] = Vector{Real}(undef,num_solutions)

#     #Re-evalutate the path at just the pareto front parameters
#     for i in 1:num_solutions
#         #Re-evalutate the path
#         # if there is a sampling that give a domain error, re-try with smaller step sizes a few times

#         trial_Ts = path_Ts[i]
#         trial_Ps = path_Ps[i]

#         #Initialize containers for the results
#         ith_results = nothing
#         objectives_dist = nothing

#         for attempts in 1:10
#             ith_results = IntrinsicDACCycle.Intrinisic_refresh_path(Base_directory, name,
#                                                     trial_Ts, trial_Ps, α)

#             objectives_dist = IntrinsicDACCycle.Intrinisic_refresh_objectives_posterior_dist(Base_directory, name,
#                                                                                     trial_Ts, trial_Ps, α,
#                                                                                     100)
#             if "Refresh_Path" in keys(ith_results)
#                 #Test if Ts_s are not NaNs (i.e. failed Close_enough test earlier)
#                 test_1 = ~any(isnan.(path_Ts[i]))
#                 #Test for domain error
#                 test_2 = any(isnan.(objectives_dist[1]))
#                 #if it passed the Close_enough test but still created a domain error 
#                 if test_1 & test_2
#                     new_Ts_start = path_Ts[i][1]
#                     new_Ts_end = path_Ts[i][end]
#                     #Create new T steps at finer resultion
#                     new_dT = (path_Ts[i][2] - path_Ts[i][1])*0.75^attempts

#                     new_Ps_start = path_Ps[i][1]
#                     new_Ps_end = path_Ps[i][end]
#                     #Create new P steps at finer resoluion
#                     new_dP = (path_Ps[i][2] - path_Ps[i][1]) *0.75^attempts

#                     new_path_T_steps = length(new_Ts_start:new_dT:new_Ts_end)
#                     new_path_P_steps = length(new_Ps_start:new_dP:new_Ps_end)

#                     new_path_steps = maximum([new_path_T_steps, new_path_P_steps])

                
#                     #Re-define the path Ts and Ps
#                     trial_Ts = collect(LinRange(new_Ts_start, new_Ts_end, new_path_steps))
#                     trial_Ps = collect(LinRange(new_Ps_start, new_Ps_end, new_path_steps))
#                     print("Uncertainty calculation lead to numerical instabilty, re-trying at finer step size.")
#                 else
#                     break
#                 end
#             end
#         end

#         if "Refresh_Path" in keys(ith_results)
#             mean_ξ = mean(objectives_dist[1])
#             std_ξ = std(objectives_dist[1])
#             mean_α = mean(objectives_dist[2])
#             std_α = std(objectives_dist[2])
            
#             mean_Δn_CO2 = mean(objectives_dist[3])
#             std_Δn_CO2 = std(objectives_dist[3])
#             mean_Δn_N2 = mean(objectives_dist[4])
#             std_Δn_N2 = std(objectives_dist[4])
            

#             #Re-package resutls
#             ith_path_dict = ith_results["Refresh_Path"]
#             ith_e_balance_dict = ith_results["E_Balance"]
#             ith_step_1_dict = ith_e_balance_dict["Step_1"]
#             ith_step_2_dict = ith_e_balance_dict["Step_2"]
#             ith_step_3_dict = ith_e_balance_dict["Step_3"]

#             Path_Dict["Temperatures"][i] = ith_path_dict["Temperatures"]
#             Path_Dict["Pressures"][i] = ith_path_dict["Pressures"]
#             Path_Dict["Betas"][i] = ith_path_dict["Betas"]
            
#             Path_Dict["Henry_CO2"][i] = ith_path_dict["Henry_CO2"]
#             Path_Dict["Henry_CO2_err"][i] = ith_path_dict["Henry_CO2_err"]
#             Path_Dict["Henry_N2"][i] = ith_path_dict["Henry_N2"]
#             Path_Dict["Henry_N2_err"][i] = ith_path_dict["Henry_N2_err"]
            
#             Path_Dict["Moles_CO2"][i] = ith_path_dict["Moles_CO2"]
#             Path_Dict["Moles_N2"][i] = ith_path_dict["Moles_N2"]
            
#             Path_Dict["Heat_of_adsorb_CO2"][i] = ith_path_dict["Heat_of_adsorb_CO2"]
#             Path_Dict["Heat_of_adsorb_CO2_err"][i] = ith_path_dict["Heat_of_adsorb_CO2_err"]
#             Path_Dict["Heat_of_adsorb_N2"][i] = ith_path_dict["Heat_of_adsorb_N2"]
#             Path_Dict["Heat_of_adsorb_N2_err"][i] = ith_path_dict["Heat_of_adsorb_N2_err"]
            
#             Path_Dict["Specific_heat_sorbent"][i] = ith_path_dict["Specific_heat_sorbent"]
#             Path_Dict["Specific_heat_sorbent_err"][i] = ith_path_dict["Specific_heat_sorbent_err"]
            
#             Step_1_Dict["Heat_to_adsorb_CO2"][i] = ith_step_1_dict["Heat_to_adsorb_CO2"]
#             Step_1_Dict["Heat_to_adsorb_N2"][i] = ith_step_1_dict["Heat_to_adsorb_N2"]
#             Step_1_Dict["Work_to_adsorb_CO2"][i] = ith_step_1_dict["Work_to_adsorb_CO2"]
#             Step_1_Dict["Work_to_adsorb_N2"][i] = ith_step_1_dict["Work_to_adsorb_N2"]

#             E_Balance_Dict["E1"][i] = ith_e_balance_dict["E1"]
            
#             Step_2_Dict["Heat_to_desorb_CO2"][i] = ith_step_2_dict["Heat_to_desorb_CO2"]
#             Step_2_Dict["Heat_to_desorb_N2"][i] = ith_step_2_dict["Heat_to_desorb_N2"]
#             Step_2_Dict["Work_to_desorb_CO2"][i] = ith_step_2_dict["Work_to_desorb_CO2"]
#             Step_2_Dict["Work_to_desorb_N2"][i] = ith_step_2_dict["Work_to_desorb_N2"]
#             Step_2_Dict["E_to_heat_adsorbed_CO2"][i] = ith_step_2_dict["E_to_heat_adsorbed_CO2"]
#             Step_2_Dict["E_to_heat_adsorbed_N2"][i] = ith_step_2_dict["E_to_heat_adsorbed_N2"]
#             Step_2_Dict["E_to_heat_sorbent"][i] = ith_step_2_dict["E_to_heat_sorbent"]
#             Step_2_Dict["E_to_change_pressure"][i] = ith_step_2_dict["E_to_change_pressure"]
            
#             E_Balance_Dict["E2"][i] = ith_e_balance_dict["E2"]
        
#             Step_3_Dict["E_recovered"][i] = ith_step_3_dict["E_recovered"]
#             E_Balance_Dict["E3"][i] = ith_e_balance_dict["E3"]
        
#             E_Balance_Dict["Total_E_of_cycle"][i] = ith_e_balance_dict["Total_E_of_cycle"]
            
#             Results_Dict["Captured_CO2"][i] = mean_Δn_CO2
#             Results_Dict["Captured_CO2_std"][i] = std_Δn_CO2
#             Results_Dict["Captured_N2"][i] = mean_Δn_N2
#             Results_Dict["Captured_N2_std"][i] = std_Δn_N2
            
#             Results_Dict["Intrinsic_capture_efficiency"][i] = mean_ξ
#             Results_Dict["Intrinsic_capture_efficiency_std"][i] = std_ξ
#             Results_Dict["Purity_captured_CO2"][i] = mean_α
#             Results_Dict["Purity_captured_CO2_std"][i] = std_α

#             Path_Dict["Flag"][i] = "no flags"
#         else
# 			Path_Dict["Flag"][i] = "No path, check Henry constant or saturation calculations"
# 		end
#     end

#     #Package all the dictionaries together
#     E_Balance_Dict["Step_1"] = Step_1_Dict
#     E_Balance_Dict["Step_2"] = Step_2_Dict
#     E_Balance_Dict["Step_3"] = Step_3_Dict

#     Results_Dict["Refresh_Path"] = Path_Dict
#     Results_Dict["E_Balance"] = E_Balance_Dict

#     return Results_Dict

# end




# ╔═╡ f155ff38-7c69-41b7-8e29-2f82f043939a
# """Function to read in the GCMC simulation results,
# and perform the full intrinsic refresh cycle analysis
# Along an arbitrary path in (Temperature,Pressure)-space.
# And return only the performance metrics."""
# function Intrinisic_refresh_objectives2(directory::String, name::String, 
#                                  Ts::AbstractArray, Ps::AbstractArray, 
#                                  α)
#     #Ts is a 1D array of Temperatures [K] along the path
#     #Ps is a 1D array of Total Pressure [Pa] along the path
#     #α is a scalar of the Concentration of CO2 in a mixture with N2 [mol/mol]
	
#     #Read in all the GCMC results
#     material, Kh_N₂, Kh_CO₂, One_atm_N₂ = IntrinsicDACCycle.read_jsons(directory, name)
#     # #Test if the Henry constant at 1 atm is close enough to the direct GCMC at 1 atm
#     # close_enough_test = Close_enough(material, Kh_N₂, One_atm_N₂)
#     # #If not close enough
#     # if close_enough_test == false
#     #     return nothing
#     # end

#     #Convert the Ts to inverse temperature β
#     βs = IntrinsicDACCycle.T_to_β.(Ts) #[mol/kJ]

#     #Test if the Henry constant at 1 atm is close enough to the direct GCMC at 1 atm
#     close_enough_test = IntrinsicDACCycle.Close_enough(material, Kh_N₂, One_atm_N₂)
#     #If not close enough
#     if close_enough_test == false
#         #set Ts and Ps and βs to NaN
#         Ts = Ts .* NaN
#         Ps = Ps .* NaN
#         βs = βs .* NaN
#     end
#     #Extrapolate Henry constants along the path
#     #Extrapolate the CO2 isotherm to the βs
#     Henry_CO2, Henry_CO2_err = IntrinsicDACCycle.Kh_extrapolate(βs, Kh_CO₂, material) #[mmol/(kg Pa)]

#     #Extrapolate the N2 isotherm to the βs
#     Henry_N2, Henry_N2_err = IntrinsicDACCycle.Kh_extrapolate(βs, Kh_N₂, material)  #[mmol/(kg Pa)]

#     #If extrapolating outside a reasonable range (for that material),
#     # the Henry Constants will un-physically increase with increasing temperature
#     #Keep only the monotonically decreasing parts of the Henry constants
#     mono_Henry_CO2, CO2_indices = IntrinsicDACCycle.keep_monotonic_decreasing(Henry_CO2)
#     mono_Henry_N2, N2_indices = IntrinsicDACCycle.keep_monotonic_decreasing(Henry_N2)
#     #Choose whichever set of indices is smaller
#     index_lenghts = [length(CO2_indices), length(N2_indices)]
# 	choice_of_indices = [CO2_indices, N2_indices]
# 	indices = choice_of_indices[argmin(index_lenghts)]

#     #Apply the index truncation to all relevant variables
#     Henry_CO2 = Henry_CO2[indices]
#     Henry_CO2_err = Henry_CO2_err[indices]
	
#     Henry_N2 = Henry_N2[indices]
#     Henry_N2_err = Henry_N2_err[indices]



#     Ts = Ts[indices]
#     Ps = Ps[indices]
#     βs = βs[indices]

#     #Compare the CO2 Henry constant to Saturation uptake of CO2
#     #Truncate to a sensible range
#     Ts, Ps, βs, Henry_CO2, Henry_CO2_err, Henry_N2, Henry_N2_err = IntrinsicDACCycle.truncate_to_saturation(directory, name, α,
#                                                                                           Ts, Ps, βs, 
#                                                                                           Henry_CO2, Henry_CO2_err, 
#                                                                                           Henry_N2, Henry_N2_err) 

	

# 	# @show size(Ts)[1]
#     #If the truncation has eliminated the path, skip the intrinisic refresh
#     if size(Ts)[1] == 0
#         objectives = [missing, missing]
#         return objectives
#     end

#     #Generate Equilibrium loadings along the path
# 	# @show maximum(Henry_CO2)
# 	# @show maximum(Henry_N2)
#     n_CO2, n_N2, d_CO2, d_N2, αs = IntrinsicDACCycle.Analytical_Henry_Generate_sorption_path(βs, Ps, α, Henry_CO2, Henry_N2) #[mmol/kg]
# 	# @show size(αs)
#     n_CO2 *= 10^-3 #convert to [mol/kg]
#     n_N2 *= 10^-3 #convert to [mol/kg]
#     d_CO2 *= 10^-3 #convert to [mol/kg]
#     d_N2 *= 10^-3 #convert to [mol/kg]
#     @show "Ah Haa0"
#     #Generate heat of adsorption along the path
# 	# @show βs
# 	# @show Kh_CO₂
# 	@show Ts
#     q_CO2, q_CO2_err = qₐ∞2(βs, Kh_CO₂) #kJ/mol of gas
#     q_CO2  *= 10^3 #[J/mol]
#     q_CO2_err  *= 10^3 #[J/mol]

	
# 	@show "Ah Haa0.5"
#     q_N2, q_N2_err = qₐ∞2(βs, Kh_N₂) #kJ/mol of gas
#     q_N2  *= 10^3 #[J/mol]
#     q_N2_err  *= 10^3 #[J/mol]

# 	@show "Ahh Haaa1!!"
#     #Generate specific heat of sorbent along the path
#     cv_s, cv_s_err =  IntrinsicDACCycle.Extrapolate_Cv(directory, name, Ts) #[J/(kg K)]
# 	@show size(cv_s)
#     #Energy balance for step 1
#     (Q_adsorb_CO2, Q_adsorb_N2, 
#     W_adsorb_CO2, W_adsorb_N2) = IntrinsicDACCycle.intrinsic_refresh_step_1(Ts, 
#                                                     n_CO2, n_N2,
#                                                     q_CO2, q_N2)
#     # [J/kg_sorb]
# 	@show "Ahh Haaa2!!"
#     E1 = Q_adsorb_CO2 +Q_adsorb_N2 + W_adsorb_CO2 + W_adsorb_N2 # [J/kg_sorb]

#     #Energy balance for step 2
#     (Q_CO2, Q_N2, 
#     W_desorb_CO2, W_desorb_N2, 
#     E_heat_ads_CO2, E_heat_ads_N2, 
#     E_heat_sorb, E_P) = IntrinsicDACCycle.intrinsic_refresh_step_2(Ts, Ps, 
#                                                  n_CO2, n_N2, d_CO2, d_N2,
#                                                  q_CO2, q_N2,
#                                                  cv_s)
#     # [J/kg_sorb]
#     @show "Ahh haaa3"
#     E2 = nansum(Q_CO2 .+ Q_N2 
#                 .+ W_desorb_CO2 .+ W_desorb_N2 
#                 .+ E_heat_ads_CO2 .+ E_heat_ads_N2 
#                 .+ E_heat_sorb .+ E_P)   # [J/kg_sorb]
#     @show "Ahh haaa4"
#     #Energy balance for step 3
#     E3 = 0 # [J/kg_sorb]

#     #Total Energy of refresh cycle
#     E = E1 + E2 + E3# [J/kg_sorb]

#     #Total captureed CO2 and N2
#     Δn_CO2 = n_CO2[1] - n_CO2[end] # [mol/kg_sorb]
#     Δn_N2 = n_N2[1] - n_N2[end] # [mol/kg_sorb]

#     #Calculate performance metrics 
#     Intrinsic_capture_efficiency = Δn_CO2/E #[mol/J]
#     Purity_captured_CO2 = Δn_CO2/(Δn_CO2 + Δn_N2) #[]
	

#     objectives = [Intrinsic_capture_efficiency, Purity_captured_CO2]
# 	@show objectives
#     return objectives
# end

# ╔═╡ b0bfa892-4f31-498f-b272-c07587544489
# """This function takes the predictions of the heat of adsorption at infinate dilution,
# and extrapolates to new thermodyanmic β."""
# function qₐ∞2(β, Kh_results)
#     #Calculate the heat of adsorption at infinate dilution
#     #Define the factor for the t-distribution for standard deivation
#     tvalue = 1.959
    
#     β₀ = Kh_results["beta"] #[mol/kJ]
# 	@show β₀
#     δᵦ = β .- β₀ #Distance in beta-space

#     Kcoeff = Kh_results["Kcoeff"] #K* coefficients  
#     power = length(Kcoeff)
#     powers = collect(range(0, power-1, power)) #the polynomial degree for each coeff
#     reps = length(δᵦ)
#     power_reps = ones((reps,1)) .* powers' #Matrix of each degree for each δᵦ
#     δᵦ_m =  δᵦ .* ones(power)' #matrix of each  δᵦ
#     # δᵦᴾ = δᵦ' .^ power_reps #each δᵦ to each power.
#     δᵦᴾ = δᵦ_m .^ power_reps #each δᵦ to each power.

#     kᵦ = δᵦᴾ .* Kcoeff'

#     denominator = sum(kᵦ, dims = 2) #Sum the terms of the polynomial for each δᵦ
#     @show size(denominator)
#     Kcoeff_1 = Kcoeff[2:end]
#     powers_1 = powers[2:end]
#     power_reps_1 = power_reps[:, 2:end]
#     δᵦ_m_1 =  δᵦ .* ones(power-1)'
#     δᵦᴾ_1 = δᵦ_m_1 .^(power_reps_1 .- 1) 
#     kᵦ_1 = δᵦᴾ_1 .* Kcoeff_1' .* powers_1'

#     numerator = sum(kᵦ_1, dims =2)
# 	@show size(numerator)
#     q_ads_∞ = 1 ./β .+ (numerator ./ denominator)

#     #Uncertainty Propagation:
#     #Add a leading zero to the k_beta for df/dKh0
    
#     kᵦ_1_0 = hcat(reshape(0 .* β', (length(β), 1)), kᵦ_1)
#     A = sum(kᵦ_1_0, dims = 2)
#     B = sum(kᵦ, dims = 2)
    
#     first_term = power_reps .* B .* δᵦ_m .^(power_reps .-1)
	
#     @show size(first_term)
#     second_term = A .* δᵦ_m .^(power_reps)
#     @show size(second_term)
#     Jacobian = (first_term .- second_term) ./(B.^2)
#     #when δᵦ = 0 use q = 1/β + Kh2/Kh1
#     #^ in which case the Jacobian becomes:
#     Jacobian_0 = zeros(power)
#     Jacobian_0[1] = -Kcoeff[2]/(Kcoeff[1]^2)
#     Jacobian_0[2] = 1/Kcoeff[1]
#     #find indexes where δ_beta == 0
#     mask_δᵦ = findall(x -> x == 0, reshape(δᵦ, (:)))
#     Jacobian[mask_δᵦ, CartesianIndex.(1:power)] .= reshape(Jacobian_0, (1,:)) 
     
# 	@show size(Jacobian)
#     Covariance = hcat(Kh_results["Kcoeff_covar"]...)
# 	@show size(Covariance)
#     q_ads_∞_var = Jacobian*Covariance*Jacobian' #The whole covariance matrix for the results
# 	@show size(q_ads_∞_var)
#     trials = Kh_results["trials"]
#     @show trials
# 	if minimum(diag(q_ads_∞_var)) < 0.0
# 		Ts_for_show = IntrinsicDACCycle.β_to_T.(β)
# 		@show β[1]
# 		@show Ts_for_show[1]
# 		@show Ts_for_show
# 		@show diag(q_ads_∞_var)
# 	end
# 	# @show maximum(diag(q_ads_∞_var))
#     q_ads_∞_std = reshape(sqrt.(Complex.(diag(q_ads_∞_var)./trials))*tvalue, (:,1)) #get the std from the covariance matrix
#     @show size(q_ads_∞_std)
#     return q_ads_∞, q_ads_∞_std #kJ/mol of gas
# end

# ╔═╡ fdc614e6-f1c4-41a2-a8b0-2d1342a05646
# """
# Define a function to try the parameters of the path
#     and return the objectives and constraint costs (gx, hx)
# """ 
# function ScorePath(parameters, dT_max, dP_max, 
#                    Base_directory, name, α,  
#                    return_flag = "objectives", i = 0)

#     T_start = parameters[1]
#     ΔT = parameters[2]
#     P_start = parameters[3]
#     ΔP = parameters[4]
    
#     outputs = try
        
#         #Calculate the number of steps from start to start+Δ with the max T step or P step size
#         T_steps = length(T_start:dT_max:T_start+ΔT)
#         P_steps = length(P_start:dP_max:P_start+ΔP)

#         #Choose the larger number of steps
#         steps = maximum([T_steps, P_steps])
#         #Create the T and P path with those steps
#         Ts = collect(LinRange(T_start, T_start+ΔT, steps))
#         Ps = collect(LinRange(P_start, P_start+ΔP, steps))

#         #Perform the Intrinsic Refresh Analysis with those steps
#         ξ, α_end =  Intrinisic_refresh_objectives2(Base_directory, name,
#                                                                 Ts, Ps, α)
# 		@show ξ
# 		@show α_end

#         (ξ, α_end, dT_max, dP_max, i)
#     catch e
#         if isa(e, DomainError) 
#             #If there is a DomainError it's likely from a numerical instabilty with IAST 
#             # So take a smaller step size and try again. Do this up to 10 times
#             if i <= 10
#                 dT_max *= 0.75
#                 dP_max *= 0.75
#                 i += 1
#                 println("Re-trying with smaller steps!")
#                 ScorePath(parameters, dT_max, dP_max, Base_directory, name, α, "outputs", i)
#             else
#                 println("Too many iterations, likely something other than stepsize is wrong")
#                 throw(e)
#             end
#         end
#     end

#     ξ = outputs[1]
#     α_end = outputs[2]
#     # Metaheuristics doesn't handle `missing`s but does handle `Nan`s
#     if ismissing(ξ)
#         ξ = NaN
#     end
#     if ismissing(α_end)
#         α_end = NaN
#     end

#     objectives = [1/ξ, 1-α_end]
#     gx = [0.0] # inequality constraints
#     hx = [0.0] # equality constraints

#     if return_flag == "outputs"
#         return outputs
#     elseif return_flag == "objectives"
#         return objectives, gx, hx
#     elseif return_flag == "steps"
#         dT_max = outputs[3]
#         dP_max = outputs[4]
#         return objectives, gx, hx, dT_max, dP_max
#     else 
#         return "Need to define return_flag"
#     end

# end



# ╔═╡ d0ddc413-6e2b-4691-ac7b-0f2a050160a9
# α = 400/1000000

# ╔═╡ e86cc213-7bdb-4dc0-a5ff-76bbbc3da97d
# Optimize_Intrinsic_Refresh_w_err(Base_directory, remaining_files[1], α)

# ╔═╡ bd59fe52-8dc8-4674-bdbf-dabdcabe43df
# begin
# 	Ts_for_extrapolation = collect(389:0.1:500)
# 	Ps_for_extrapolation = range(101325.0, 0.5*101325.0, length(Ts_for_extrapolation))

# 	material, Kh_N₂, Kh_CO₂, One_atm_N₂ = IntrinsicDACCycle.read_jsons(Base_directory, remaining_files[2])

	
# 	βs_for_extrapolation = IntrinsicDACCycle.T_to_β.(Ts_for_extrapolation)

# 	Henry_CO2, Henry_CO2_err = IntrinsicDACCycle.Kh_extrapolate(βs_for_extrapolation, Kh_CO₂, material)

# 	Henry_N2, Henry_N2_err = IntrinsicDACCycle.Kh_extrapolate(βs_for_extrapolation, Kh_N₂, material)

# 	mono_Henry_CO2, CO2_indices = IntrinsicDACCycle.keep_monotonic_decreasing(Henry_CO2)
#     mono_Henry_N2, N2_indices = IntrinsicDACCycle.keep_monotonic_decreasing(Henry_N2)
#     #Choose whichever set of indices is smaller
#     index_lenghts = [length(CO2_indices), length(N2_indices)]
# 	choice_of_indices = [CO2_indices, N2_indices]
# 	indices = choice_of_indices[argmin(index_lenghts)]

#     #Apply the index truncation to all relevant variables
#     Henry_CO2 = Henry_CO2[indices]
#     Henry_CO2_err = Henry_CO2_err[indices]

#     Henry_N2 = Henry_N2[indices]
#     Henry_N2_err = Henry_N2_err[indices]
	
#     Ts_for_extrapolation = Ts_for_extrapolation[indices]
#     Ps_for_extrapolation = Ps_for_extrapolation[indices]
#     βs_for_extrapolation = βs_for_extrapolation[indices]

#     #Compare the CO2 Henry constant to Saturation uptake of CO2
#     #Truncate to a sensible range
#     Ts_for_extrapolation, Ps_for_extrapolation, βs_for_extrapolation, Henry_CO2, Henry_CO2_err, Henry_N2, Henry_N2_err = IntrinsicDACCycle.truncate_to_saturation(Base_directory, remaining_files[2], α,
#                                                                                           Ts_for_extrapolation, Ps_for_extrapolation, βs_for_extrapolation, 
#                                                                                           Henry_CO2, Henry_CO2_err, 
#                                                                                           Henry_N2, Henry_N2_err) 

# 	# q_CO2, q_CO2_err = IntrinsicDACCycle.qₐ∞(βs_for_extrapolation, Kh_CO₂) #kJ/mol of gas
# 	# q_N2, q_N2_err = IntrinsicDACCycle.qₐ∞(βs_for_extrapolation, Kh_N₂) #kJ/mol of gas
# end
	

# ╔═╡ 33c60c80-bdb8-490d-8521-f48a4a0e964b
# remaining_files[2]

# ╔═╡ 371db837-bd6f-4fe7-83de-11d60794fed7
# begin
# 	plot(Ts_for_extrapolation, reshape(Henry_CO2_err, :))
# end

# ╔═╡ 436ae9a5-075c-48fb-9540-e5b35036c8e5
# begin
# 	plot(Ts_for_extrapolation, reshape(Henry_CO2, :))
# 	plot!(Ts_for_extrapolation, Henry_CO2 + 2 .* Henry_CO2_err,
# 		fillrange = Henry_CO2 - 2 .* Henry_CO2_err,
# 		fillalpha = 0.35, c = 1)

# 	plot!(xlabel = "Temperature (K)", 
# 	  ylabel ="Kh (mmol/(kg Pa))")
# end

# ╔═╡ 05346ffa-f0e8-45e5-b02e-875f9b02121c
# begin
# 	q_CO2, q_CO2_err = qₐ∞2(βs_for_extrapolation, Kh_CO₂) #kJ/mol of gas
# 	q_N2, q_N2_err = qₐ∞2(βs_for_extrapolation, Kh_N₂) #kJ/mol of gas
# end

# ╔═╡ b5db5f9e-f0c9-4c94-9134-f04ac3d5a81d
# plot(Ts_for_extrapolation, q_CO2)

# ╔═╡ e5a460c5-dcd3-417f-8371-be2790ff54ba
# minimum(q_CO2)

# ╔═╡ f3e4e84d-81ff-4690-ab90-b2cfdcda1705
# begin
# 	plot(Ts_for_extrapolation, real.(q_CO2_err), label= "Real part")
# 	plot!(Ts_for_extrapolation, imag.(q_CO2_err), label= "Imag. part")
# 	plot!(xlabel = "Temperature (K)", 
# 		  ylabel ="σ_Q (kJ/mol)")
# end

# ╔═╡ 8a7e11cc-4955-4f3e-a6ad-025bfb5e4040
# begin
# 	plot(Ts_for_extrapolation, q_CO2)
# 	plot!(Ts_for_extrapolation, q_CO2 + 2 .* real.(q_CO2_err),
# 			fillrange = q_CO2 - 2 .* real.(q_CO2_err),
# 			fillalpha = 0.35, c = 1)

# 	plot!(xlabel = "Temperature (K)", 
# 		  ylabel ="Q (kJ/mol)")
# end

# ╔═╡ c6cf48cb-350c-43d6-a41a-e9f4f779f7b7
# begin
# 	plot(Ts_for_extrapolation, q_CO2)
# 	plot!(Ts_for_extrapolation, q_CO2 + 2 .* real.(q_CO2_err),
# 			fillrange = q_CO2 - 2 .* real.(q_CO2_err),
# 			fillalpha = 0.35, c = 1)
# 	plot!(Ts_for_extrapolation, q_CO2 + 2 .* imag.(q_CO2_err),
# 			fillrange = q_CO2 - 2 .* imag.(q_CO2_err),
# 			fillalpha = 0.35, c = 2)

# 	plot!(xlabel = "Temperature (K)", 
# 		  ylabel ="Q (kJ/mol)",
# 		  xlims = (433, 435),
# 		  ylims = (75, 80)
# 	)
# end

# ╔═╡ 22a56b2e-39b5-444b-bf62-85308c4ef04d
# minimum(minimum(Kh_CO₂["Kcoeff_covar"]))

# ╔═╡ Cell order:
# ╠═fda35b50-71d7-11ef-23da-850f3a6532ea
# ╠═63ff5602-89e2-45ee-925c-543db2bc8e3e
# ╠═b8f5f181-367d-41ca-b463-cc3dc8f3b3b6
# ╠═bca614e9-ef9d-4eef-9e34-7b7f26c7b6cf
# ╠═d866129f-b8d0-4996-ad83-c55569749fc8
# ╠═d139dcb6-66a4-402b-b46c-8e519909fe0d
# ╠═7ce56fb6-a859-4678-b1ae-57a179e648cd
# ╠═bb9814ac-5a96-4022-9011-22afee662fdb
# ╠═c8d4c5ec-5eb1-4856-9590-01efaeb97f98
# ╠═10515cb5-4c36-455c-9f5a-e5db0ed82e63
# ╠═4fa763cd-df93-4e58-845c-0a98d45ac858
# ╠═db3c9b0c-da6c-4fdb-9813-ac949c1b096c
# ╠═01cf8b49-387e-46aa-bf6c-9a6d828f33a5
# ╠═70230652-6811-4911-ba9b-0d0113144f21
# ╠═3517897b-dc70-4db4-8e39-4e1d471fda74
# ╠═dfe4e04f-7b76-4691-abc6-bf9abc0f096c
# ╠═525abae2-dd07-48be-8aaf-ebb37389a628
# ╠═1084bedb-a0f2-4c6e-a1ec-330bdb339613
# ╠═92718c52-6d3d-49a9-b0a2-94c18655200b
# ╠═279e087b-f1ca-4bdd-b1aa-4ee7121ac751
# ╠═07d997f5-0758-48c7-b19e-40c57591bd21
# ╠═ea14a359-fd6d-4401-9b68-4634fec3582b
# ╠═9bed4fa9-df5b-4c18-bfd4-d027c4c46a64
# ╠═9ea586d8-9e96-4f00-a4ed-005b898e6063
# ╠═418fa3b5-58ce-42d8-aa18-e4e009b1af21
# ╠═93b65e57-c80a-49e0-ab4e-56a29010d1d4
# ╠═6f931e86-767c-4113-a30a-47dd81af0159
# ╠═546d3951-ce1a-434a-908a-a6930cb39482
# ╠═96fe7099-2607-45b4-8cd0-b84cfd33a38d
# ╠═46bdd140-0839-4133-b87d-28fecf144ec4
# ╠═08f327a4-1f42-4a1a-b71e-0597142f1382
# ╠═3a767149-49e6-4c1c-af7d-34c230563ce8
# ╠═7f5ec415-73d4-49db-a22e-cf3121134767
# ╠═f155ff38-7c69-41b7-8e29-2f82f043939a
# ╠═b0bfa892-4f31-498f-b272-c07587544489
# ╠═fdc614e6-f1c4-41a2-a8b0-2d1342a05646
# ╠═d0ddc413-6e2b-4691-ac7b-0f2a050160a9
# ╠═e86cc213-7bdb-4dc0-a5ff-76bbbc3da97d
# ╠═bd59fe52-8dc8-4674-bdbf-dabdcabe43df
# ╠═33c60c80-bdb8-490d-8521-f48a4a0e964b
# ╠═371db837-bd6f-4fe7-83de-11d60794fed7
# ╠═436ae9a5-075c-48fb-9540-e5b35036c8e5
# ╠═05346ffa-f0e8-45e5-b02e-875f9b02121c
# ╠═b5db5f9e-f0c9-4c94-9134-f04ac3d5a81d
# ╠═e5a460c5-dcd3-417f-8371-be2790ff54ba
# ╠═f3e4e84d-81ff-4690-ab90-b2cfdcda1705
# ╠═8a7e11cc-4955-4f3e-a6ad-025bfb5e4040
# ╠═c6cf48cb-350c-43d6-a41a-e9f4f779f7b7
# ╠═22a56b2e-39b5-444b-bf62-85308c4ef04d
