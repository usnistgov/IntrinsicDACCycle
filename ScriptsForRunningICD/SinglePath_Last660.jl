### A Pluto.jl notebook ###
# v0.19.47

using Markdown
using InteractiveUtils

# ╔═╡ 7207e120-a5bd-11ef-3d58-f515a2436be8
using Pkg

# ╔═╡ 8ac3c33e-6f74-4bfe-afd3-53d4fa0c3f23
Pkg.activate("/users/asm6/Julia_scripts/IntrinsicDACCycle")


# ╔═╡ 94d20fb8-5cde-44f3-98c1-7dc7c8a0d07c
begin
	using Plots
	using Distributed
	using JSON
	using DataFrames

end


# ╔═╡ ebfa4aa1-d7bd-40fe-8d04-729451eeca2e
using LinearAlgebra

# ╔═╡ 78cff125-3b07-413a-b81e-9ecadeba41e5
using NaNStatistics

# ╔═╡ 2f76bcd2-5f12-4e1a-97ae-eb8beff103b4
using Statistics

# ╔═╡ d776a7a1-58b8-4e27-b79a-28ef7e1870b6
using Metaheuristics 

# ╔═╡ 9140caab-cdfb-4b16-aa7b-5942f8d9ebd5
using IntrinsicDACCycle

# ╔═╡ 5799d8c0-254b-4dbb-a499-1c0553878d6c
Base_directory = "/users/asm6/DAC_data"

# ╔═╡ 85763539-6c03-486a-b0b3-d7aeba555cda
begin
	#get all the material files
	list_of_material_files = filter(x -> occursin.(".json",x), readdir(Base_directory*"/CSD_FEASST_Materials/Materials/"))
	#strip off the .json tag
    list_of_materials = replace.(list_of_material_files, ".json" => "")
	#filter for _clean matierals
	# list_of_clean_materials = filter(x -> occursin.("_clean", x), list_of_materials)
end


# ╔═╡ c0ff92b9-ef15-40af-8bda-e502d12c796d
length(list_of_materials)

# ╔═╡ f9fc4d9f-f2c2-46d9-9222-7ff2d62ec8e9
begin
	#get all the Optimized Intrinsic Cycle files
	list_of_SinglePath_files = filter(x -> occursin.(".json",x), readdir(Base_directory*"/SinglePath_Intrinsic_Cycle/"))
	#strip off the .json tag
    list_of_SinglePath = replace.(list_of_SinglePath_files, ".json" => "")
	list_of_SinglePath = replace.(list_of_SinglePath, "SinglePathIDC_" => "")
end

# ╔═╡ a35c47da-1038-4237-9385-d398fafb9b84
# begin
# 	for test in list_of_SinglePath
# 		dict = JSON.parsefile(Base_directory*"/SinglePath_Intrinsic_Cycle/SinglePathIDC_"*test*".json")
	
# 		oh_no = dict["Intrinsic_capture_efficiency"] == dict["Purity_captured_CO2"]
# 		if oh_no
# 			@show "Removing:"
# 			@show "test"
# 			rm(Base_directory*"/SinglePath_Intrinsic_Cycle/SinglePathIDC_"*test*".json")
# 		end
# 	end
# end

		

# ╔═╡ c07f8da9-4b11-4957-ba59-0e52953d5cdc
length(list_of_SinglePath)

# ╔═╡ 4dce02f3-43e1-4779-867c-ff9e78528e6e
begin
	remaining_files = setdiff(list_of_materials, list_of_SinglePath)
end

# ╔═╡ f841825d-2cc3-4f28-89d7-95275ec87346
function Run_a_SinglePath_IntrinsicDACCycle(name)
    T_start = 250.0 #K
    T_end = 375.0 #K
    P_start = 101325.0 #Pa
    P_end = 0.4 * P_start	

    steps = 500
    #Generate Path 1: Change T, then P
    path1_T_A = LinRange(T_start, T_end, steps)
    path1_T_B = LinRange(T_end, T_end, steps)
    path1_T = vcat(path1_T_A, path1_T_B)

    path1_P_A = LinRange(P_start, P_start, steps)
    path1_P_B = LinRange(P_start, P_end, steps)
    path1_P = vcat(path1_P_A, path1_P_B) 

    directory = "/users/asm6/DAC_data/"
    alpha = 400/1000000
    results = IntrinsicDACCycle.Intrinisic_refresh_path(directory, name, 
						 path1_T, path1_P, 
					  	 alpha)

    objectives = IntrinsicDACCycle.Intrinisic_refresh_objectives_posterior_dist(directory, name, 
								       path1_T, path1_P, 
								       alpha, 100)
    results["Intrinsic_capture_efficiency"] = mean(objectives[1])
    results["Intrinsic_capture_efficiency_std"] = std(objectives[1])
    results["Purity_captured_CO2"] = mean(objectives[1])
    results["Purity_captured_CO2_std"] = std(objectives[1])


    filename = directory*"/SinglePath_Intrinsic_Cycle/SinglePathIDC_"*name*".json"
    open(filename,"w") do f
        JSON.print(f, results, 4)
    end
end


# ╔═╡ cfc7793c-bbf8-40ad-bd2f-1a275bd9b565
for file in remaining_files
	@show "running", file
	try
		Run_an_Opitimize_IntrinsicDACCycle(file)
		@show "done!"
	catch
		@show "Skipping"
	end
	
end

# ╔═╡ 2a634690-a7d3-4a9b-b50b-bf6d77ea24f3
begin
	directory = "/users/asm6/DAC_data/"
    α = 400/1000000

    T_start = 250.0 #K
    T_end = 375.0 #K
    P_start = 101325.0 #Pa
    P_end = 0.4 * P_start	

    steps = 500
    #Generate Path 1: Change T, then P
    path1_T_A = LinRange(T_start, T_end, steps)
    path1_T_B = LinRange(T_end, T_end, steps)
    path1_T = vcat(path1_T_A, path1_T_B)

    path1_P_A = LinRange(P_start, P_start, steps)
    path1_P_B = LinRange(P_start, P_end, steps)
    path1_P = vcat(path1_P_A, path1_P_B) 
	
	#For each remaining file
	for name in remaining_files
		# complete_flag = 0
	# name = remaining_files[1]

		results = IntrinsicDACCycle.Intrinisic_refresh_path(directory, name, 
						 path1_T, path1_P, 
						 α)
	
		objectives = IntrinsicDACCycle.Intrinisic_refresh_objectives_posterior_dist(directory, name, 
									   path1_T, path1_P, 
									   α, 100)
	
		if "Refresh_Path" in results.keys
		    results["Intrinsic_capture_efficiency"] = mean(objectives[1])
		    results["Intrinsic_capture_efficiency_std"] = std(objectives[1])
		    results["Purity_captured_CO2"] = mean(objectives[2])
		    results["Purity_captured_CO2_std"] = std(objectives[2])
		end
	
		# end
	
			
		# 	#try 10 times 
		# 	for i in 1:10 
		# 		@show complete_flag
				
		# 		if complete_flag <= 10
		# 			#calculate the optimized cycle
		# 			global results = IntrinsicDACCycle.Optimize_Intrinsic_Refresh_w_err(directory, name, α)
		# 			#iterate the flag
		# 			complete_flag += 1
		# 			#flag if a good calc
		# 			if isassigned(results["Captured_CO2"])
		# 				complete_flag = 1000
		# 			end
		# 		end
		# 	end
			
		# 	#Clean up the undefined results	
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
			
		# 	if ~isassigned(results["Purity_captured_CO2_std"])
		# 		results["Purity_captured_CO2_std"] = []
		# 	end
		
		# 	if ~isassigned(results["Captured_CO2"])
		# 		results["Refresh_Path"] = []
		# 	end
		
		# 	if ~isassigned(results["Captured_CO2"])
		# 		results["E_Balance"] = []
		# 	end
		
		#Write the file to disc
		filename = directory*"/SinglePath_Intrinsic_Cycle/SinglePathIDC_"*name*".json"
		open(filename,"w") do f
			JSON.print(f, results, 4)
		end
	end
end
		

# ╔═╡ 275b28d1-9437-4e9e-a999-f02efb2028a2
results

# ╔═╡ 786ec015-0828-49b1-9571-ffb4b782c8e5
~ ("Refresh_Path" in results.keys )

# ╔═╡ Cell order:
# ╠═7207e120-a5bd-11ef-3d58-f515a2436be8
# ╠═8ac3c33e-6f74-4bfe-afd3-53d4fa0c3f23
# ╠═94d20fb8-5cde-44f3-98c1-7dc7c8a0d07c
# ╠═ebfa4aa1-d7bd-40fe-8d04-729451eeca2e
# ╠═78cff125-3b07-413a-b81e-9ecadeba41e5
# ╠═2f76bcd2-5f12-4e1a-97ae-eb8beff103b4
# ╠═d776a7a1-58b8-4e27-b79a-28ef7e1870b6
# ╠═9140caab-cdfb-4b16-aa7b-5942f8d9ebd5
# ╠═5799d8c0-254b-4dbb-a499-1c0553878d6c
# ╠═85763539-6c03-486a-b0b3-d7aeba555cda
# ╠═c0ff92b9-ef15-40af-8bda-e502d12c796d
# ╠═f9fc4d9f-f2c2-46d9-9222-7ff2d62ec8e9
# ╠═a35c47da-1038-4237-9385-d398fafb9b84
# ╠═c07f8da9-4b11-4957-ba59-0e52953d5cdc
# ╠═4dce02f3-43e1-4779-867c-ff9e78528e6e
# ╠═f841825d-2cc3-4f28-89d7-95275ec87346
# ╠═cfc7793c-bbf8-40ad-bd2f-1a275bd9b565
# ╠═2a634690-a7d3-4a9b-b50b-bf6d77ea24f3
# ╠═275b28d1-9437-4e9e-a999-f02efb2028a2
# ╠═786ec015-0828-49b1-9571-ffb4b782c8e5
