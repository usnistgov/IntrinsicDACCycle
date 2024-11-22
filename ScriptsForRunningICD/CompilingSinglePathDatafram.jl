### A Pluto.jl notebook ###
# v0.19.47

using Markdown
using InteractiveUtils

# ╔═╡ a75dac14-a2ce-11ef-2ab1-d7e16a195dbf
begin 
	using Pkg
	Pkg.activate("/users/asm6/Julia_scripts/IntrinsicDACCycle/")
end


# ╔═╡ 910cd017-f0fc-41e8-b7a3-b16e6783b191
begin
	using JSON
	using DataFrames
end
	

# ╔═╡ 9b96b016-c57f-47c9-bdbf-de5cf70a13c8
using IntrinsicDACCycle

# ╔═╡ 7cc7e6b7-05bd-490e-a905-23ce4fc7e1cf
using Statistics

# ╔═╡ 3ba2c4db-bf68-4ea3-9909-32b2939249ab
begin
	using Plots
end

# ╔═╡ 3dabbfbd-7923-4891-b4a3-5f754c69963e
begin
	Base_directory = "/users/asm6/DAC_data/"
	#get all completed 
	list_of_material_files = filter(x -> occursin.(".json",x), readdir(Base_directory*"SinglePath_Intrinsic_Cycle/"))
	#strip off the .json tag
    list_of_materials = replace.(list_of_material_files, ".json" => "")
	#strip off the filename prefix
	list_of_materials = replace.(list_of_materials, "SinglePathIDC_" => "")
end


# ╔═╡ 321c35d2-3932-4363-9231-6299bea16355
begin
	Non_adsorbing = []
	Not_close_enough = []
	truncated_to_nothing = []
	optimizer_failed = []
	materials_that_completed = []

	
	for name in list_of_materials
		@show name
		cycle_dict = JSON.parsefile(Base_directory*"SinglePath_Intrinsic_Cycle/SinglePathIDC_"*name*".json")

		cycle_df = DataFrame(cycle_dict)

		#if there is no key "Refresh_Path"
		if ~ ("Refresh_Path" in names(cycle_df))
			push!(truncated_to_nothing, name)
		else
			#If the purity is not a float64
			if any(typeof.(cycle_dict["Purity_captured_CO2"]) .!= Float64)
				push!(optimizer_failed, name)
			else
				push!(materials_that_completed, name)
			end
		end
				
	end
end
	

# ╔═╡ 5e4351fb-0e43-466b-b157-e5b27d847340
begin
	@show length(Non_adsorbing)
	@show length(Not_close_enough)
	@show length(truncated_to_nothing)
	@show length(optimizer_failed)
	@show length(materials_that_completed)
end

# ╔═╡ 1f02d66f-a4cd-47bf-8c45-103d0274a148
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

	#For each file that completed CALCULATE THE UNCERTAINTY IN CAPTURED CO2 AND N2
	for name in materials_that_completed

		results = JSON.parsefile(Base_directory*"SinglePath_Intrinsic_Cycle/SinglePathIDC_"*name*".json")
		results_df = DataFrame(results)

		#if there is no key "Captured_CO2_std"
		if ~ ("Captured_CO2_std" in names(results_df))
			objectives = IntrinsicDACCycle.Intrinisic_refresh_objectives_posterior_dist(directory, name, 
										   path1_T, path1_P, 
										   α, 100)

		
		
			if "Refresh_Path" in names(results_df)
				results["Intrinsic_capture_efficiency"] = mean(objectives[1])
		    	results["Intrinsic_capture_efficiency_std"] = std(objectives[1])
		    	results["Purity_captured_CO2"] = mean(objectives[2])
		    	results["Purity_captured_CO2_std"] = std(objectives[2])
				
			    results["Captured_CO2"] = mean(objectives[3])
			    results["Captured_CO2_std"] = std(objectives[3])
			    results["Captured_N2"] = mean(objectives[4])
			    results["Captured_N2_std"] = std(objectives[4])
			end

			#Write the file to disc
			filename = directory*"/SinglePath_Intrinsic_Cycle/SinglePathIDC_"*name*".json"
			open(filename,"w") do f
				JSON.print(f, results, 4)
			end
		end
	end
end

# ╔═╡ 477883e2-f6de-4b9f-88c7-106c242a3ba0
begin
	SinglePath_dataframe = DataFrame(Name =[], 
						  Intrinsic_capture_efficiency = [],
						  Intrinsic_capture_efficiency_std = [],
						  Purity_Captured_CO2 = [],
						  Purity_Captured_CO2_std = [],
						  Captured_CO2 = [],
						  Captured_CO2_std = [],
						  Captured_N2 = [],
						  Captured_N2_std = [])
	
	for name in materials_that_completed
		results_dict = JSON.parsefile(Base_directory*"SinglePath_Intrinsic_Cycle/SinglePathIDC_"*name*".json")
		Captured_CO2 = results_dict["Captured_CO2"]
		Captured_CO2_std = results_dict["Captured_CO2_std"]
		Captured_N2 = results_dict["Captured_N2"]
		Captured_N2_std = results_dict["Captured_N2_std"]
	
		Purity = results_dict["Purity_captured_CO2"]
		Purity_std = results_dict["Purity_captured_CO2_std"]
	
		Intrinsic_Capture_efficiency = results_dict["Intrinsic_capture_efficiency"]
		Intrinsic_Capture_efficiency_std = results_dict["Intrinsic_capture_efficiency_std"]
	
	
		for (ICE, 
			 ICE_std, 
			 Pure, 
			 Pure_std, 
			 CO2, 
			 CO2_std, 
			 N2, 
			 N2_std) in zip(Intrinsic_Capture_efficiency,
							   Intrinsic_Capture_efficiency_std,
							   Purity,
							   Purity_std,
							   Captured_CO2,
							   Captured_CO2_std,
							   Captured_N2,
							   Captured_N2_std)
												 
			push!(SinglePath_dataframe, (name,
								 ICE, 
								 ICE_std, 
								 Pure, 
								 Pure_std, 
								 CO2, 
								 CO2_std, 
								 N2, 
								 N2_std))
		end
	end
	
end

# ╔═╡ 297c501e-065c-4287-bca1-aeb51f7b27f0
SinglePath_dataframe

# ╔═╡ 1a3c7613-9cfb-434e-8ab7-e21564e2b3c4
begin
	histogram2d(SinglePath_dataframe[:,"Purity_Captured_CO2"], 
				SinglePath_dataframe[:,"Intrinsic_capture_efficiency"],
				bins = (100,100),
				# show_empty_bins =true,
				normalize=:none,
				# color = cgrad(:plasma,scale=log))
				color = :plasma,
				colorbar_title = "Counts"
				)
	xlabel!("Purity of Captured CO2")
	ylabel!("Intrinsic Capture Efficency (mol/J)")
end

# ╔═╡ 09d8af7d-4f12-426f-803a-0df40282fc61
begin
	SinglePath_pareto = DataFrame(Name =[], 
					  Intrinsic_capture_efficiency = [],
					  Intrinsic_capture_efficiency_std = [],
					  Purity_Captured_CO2 = [],
					  Purity_Captured_CO2_std = [],
					  Captured_CO2 = [],
					  Captured_CO2_std = [],
					  Captured_N2 = [],
					  Captured_N2_std = [])

	sort!(SinglePath_dataframe, [order("Intrinsic_capture_efficiency", rev = true), order("Purity_Captured_CO2", rev = true)])

	push!(SinglePath_pareto, SinglePath_dataframe[1,:])

	for entry in eachrow(SinglePath_dataframe[2:end, :])
		if entry["Purity_Captured_CO2"] >= SinglePath_pareto[end, "Purity_Captured_CO2"]
			push!(SinglePath_pareto, entry)
		end
	end
end

# ╔═╡ b8efcd51-877b-4ba2-bc06-92ebb75b082f
SinglePath_pareto

# ╔═╡ 91a3d227-7163-44b5-af88-8aebd60f4b43
begin
	#Write the file to disc
	filename = "/users/asm6/Julia_scripts/SinglePath_IDC_pareto.json"
	open(filename,"w") do f
		JSON.print(f, SinglePath_pareto, 4)
	end
end

# ╔═╡ f70f4153-4f38-47b3-aad2-e4a27662d914
begin
	# scatter(Opt_dataframe.Purity_Captured_CO2), 					  Opt_dataframe.Intrinsic_capture_efficiency)

	scatter(SinglePath_pareto[:,"Purity_Captured_CO2"], 
			SinglePath_pareto[:,"Intrinsic_capture_efficiency"],
			xerror = SinglePath_pareto[:,"Purity_Captured_CO2_std"],
			yerror = SinglePath_pareto[:,"Intrinsic_capture_efficiency_std"],
			label = "Pareto Front")

	xlabel!("Purity Captured CO₂")
	ylabel!("Intrinsic Capture Efficiency (mol/J)")
	savefig("/users/asm6/Julia_scripts/SinglePath_pareto.svg")
end

# ╔═╡ 17e86dd1-3c1a-45a4-bd04-b03d39f4fbb8
Base_directory

# ╔═╡ e2c6ed6f-df62-49fd-afc6-7578f045cc22
begin
	plot1 = 	histogram2d(SinglePath_dataframe[:,"Purity_Captured_CO2"], 
				SinglePath_dataframe[:,"Intrinsic_capture_efficiency"],
				bins = (100,100),
				# show_empty_bins =true,
				normalize=:none,
				# color = cgrad(:plasma,scale=log))
				color = :plasma,
				colorbar_title = "Counts"
				)
	plot1 = xlabel!("Purity of Captured CO₂")
	plot1 = ylabel!("Intrinsic Capture Efficency (mol/J)")	
	plot1 = title!("a)", titlelocation = :left)

	plot2 = scatter(SinglePath_pareto[:,"Purity_Captured_CO2"], 
			SinglePath_pareto[:,"Intrinsic_capture_efficiency"],
			xerror = SinglePath_pareto[:,"Purity_Captured_CO2_std"],
			yerror = SinglePath_pareto[:,"Intrinsic_capture_efficiency_std"],
			label = "Pareto Front")

	plot2 = xlabel!("Purity Captured CO₂")
	plot2 = ylabel!("Intrinsic Capture Efficiency (mol/J)")
	plot2 = title!("b)", titlelocation = :left)

	plot(plot1, plot2, layout = (1,2))
	plot!(size=(900,400), left_margin = 15*Plots.mm, bottom_margin = 10*Plots.mm)
	savefig("/users/asm6/Julia_scripts/SinglePath_performance.svg")
end

# ╔═╡ Cell order:
# ╠═a75dac14-a2ce-11ef-2ab1-d7e16a195dbf
# ╠═910cd017-f0fc-41e8-b7a3-b16e6783b191
# ╠═9b96b016-c57f-47c9-bdbf-de5cf70a13c8
# ╠═7cc7e6b7-05bd-490e-a905-23ce4fc7e1cf
# ╠═3ba2c4db-bf68-4ea3-9909-32b2939249ab
# ╠═3dabbfbd-7923-4891-b4a3-5f754c69963e
# ╠═321c35d2-3932-4363-9231-6299bea16355
# ╠═5e4351fb-0e43-466b-b157-e5b27d847340
# ╠═1f02d66f-a4cd-47bf-8c45-103d0274a148
# ╠═477883e2-f6de-4b9f-88c7-106c242a3ba0
# ╠═297c501e-065c-4287-bca1-aeb51f7b27f0
# ╠═1a3c7613-9cfb-434e-8ab7-e21564e2b3c4
# ╠═09d8af7d-4f12-426f-803a-0df40282fc61
# ╠═b8efcd51-877b-4ba2-bc06-92ebb75b082f
# ╠═91a3d227-7163-44b5-af88-8aebd60f4b43
# ╠═f70f4153-4f38-47b3-aad2-e4a27662d914
# ╠═17e86dd1-3c1a-45a4-bd04-b03d39f4fbb8
# ╠═e2c6ed6f-df62-49fd-afc6-7578f045cc22
