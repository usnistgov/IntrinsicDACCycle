### A Pluto.jl notebook ###
# v0.19.47

using Markdown
using InteractiveUtils

# ╔═╡ e0d95bff-efe1-4e6a-be09-cd54644d9df2
begin 
	using Pkg
	Pkg.activate("/users/asm6/Julia_scripts/IntrinsicDACCycle/")
end


# ╔═╡ 10d53dca-9147-11ef-3a6d-951ef039c737
begin
	using JSON
	using DataFrames
end
	

# ╔═╡ 420a9a15-2ee7-40a3-b14e-33384165bc40
begin
	using Plots
end

# ╔═╡ ae79d9cd-f16a-4275-a15a-08e710ea8242
begin
	Base_directory = "/users/asm6/DAC_data/"
	#get all completed 
	list_of_material_files = filter(x -> occursin.(".json",x), readdir(Base_directory*"Optimized_Intrinsic_Cycle/"))
	#strip off the .json tag
    list_of_materials = replace.(list_of_material_files, ".json" => "")
	#strip off the filename prefix
	list_of_materials = replace.(list_of_materials, "OptIDC_" => "")
end


# ╔═╡ a0cb4aa3-6180-44f4-b250-ad1ed3603b70
begin
	Non_adsorbing = []
	Not_close_enough = []
	truncated_to_nothing = []
	optimizer_failed = []
	materials_that_optimized = []
	
	for name in list_of_materials
		
		try
			Opt_cycle_dict = JSON.parsefile(Base_directory*"Optimized_Intrinsic_Cycle/OptIDC_"*name*".json")
			
			if Opt_cycle_dict["Saturation Adsorb Any CO2"] == false
				push!(Non_adsorbing, name)
			elseif Opt_cycle_dict["Close Enough to Linear Adsorption"] == false
				push!(Not_close_enough, name)
			elseif length(Opt_cycle_dict["Captured_CO2"]) == 0
				push!(truncated_to_nothing, name)
			elseif any(typeof.(Opt_cycle_dict["Purity_captured_CO2"]) .!= Float64)
				push!(optimizer_failed, name)
			else
				push!(materials_that_optimized, name)
			end
		catch
			# @show "Removing:"
			# @show name
			# rm(Base_directory*"Optimized_Intrinsic_Cycle/OptIDC_"*name*".json")
		end
	
	end
	
end

# ╔═╡ 54f22445-f4bb-46eb-ae95-3a8aee463a8a


# ╔═╡ a2cb1ac0-289b-4008-ade8-d15720ce3baa
begin
	#Find the materials where the extraploation of Kh failed - likely due to insufficient statistics (even after all the GCMC trials)
	list_of_all_material_files = filter(x -> occursin.(".json",x), readdir(Base_directory*"CSD_FEASST_Materials/Materials/"))
	#strip off the .json tag
    list_of_all_materials = replace.(list_of_all_material_files, ".json" => "")
	
	
	failed_Kh_extrapolation = collect(setdiff(Set(list_of_all_materials), Set(list_of_materials)))
	
	
end

# ╔═╡ 39580cbf-7519-4e6f-b9d8-9da07557e175
optimizer_failed

# ╔═╡ ba133f32-5d9e-4307-9c10-03beaf778c61
begin
	@show length(failed_Kh_extrapolation)
	@show length(truncated_to_nothing)
	@show length(Non_adsorbing)
	@show length(Not_close_enough)
	@show length(optimizer_failed)

	@show length(materials_that_optimized)

	@show sum([length(failed_Kh_extrapolation),
			length(truncated_to_nothing),
			length(Non_adsorbing),
			length(Not_close_enough),
			length(optimizer_failed),
			length(materials_that_optimized)])
end

# ╔═╡ cefa7523-49db-48bf-8f90-20befd802b9f
begin
	test1_name = materials_that_optimized[209]
	test2_name = materials_that_optimized[212]

	test1_file = JSON.parsefile(Base_directory*"Optimized_Intrinsic_Cycle/OptIDC_"*test1_name*".json")

	test2_file = JSON.parsefile(Base_directory*"Optimized_Intrinsic_Cycle/OptIDC_"*test2_name*".json")
end

# ╔═╡ 8a94b9aa-251a-43f8-83f9-ecfe4b6dca2e
test1_file

# ╔═╡ 6d7703aa-b090-4f7e-840a-93160cea2855
begin
	Opt_dataframe = DataFrame(Name =[], 
						  Intrinsic_capture_efficiency = [],
						  Intrinsic_capture_efficiency_std = [],
						  Purity_Captured_CO2 = [],
						  Purity_Captured_CO2_std = [],
						  Captured_CO2 = [],
						  Captured_CO2_std = [],
						  Captured_N2 = [],
						  Captured_N2_std = [],
						  Temp_start = [],
						  Temp_end = [],
						  Press_start = [],
					      Press_end = [])
	
	for name in materials_that_optimized
		results_dict = JSON.parsefile(Base_directory*"Optimized_Intrinsic_Cycle/OptIDC_"*name*".json")
		
		Captured_CO2 = results_dict["Captured_CO2"]
		Captured_CO2_std = results_dict["Captured_CO2_std"]
		Captured_N2 = results_dict["Captured_N2"]
		Captured_N2_std = results_dict["Captured_N2_std"]
	
		Purity = results_dict["Purity_captured_CO2"]
		Purity_std = results_dict["Purity_captured_CO2_std"]
	
		Intrinsic_Capture_efficiency = results_dict["Intrinsic_capture_efficiency"]
		Intrinsic_Capture_efficiency_std = results_dict["Intrinsic_capture_efficiency_std"]
		
		
		
		
		Temp_starts = []
		Temp_ends = []
		# minimum(test2_file["Refresh_Path"]["Temperatures"], dims = 2)
		for Temp in results_dict["Refresh_Path"]["Temperatures"]
			Temp_start = minimum(Temp)
			push!(Temp_starts, Temp_start)
	
			Temp_end = maximum(Temp)
			push!(Temp_ends, Temp_end)
		end
	
		Press_starts = []
		Press_ends = []
		# minimum(test2_file["Refresh_Path"]["Temperatures"], dims = 2)
		for Press in results_dict["Refresh_Path"]["Pressures"]
			
			Press_start = maximum(Press)
			push!(Press_starts, Press_start)
	
			Press_end = minimum(Press)
			push!(Press_ends, Press_end)
		end
	
		for (ICE, 
			 ICE_std, 
			 Pure, 
			 Pure_std, 
			 CO2, 
			 CO2_std, 
			 N2, 
			 N2_std, 
			 Temp_start, 
			 Temp_end, 
			 Press_start, 
			 Press_end) in zip(Intrinsic_Capture_efficiency,
							   Intrinsic_Capture_efficiency_std,
							   Purity,
							   Purity_std,
							   Captured_CO2,
							   Captured_CO2_std,
							   Captured_N2,
							   Captured_N2_std,
							   Temp_starts,
							   Temp_ends,
							   Press_starts,
							   Press_ends)
												 
			push!(Opt_dataframe, (name,
								 ICE, 
								 ICE_std, 
								 Pure, 
								 Pure_std, 
								 CO2, 
								 CO2_std, 
								 N2, 
								 N2_std, 
								 Temp_start, 
								 Temp_end, 
								 Press_start, 
								 Press_end))
		end
	end
	
end

# ╔═╡ 2b8b8219-dcbf-476f-9b8c-59bdf830ad86
Opt_dataframe

# ╔═╡ f6d2315c-a85d-4e4a-a1a8-5937dc2221ff
begin
	# scatter(Opt_dataframe.Purity_Captured_CO2), 					  Opt_dataframe.Intrinsic_capture_efficiency)

	scatter(Opt_dataframe[:,"Purity_Captured_CO2"], 
	Opt_dataframe[:,"Intrinsic_capture_efficiency"])
end

# ╔═╡ e2782249-ff89-40c5-8d0f-f33f4be960aa


# ╔═╡ 2ba736b8-7b72-4672-af7a-6f68db19e487
begin
	histogram2d(Opt_dataframe[:,"Purity_Captured_CO2"], 
				Opt_dataframe[:,"Intrinsic_capture_efficiency"],
				bins = (100,100),
				# show_empty_bins =true,
				normalize=:none,
				# color = cgrad(:plasma,scale=log))
				color = :plasma,
				# colorbar_title = "Counts"
				)
	xlabel!("Purity of Captured CO2")
	ylabel!("Intrinsic Capture Efficency (mol/J)")
end

# ╔═╡ b4d94a59-f02c-4f05-87a1-96af18c48de5
begin
	histogram2d(Opt_dataframe[:,"Temp_start"],
				Opt_dataframe[:,"Intrinsic_capture_efficiency"], 
				bins = (50,50),
				# show_empty_bins =true,
				normalize=:none,
				# color = cgrad(:plasma,scale=log))
				color = :plasma)
	xlabel!("Start Temperature (K)")
	ylabel!("Intrinsic Capture Efficency (mol/J)")
end

# ╔═╡ 443c7801-e66a-43af-bed5-48a450e8274c
begin
	histogram2d(Opt_dataframe[:,"Temp_end"],
				Opt_dataframe[:,"Intrinsic_capture_efficiency"], 
				bins = (50,50),
				# show_empty_bins =true,
				normalize=:none,
				# color = cgrad(:plasma,scale=log))
				color = :plasma)
	xlabel!("End Temperature (K)")
	ylabel!("Intrinsic Capture Efficency (mol/J)")
end

# ╔═╡ 62c74a4e-6dad-4b70-a5c2-50075a856f5e
begin
	histogram2d(Opt_dataframe[:,"Press_start"],
				Opt_dataframe[:,"Intrinsic_capture_efficiency"], 
				bins = (50,50),
				# show_empty_bins =true,
				normalize=:none,
				# color = cgrad(:plasma,scale=log))
				color = :plasma)
	xlabel!("Start Pressure (Pa)")
	ylabel!("Intrinsic Capture Efficency (mol/J)")
end

# ╔═╡ 37c70cbc-bad7-44d0-abe8-75ba6fb0410c
begin
	histogram2d(Opt_dataframe[:,"Press_end"],
				Opt_dataframe[:,"Intrinsic_capture_efficiency"], 
				bins = (50,50),
				# show_empty_bins =true,
				normalize=:none,
				# color = cgrad(:plasma,scale=log))
				color = :plasma)
	xlabel!("End Pressure (Pa)")
	ylabel!("Intrinsic Capture Efficency (mol/J)")
end

# ╔═╡ f0a5f4be-d152-4218-8e79-be63f2efbca8
begin
	histogram2d(Opt_dataframe[:,"Temp_start"],
				Opt_dataframe[:,"Purity_Captured_CO2"], 
				bins = (50,50),
				# show_empty_bins =true,
				normalize=:none,
				# color = cgrad(:plasma,scale=log))
				color = :plasma)
	xlabel!("Start Temperature (K)")
	ylabel!("Purity Captured CO2")
end

# ╔═╡ 973addd2-8c3e-408a-a6a7-90ef2d1be7bc
begin
	histogram2d(Opt_dataframe[:,"Temp_end"],
				Opt_dataframe[:,"Purity_Captured_CO2"], 
				bins = (50,50),
				# show_empty_bins =true,
				normalize=:none,
				# color = cgrad(:plasma,scale=log))
				color = :plasma)
	xlabel!("End Temperature (K)")
	ylabel!("Purity Captured CO2")
end

# ╔═╡ 3f3490b6-725d-4ca2-ad1c-518b72dae35b
begin
	histogram2d(Opt_dataframe[:,"Press_start"],
				Opt_dataframe[:,"Purity_Captured_CO2"], 
				bins = (50,50),
				# show_empty_bins =true,
				normalize=:none,
				# color = cgrad(:plasma,scale=log))
				color = :plasma)
	xlabel!("Start Pressure (Pa)")
	ylabel!("Purity Captured CO2")
end

# ╔═╡ b09ee4e9-6086-4ef1-8d93-04a03be09c6b
begin
	histogram2d(Opt_dataframe[:,"Press_end"],
				Opt_dataframe[:,"Purity_Captured_CO2"], 
				bins = (50,50),
				# show_empty_bins =true,
				normalize=:none,
				# color = cgrad(:plasma,scale=log))
				color = :plasma)
	xlabel!("End Pressure (Pa)")
	ylabel!("Purity Captured CO2")
end

# ╔═╡ 9ea02bcf-9e8e-4a7c-9d26-ff87dcbf3e06
begin
	histogram(Opt_dataframe[:,"Purity_Captured_CO2"], 
				# bins = (50,50),
				normalize=:none,)
	xlabel!("Purity Captured CO2")
	ylabel!("Counts")
end

# ╔═╡ 2005253e-4cc2-498d-86a4-2fde787a266c
begin
	histogram(Opt_dataframe[:,"Intrinsic_capture_efficiency"], 
				# bins = (50,50),
				normalize=:none,)
	xlabel!("Intrinsic Capture Efficiency (mol/J)")
	ylabel!("Counts")
end

# ╔═╡ acf56879-d065-4127-be70-fc4e69538345
begin
	histogram(Opt_dataframe[:,"Temp_start"], 
				# bins = (50,50),
				normalize=:none,
				label = :none)
	xlabel!("Start Temperature (K)")
	ylabel!("Counts")
	
	savefig("/users/asm6/Julia_scripts/OptIDC_TempHist.svg")
end

# ╔═╡ c804ea9d-b399-47f0-8dd5-d9388f6b6c70
begin
	histogram(Opt_dataframe[:,"Temp_end"], 
				# bins = (50,50),
				normalize=:none,)
	xlabel!("End Temperature (K)")
	ylabel!("Counts")
end

# ╔═╡ add66865-896b-4ab3-a45d-0e39186883d0
begin
	histogram(Opt_dataframe[:,"Press_start"], 
				# bins = (50,50),
				normalize=:none,)
	xlabel!("Start Pressure (Pa)")
	ylabel!("Counts")
end

# ╔═╡ 47d70bb2-ff12-49ae-8f85-50bb160b2fa4
begin
	histogram(Opt_dataframe[:,"Press_end"], 
				# bins = (50,50),
				normalize=:none,)
	xlabel!("End Pressure (Pa)")
	ylabel!("Counts")
end

# ╔═╡ 81e42773-9748-4dd9-91b4-f412abae1655
# begin 
# 	l = @layout grid(6,6)
# 	p1 = histogram2d(Opt_dataframe[:,"Purity_Captured_CO2"], 
# 				Opt_dataframe[:,"Intrinsic_capture_efficiency"],
# 				bins = (100,100),
# 				# show_empty_bins =true,
# 				normalize=:none,
# 				# color = cgrad(:plasma,scale=log))
# 				color = :plasma)

# 	p2 = histogram2d(Opt_dataframe[:,"Purity_Captured_CO2"], 
# 				Opt_dataframe[:,"Temp_start"],
# 				bins = (100,100),
# 				# show_empty_bins =true,
# 				normalize=:none,
# 				# color = cgrad(:plasma,scale=log))
# 				color = :plasma)
# 	p3 = histogram2d(Opt_dataframe[:,"Purity_Captured_CO2"], 
# 				Opt_dataframe[:,"Temp_end"],
# 				bins = (100,100),
# 				# show_empty_bins =true,
# 				normalize=:none,
# 				# color = cgrad(:plasma,scale=log))
# 				color = :plasma)
# 	p4 = histogram2d(Opt_dataframe[:,"Purity_Captured_CO2"], 
# 				Opt_dataframe[:,"Press_start"],
# 				bins = (100,100),
# 				# show_empty_bins =true,
# 				normalize=:none,
# 				# color = cgrad(:plasma,scale=log))
# 				color = :plasma)
# 	p5 = histogram2d(Opt_dataframe[:,"Purity_Captured_CO2"], 
# 				Opt_dataframe[:,"Press_end"],
# 				bins = (100,100),
# 				# show_empty_bins =true,
# 				normalize=:none,
# 				# color = cgrad(:plasma,scale=log))
# 				color = :plasma)

# 	p6 = histogram2d(Opt_dataframe[:,"Intrinsic_capture_efficiency"], 
# 				Opt_dataframe[:,"Temp_start"],
# 				bins = (100,100),
# 				# show_empty_bins =true,
# 				normalize=:none,
# 				# color = cgrad(:plasma,scale=log))
# 				color = :plasma)
# 	p7 = histogram2d(Opt_dataframe[:,"Intrinsic_capture_efficiency"], 
# 				Opt_dataframe[:,"Temp_end"],
# 				bins = (100,100),
# 				# show_empty_bins =true,
# 				normalize=:none,
# 				# color = cgrad(:plasma,scale=log))
# 				color = :plasma)
# 	p8 = histogram2d(Opt_dataframe[:,"Intrinsic_capture_efficiency"], 
# 				Opt_dataframe[:,"Press_start"],
# 				bins = (100,100),
# 				# show_empty_bins =true,
# 				normalize=:none,
# 				# color = cgrad(:plasma,scale=log))
# 				color = :plasma)
# 	p9 = histogram2d(Opt_dataframe[:,"Intrinsic_capture_efficiency"], 
# 				Opt_dataframe[:,"Press_end"],
# 				bins = (100,100),
# 				# show_empty_bins =true,
# 				normalize=:none,
# 				# color = cgrad(:plasma,scale=log))
# 				color = :plasma)

# 	plot(_, p1, p2, p3, p4, p5, 
# 		 _, _, p6, p7, p8, p9,
# 		_, _, _, _, _, _, 
# 		_, _, _, _, _, _, 
# 		_, _, _, _, _, _, 
# 		_, _, _, _, _, _, 
# 		layout = l)
# end

# ╔═╡ eb926db5-b130-4c3b-888e-a8fc542b1419


# ╔═╡ 2d3558c7-635d-4acd-9edc-199275ea97cf
begin
	# w = log10.(Opt_dataframe[:,"Purity_Captured_CO2"])
	histogram2d(log10.(Opt_dataframe[:,"Purity_Captured_CO2"]), 
				log10.(1 ./(Opt_dataframe[:,"Intrinsic_capture_efficiency"])),
				bins = (50,100),
				# weights = w,
				# show_empty_bins =true,
				normalize=:none,
				# color = cgrad(:plasma,scale=log))
				color = :plasma)
	
end

# ╔═╡ 1a010f1a-240f-41b9-8c02-0371acac3360
begin
	sort!(Opt_dataframe, [order("Intrinsic_capture_efficiency", rev = false), order("Purity_Captured_CO2", rev = false)])
	scatter(Opt_dataframe[:,"Purity_Captured_CO2"], 
			Opt_dataframe[:,"Intrinsic_capture_efficiency"],
			markerstrokewidth=0,
			# xerror = Opt_dataframe[:,"Purity_Captured_CO2_std"],
			# yerror = Opt_dataframe[:,"Intrinsic_capture_efficiency_std"],
			zcolor = Opt_dataframe[:,"Temp_start"],
			color =:coolwarm,
			colorbar_title = "Start Temperature (K)")
	title!("Start Temperature trend")
	xlabel!("Purity Captured CO2")
	ylabel!("Intrinsic Capture Efficiency (mol/J)")
end

# ╔═╡ 7e27220c-12c9-4fdb-a877-526308708de5


# ╔═╡ f8fdb476-5d3d-4706-952a-3a07c8129833


# ╔═╡ 9226b1f0-de3b-4fcb-a9e3-ec4423815e56
begin
	Opt_pareto = DataFrame(Name =[], 
					  Intrinsic_capture_efficiency = [],
					  Intrinsic_capture_efficiency_std = [],
					  Purity_Captured_CO2 = [],
					  Purity_Captured_CO2_std = [],
					  Captured_CO2 = [],
					  Captured_CO2_std = [],
					  Captured_N2 = [],
					  Captured_N2_std = [],
					  Temp_start = [],
					  Temp_end = [],
					  Press_start = [],
					  Press_end = [])

	sort!(Opt_dataframe, [order("Intrinsic_capture_efficiency", rev = true), order("Purity_Captured_CO2", rev = true)])

	push!(Opt_pareto, Opt_dataframe[1,:])

	for entry in eachrow(Opt_dataframe[2:end, :])
		if entry["Purity_Captured_CO2"] >= Opt_pareto[end, "Purity_Captured_CO2"]
			push!(Opt_pareto, entry)
		end
	end
end

# ╔═╡ 4f3441ac-97ad-451f-8d84-16e20fe3f79a
Opt_pareto

# ╔═╡ 095c2897-0461-4a90-b051-10f2ee360dd6
begin
	#Write the file to disc
	filename = "/users/asm6/Julia_scripts/Optimized_IDC_pareto.json"
	open(filename,"w") do f
		JSON.print(f, Opt_pareto, 4)
	end
end

# ╔═╡ c1802599-c878-47e2-a778-bc3aa4aee606
begin
	# scatter(Opt_dataframe.Purity_Captured_CO2), 					  Opt_dataframe.Intrinsic_capture_efficiency)

	scatter(Opt_pareto[:,"Purity_Captured_CO2"], 
			Opt_pareto[:,"Intrinsic_capture_efficiency"],
			xerror = Opt_pareto[:,"Purity_Captured_CO2_std"],
			yerror = Opt_pareto[:,"Intrinsic_capture_efficiency_std"],
			label = "Pareto Front")

	xlabel!("Purity Captured CO2")
	ylabel!("Intrinsic Capture Efficiency (mol/J)")
end

# ╔═╡ 38dfc2a5-f068-4e5a-8827-21ea82aa44b6
begin
	plot1 = histogram2d(Opt_dataframe[:,"Purity_Captured_CO2"], 
				Opt_dataframe[:,"Intrinsic_capture_efficiency"],
				bins = (100,100),
				# show_empty_bins =true,
				normalize=:none,
				# color = cgrad(:plasma,scale=log))
				color = :plasma,
				# colorbar_title = "Counts"
				)
	plot1 = xlabel!("Purity of Captured CO₂")
	plot1 = ylabel!("Intrinsic Capture Efficency (mol/J)")
	plot1 = title!("a)", titlelocation = :left)

	plot2 = scatter(Opt_pareto[:,"Purity_Captured_CO2"], 
			Opt_pareto[:,"Intrinsic_capture_efficiency"],
			xerror = Opt_pareto[:,"Purity_Captured_CO2_std"],
			yerror = Opt_pareto[:,"Intrinsic_capture_efficiency_std"],
			label = "Pareto Front")

	plot2 = xlabel!("Purity Captured CO₂")
	plot2 = ylabel!("Intrinsic Capture Efficiency (mol/J)")
	plot2 = title!("b)", titlelocation = :left)

	plot(plot1, plot2, layout = (1, 2))
	plot!(size=(900,400), left_margin = 15*Plots.mm, bottom_margin = 10*Plots.mm)
	savefig("/users/asm6/Julia_scripts/OptIDC_HistAndPareto.svg")
end

# ╔═╡ 298af6ba-6364-4761-9c51-48bd037d6306
begin
	# scatter(Opt_dataframe.Purity_Captured_CO2), 					  Opt_dataframe.Intrinsic_capture_efficiency)

	scatter(Opt_pareto[:,"Purity_Captured_CO2"], 
			Opt_pareto[:,"Intrinsic_capture_efficiency"],
			xerror = Opt_pareto[:,"Purity_Captured_CO2_std"],
			yerror = Opt_pareto[:,"Intrinsic_capture_efficiency_std"],
			zcolor = Opt_pareto[:,"Temp_start"],
			color =:coolwarm,
			colorbar_title = "Start Temperature (K)")
	title!("Pareto Start Temperature")
	xlabel!("Purity Captured CO2")
	ylabel!("Intrinsic Capture Efficiency (mol/J)")
end

# ╔═╡ 1cb8e1ec-f8be-46ce-a6cc-1287191c1f21
begin
	# scatter(Opt_dataframe.Purity_Captured_CO2), 					  Opt_dataframe.Intrinsic_capture_efficiency)

	scatter(Opt_pareto[:,"Purity_Captured_CO2"], 
			Opt_pareto[:,"Intrinsic_capture_efficiency"],
			xerror = Opt_pareto[:,"Purity_Captured_CO2_std"],
			yerror = Opt_pareto[:,"Intrinsic_capture_efficiency_std"],
			zcolor = Opt_pareto[:,"Temp_end"],
			color =:coolwarm,
			colorbar_title = "End Temperature (K)")
	title!("Pareto End Temperature")
	xlabel!("Purity Captured CO2")
	ylabel!("Intrinsic Capture Efficiency (mol/J)")
end

# ╔═╡ 032188e5-fea2-4cda-905e-6b03750798c2
begin
	# scatter(Opt_dataframe.Purity_Captured_CO2), 					  Opt_dataframe.Intrinsic_capture_efficiency)

	scatter(Opt_pareto[:,"Purity_Captured_CO2"], 
			Opt_pareto[:,"Intrinsic_capture_efficiency"],
			xerror = Opt_pareto[:,"Purity_Captured_CO2_std"],
			yerror = Opt_pareto[:,"Intrinsic_capture_efficiency_std"],
			zcolor = Opt_pareto[:,"Press_start"],
			color =:coolwarm,
			colorbar_title = "Start Pressure (Pa)")
	title!("Pareto Start Pressure")
	xlabel!("Purity Captured CO2")
	ylabel!("Intrinsic Capture Efficiency (mol/J)")
end

# ╔═╡ 924baccd-0360-41b0-91f1-c7ab344332df
begin
	# scatter(Opt_dataframe.Purity_Captured_CO2), 					  Opt_dataframe.Intrinsic_capture_efficiency)

	scatter(Opt_pareto[:,"Purity_Captured_CO2"], 
			Opt_pareto[:,"Intrinsic_capture_efficiency"],
			xerror = Opt_pareto[:,"Purity_Captured_CO2_std"],
			yerror = Opt_pareto[:,"Intrinsic_capture_efficiency_std"],
			zcolor = Opt_pareto[:,"Press_end"],
			color =:coolwarm,
			colorbar_title = "End Pressure (Pa)")
	title!("Pareto End Pressure")
	xlabel!("Purity Captured CO2")
	ylabel!("Intrinsic Capture Efficiency (mol/J)")
end

# ╔═╡ 07048db3-73f6-46d7-9f27-e92aea6b3524


# ╔═╡ 851bf844-0666-4771-a1e8-ffdbdded8cae


# ╔═╡ b2ebbe4f-ac85-4982-9e76-93d4c8285ec6


# ╔═╡ b2a482d9-b5ae-45e8-b026-907bca21ce49


# ╔═╡ f5165b12-fbf8-462b-a21f-2a1c50f9d938


# ╔═╡ 76bb12ba-e695-4924-b919-f638f7988865


# ╔═╡ 35759f0a-1b96-44e0-a4cb-3b1ace3c6b84


# ╔═╡ 1d633f04-3d86-4ac1-a2f2-c1410c86021f
begin
	Opt_pareto2 = DataFrame(Name =[], 
					  Intrinsic_capture_efficiency = [],
					  Intrinsic_capture_efficiency_std = [],
					  Purity_Captured_CO2 = [],
					  Purity_Captured_CO2_std = [],
					  Captured_CO2 = [],
					  Captured_CO2_std = [],
					  Captured_N2 = [],
					  Captured_N2_std = [],
					  Temp_start = [],
					  Temp_end = [],
					  Press_start = [],
					  Press_end = [])

	sort!(Opt_dataframe, [order("Purity_Captured_CO2", rev = true), order("Intrinsic_capture_efficiency", rev = true)])

	push!(Opt_pareto2, Opt_dataframe[1,:])

	for entry in eachrow(Opt_dataframe[2:end, :])
		if entry["Intrinsic_capture_efficiency"] >= Opt_pareto2[end, "Intrinsic_capture_efficiency"]
			push!(Opt_pareto2, entry)
		end
	end
end

# ╔═╡ fc8bd191-5dfb-4237-9105-a3df83bc9557
begin
	# scatter(Opt_dataframe.Purity_Captured_CO2), 					  Opt_dataframe.Intrinsic_capture_efficiency)

	scatter(Opt_pareto2[:,"Purity_Captured_CO2"], 
			Opt_pareto2[:,"Intrinsic_capture_efficiency"],
			xerror = Opt_pareto2[:,"Purity_Captured_CO2_std"],
			yerror = Opt_pareto2[:,"Intrinsic_capture_efficiency_std"])
end

# ╔═╡ 3e59e5b1-df7c-48d9-982e-510c01ab2f33
begin
	Fuzzy_pareto = DataFrame(Name =[], 
					  Intrinsic_capture_efficiency = [],
					  Intrinsic_capture_efficiency_std = [],
					  Purity_Captured_CO2 = [],
					  Purity_Captured_CO2_std = [],
					  Captured_CO2 = [],
					  Captured_CO2_std = [],
					  Captured_N2 = [],
					  Captured_N2_std = [],
					  Temp_start = [],
					  Temp_end = [],
					  Press_start = [],
					  Press_end = [])

	sort!(Opt_dataframe, [order("Intrinsic_capture_efficiency", rev = true), order("Purity_Captured_CO2", rev = true)])

	push!(Fuzzy_pareto, Opt_dataframe[1,:])

	for entry in eachrow(Opt_dataframe[2:end, :])
		if entry["Purity_Captured_CO2"] >= 0.99 * Fuzzy_pareto[end, "Purity_Captured_CO2"]
			push!(Fuzzy_pareto, entry)
		end
	end
end

# ╔═╡ 6e86e128-15dc-4201-95a1-0a6336794984
begin
	# scatter(Opt_dataframe.Purity_Captured_CO2), 					  Opt_dataframe.Intrinsic_capture_efficiency)

	scatter(Fuzzy_pareto[:,"Purity_Captured_CO2"], 
			Fuzzy_pareto[:,"Intrinsic_capture_efficiency"],
			xerror = Fuzzy_pareto[:,"Purity_Captured_CO2_std"],
			yerror = Fuzzy_pareto[:,"Intrinsic_capture_efficiency_std"])
end

# ╔═╡ 8d47e845-916c-44df-b396-bd30abf611bb
begin
	Fuzzy_pareto2 = DataFrame(Name =[], 
					  Intrinsic_capture_efficiency = [],
					  Intrinsic_capture_efficiency_std = [],
					  Purity_Captured_CO2 = [],
					  Purity_Captured_CO2_std = [],
					  Captured_CO2 = [],
					  Captured_CO2_std = [],
					  Captured_N2 = [],
					  Captured_N2_std = [],
					  Temp_start = [],
					  Temp_end = [],
					  Press_start = [],
					  Press_end = [])

	sort!(Opt_dataframe, [order("Purity_Captured_CO2", rev = true), order("Intrinsic_capture_efficiency", rev = true)])

	push!(Fuzzy_pareto2, Opt_dataframe[1,:])

	for entry in eachrow(Opt_dataframe[2:end, :])
		if log(entry["Intrinsic_capture_efficiency"]) >= 1.0040 * log(Fuzzy_pareto2[end, "Intrinsic_capture_efficiency"])
			push!(Fuzzy_pareto2, entry)
		end
	end
end

# ╔═╡ 56978741-3edf-48c2-819c-e1cb9ad8e88d
begin
	# scatter(Opt_dataframe.Purity_Captured_CO2), 					  Opt_dataframe.Intrinsic_capture_efficiency)

	scatter(Fuzzy_pareto2[:,"Purity_Captured_CO2"], 
			Fuzzy_pareto2[:,"Intrinsic_capture_efficiency"],
			xerror = Fuzzy_pareto2[:,"Purity_Captured_CO2_std"],
			yerror = Fuzzy_pareto2[:,"Intrinsic_capture_efficiency_std"])
end

# ╔═╡ ea1f2731-3bc0-4722-93e1-357176ffd33c
begin
	@show length(Fuzzy_pareto.Name)
	@show length(Fuzzy_pareto2.Name)
end

# ╔═╡ c362758d-6a6e-456f-b011-7e21179fa6dd
begin
	for entry in eachrow(Fuzzy_pareto2)
		push!(Fuzzy_pareto, entry)
	end
	unique(Fuzzy_pareto)
end

# ╔═╡ 694ab0ad-be85-4ddb-805b-6c7cef62290a
begin
	# scatter(Opt_dataframe.Purity_Captured_CO2), 					  Opt_dataframe.Intrinsic_capture_efficiency)

	scatter(Fuzzy_pareto[:,"Purity_Captured_CO2"], 
			Fuzzy_pareto[:,"Intrinsic_capture_efficiency"],
			xerror = Fuzzy_pareto[:,"Purity_Captured_CO2_std"],
			yerror = Fuzzy_pareto[:,"Intrinsic_capture_efficiency_std"])
end

# ╔═╡ 7f81c082-ab14-4756-88c5-459bbcf86723


# ╔═╡ c5400971-ffcb-4127-8eb6-9775d80f8514


# ╔═╡ faee1450-1c78-44de-9445-60c0495448f2


# ╔═╡ 43f3ea60-8fea-40e3-8b7f-a0525f5e0a32


# ╔═╡ cfb00094-1586-428e-8303-d69cfcf21d38


# ╔═╡ fede1fa1-b6f9-419c-a360-801b63ca2c01


# ╔═╡ 0fa26d32-6f31-48ff-a409-f62be1f66ccb


# ╔═╡ 39e62b38-2e27-46e1-98db-14027f1ed413


# ╔═╡ 5f64c655-f6a0-454e-b4a8-84e5ebb5a5ee


# ╔═╡ dfff353a-f94d-418d-8b41-32d7394b75b8


# ╔═╡ 8e87c955-3b4d-4ede-98ff-d4618a701942


# ╔═╡ 26e4b492-0bfe-4ff5-af57-f73de5feb7a4
Opt_dataframe[typeof.(Opt_dataframe.Purity_Captured_CO2) .!= Float64, :]

# ╔═╡ 4917a664-22c0-4ace-b35f-3b99068df515
begin
	wierd_dict = JSON.parsefile(Base_directory*"Optimized_Intrinsic_Cycle/OptIDC_"*"AWAGOV_clean"*".json")

	print("Looks like something went wrong in: objectives_dist = Intrinisic_refresh_objectives_posterior_dist")
end

# ╔═╡ e0b52c3d-2290-4cf3-9d12-deb43b6a3b16


# ╔═╡ 13ef1b83-c736-4ad7-b9c2-bab94e16a286
wierd_dict["Refresh_Path"]["Moles_CO2"][2][1] - wierd_dict["Refresh_Path"]["Moles_CO2"][2][end] 

# ╔═╡ caff2846-852e-424a-b231-f4d3278e7588
wierd_dict["Refresh_Path"]["Moles_N2"][2][1] - wierd_dict["Refresh_Path"]["Moles_N2"][2][end] 

# ╔═╡ b24472ce-be60-4b72-8c6c-3e6ff96df95d
begin
	for thing in wierd_dict["E_Balance"]["Step_2"]["E_to_heat_sorbent"]
		@show minimum(thing)
	end
end

# ╔═╡ Cell order:
# ╠═e0d95bff-efe1-4e6a-be09-cd54644d9df2
# ╠═10d53dca-9147-11ef-3a6d-951ef039c737
# ╠═420a9a15-2ee7-40a3-b14e-33384165bc40
# ╠═ae79d9cd-f16a-4275-a15a-08e710ea8242
# ╠═a0cb4aa3-6180-44f4-b250-ad1ed3603b70
# ╠═54f22445-f4bb-46eb-ae95-3a8aee463a8a
# ╠═a2cb1ac0-289b-4008-ade8-d15720ce3baa
# ╠═39580cbf-7519-4e6f-b9d8-9da07557e175
# ╠═ba133f32-5d9e-4307-9c10-03beaf778c61
# ╠═cefa7523-49db-48bf-8f90-20befd802b9f
# ╠═8a94b9aa-251a-43f8-83f9-ecfe4b6dca2e
# ╠═6d7703aa-b090-4f7e-840a-93160cea2855
# ╠═2b8b8219-dcbf-476f-9b8c-59bdf830ad86
# ╠═f6d2315c-a85d-4e4a-a1a8-5937dc2221ff
# ╠═e2782249-ff89-40c5-8d0f-f33f4be960aa
# ╠═2ba736b8-7b72-4672-af7a-6f68db19e487
# ╠═b4d94a59-f02c-4f05-87a1-96af18c48de5
# ╠═443c7801-e66a-43af-bed5-48a450e8274c
# ╠═62c74a4e-6dad-4b70-a5c2-50075a856f5e
# ╠═37c70cbc-bad7-44d0-abe8-75ba6fb0410c
# ╠═f0a5f4be-d152-4218-8e79-be63f2efbca8
# ╠═973addd2-8c3e-408a-a6a7-90ef2d1be7bc
# ╠═3f3490b6-725d-4ca2-ad1c-518b72dae35b
# ╠═b09ee4e9-6086-4ef1-8d93-04a03be09c6b
# ╠═9ea02bcf-9e8e-4a7c-9d26-ff87dcbf3e06
# ╠═2005253e-4cc2-498d-86a4-2fde787a266c
# ╠═acf56879-d065-4127-be70-fc4e69538345
# ╠═c804ea9d-b399-47f0-8dd5-d9388f6b6c70
# ╠═add66865-896b-4ab3-a45d-0e39186883d0
# ╠═47d70bb2-ff12-49ae-8f85-50bb160b2fa4
# ╠═81e42773-9748-4dd9-91b4-f412abae1655
# ╠═eb926db5-b130-4c3b-888e-a8fc542b1419
# ╠═2d3558c7-635d-4acd-9edc-199275ea97cf
# ╠═1a010f1a-240f-41b9-8c02-0371acac3360
# ╠═7e27220c-12c9-4fdb-a877-526308708de5
# ╠═f8fdb476-5d3d-4706-952a-3a07c8129833
# ╠═9226b1f0-de3b-4fcb-a9e3-ec4423815e56
# ╠═4f3441ac-97ad-451f-8d84-16e20fe3f79a
# ╠═095c2897-0461-4a90-b051-10f2ee360dd6
# ╠═c1802599-c878-47e2-a778-bc3aa4aee606
# ╠═38dfc2a5-f068-4e5a-8827-21ea82aa44b6
# ╠═298af6ba-6364-4761-9c51-48bd037d6306
# ╠═1cb8e1ec-f8be-46ce-a6cc-1287191c1f21
# ╠═032188e5-fea2-4cda-905e-6b03750798c2
# ╠═924baccd-0360-41b0-91f1-c7ab344332df
# ╠═07048db3-73f6-46d7-9f27-e92aea6b3524
# ╠═851bf844-0666-4771-a1e8-ffdbdded8cae
# ╠═b2ebbe4f-ac85-4982-9e76-93d4c8285ec6
# ╠═b2a482d9-b5ae-45e8-b026-907bca21ce49
# ╠═f5165b12-fbf8-462b-a21f-2a1c50f9d938
# ╠═76bb12ba-e695-4924-b919-f638f7988865
# ╠═35759f0a-1b96-44e0-a4cb-3b1ace3c6b84
# ╠═1d633f04-3d86-4ac1-a2f2-c1410c86021f
# ╠═fc8bd191-5dfb-4237-9105-a3df83bc9557
# ╠═3e59e5b1-df7c-48d9-982e-510c01ab2f33
# ╠═6e86e128-15dc-4201-95a1-0a6336794984
# ╠═8d47e845-916c-44df-b396-bd30abf611bb
# ╠═56978741-3edf-48c2-819c-e1cb9ad8e88d
# ╠═ea1f2731-3bc0-4722-93e1-357176ffd33c
# ╠═c362758d-6a6e-456f-b011-7e21179fa6dd
# ╠═694ab0ad-be85-4ddb-805b-6c7cef62290a
# ╠═7f81c082-ab14-4756-88c5-459bbcf86723
# ╠═c5400971-ffcb-4127-8eb6-9775d80f8514
# ╠═faee1450-1c78-44de-9445-60c0495448f2
# ╠═43f3ea60-8fea-40e3-8b7f-a0525f5e0a32
# ╠═cfb00094-1586-428e-8303-d69cfcf21d38
# ╠═fede1fa1-b6f9-419c-a360-801b63ca2c01
# ╠═0fa26d32-6f31-48ff-a409-f62be1f66ccb
# ╠═39e62b38-2e27-46e1-98db-14027f1ed413
# ╠═5f64c655-f6a0-454e-b4a8-84e5ebb5a5ee
# ╠═dfff353a-f94d-418d-8b41-32d7394b75b8
# ╠═8e87c955-3b4d-4ede-98ff-d4618a701942
# ╠═26e4b492-0bfe-4ff5-af57-f73de5feb7a4
# ╠═4917a664-22c0-4ace-b35f-3b99068df515
# ╠═e0b52c3d-2290-4cf3-9d12-deb43b6a3b16
# ╠═13ef1b83-c736-4ad7-b9c2-bab94e16a286
# ╠═caff2846-852e-424a-b231-f4d3278e7588
# ╠═b24472ce-be60-4b72-8c6c-3e6ff96df95d
