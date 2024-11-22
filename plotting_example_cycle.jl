### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# ╔═╡ 4faa61c2-b62c-11ee-1952-6b8afe096c06
using Pkg

# ╔═╡ 4990f702-0416-4be9-86b6-7afbbffc160c
Pkg.activate("/users/asm6/Julia_scripts/IntrinsicDACCycle")

# ╔═╡ b3ed0739-c814-4663-9853-506d9124b188
using Revise

# ╔═╡ a78c995d-5bba-4361-a27e-be10030f1a58
begin
	using Plots
	using Distributed
	using JSON
	using DataFrames

end

# ╔═╡ 9517e962-fcdb-4fe3-a077-5799ebaa61b7
using IntrinsicDACCycle

# ╔═╡ 4d15a7c0-47c6-48b7-93dc-cbd0ad680660
Base_directory = "/users/asm6/DAC_data"

# ╔═╡ 544bd32c-4b01-4d73-85c6-ccb2a19bcb14
begin
	#get all the material files
	list_of_material_files = filter(x -> occursin.(".json",x), readdir(Base_directory*"/Intrinsic_cycle/"))
	#strip off the .json tag
    list_of_materials = replace.(list_of_material_files, ".json" => "")
	#filter for _clean matierals
	# list_of_clean_materials = filter(x -> occursin.("_clean", x), list_of_materials)
end

# ╔═╡ 9496c6a4-ba82-4cff-a342-a6c4affad650
begin
	# name = list_of_materials[1]
	name = "Intrinsic_cyle_ZIJSAO_clean"
	# name = "Intrinsic_cyle_BUSQIQ_clean"
	full_file = Base_directory*"/Intrinsic_cycle/"*name*".json"
	
	
	results_dict = JSON.parsefile(full_file)
end

# ╔═╡ 228faca3-168f-452e-bebe-bd819321733e


# ╔═╡ 8bc5cdf8-ae08-4606-a8e2-e6527a3e520b
begin
	path_Ts = results_dict["Refresh_Path"]["Temperatures"]
	path_Ps = results_dict["Refresh_Path"]["Pressures"]

	step1_T = [path_Ts[1], path_Ts[1]]
	step1_P = [path_Ps[end], path_Ps[1]]

	step3_T = [path_Ts[end], path_Ts[1]]
	step3_P = [path_Ps[end], path_Ps[end]]
	
	path_plot = plot(step1_T, step1_P, marker =(1, 5, 0.8 ), label = "Step 1" )
	path_plot = plot!(path_Ts, path_Ps, marker =(1, 5, 0.8 ), label = "Step 2")
	path_plot = plot!(step3_T, step3_P, marker =(1, 5, 0.8 ), label = "Step 3")
	path_plot = plot!(xlabel = "Temperature (K)")
	path_plot = plot!(ylabel = "Pressure (Pa)")
	path_plot = plot!(legend = :inside)
	path_plot = plot!(xlim = (295, 355))
	path_plot = plot!(ylim = (4.8e4, 1.03e5))
	
end

# ╔═╡ f839cf32-15f2-4a1c-8d70-c988659b5562


# ╔═╡ 4bcbc1ae-42d7-4f1f-b2d9-cf4dac7378ae
begin
	n_CO2 = results_dict["Refresh_Path"]["Moles_CO2"]
	n_N2 = results_dict["Refresh_Path"]["Moles_N2"]
	
	uptake_plot = plot(n_CO2, label = "CO2", marker =(1, 5, 0.8), markerstrokewidth=0.5)
	uptake_plot = plot!(n_N2, label = "N2", marker =(1, 5, 0.8 ), markerstrokewidth=0.5)
	uptake_plot = plot!(yaxis = :log10)
	uptake_plot = plot!(xlabel = "Path Index")
	uptake_plot = plot!(ylabel = "Uptake (mol/kg)")
end

# ╔═╡ 654fcf7f-91b1-4337-b1ca-1406e9274647


# ╔═╡ 0f6cac94-d9b9-4802-8042-b97060c6c0f1


# ╔═╡ 8853d631-eb8b-4fcf-bfca-1c77d586505c


# ╔═╡ 087fe079-6e6f-42ea-92ae-dca11587b311


# ╔═╡ 0bb992e1-69fe-4743-a74a-1ded9ef42159


# ╔═╡ aa07906e-809a-4305-864a-a04921498773
begin
	E_to_heat_sorbent = results_dict["E_Balance"]["Step_2"]["E_to_heat_sorbent"][1]
	
	Heat_to_desorb_N2 = results_dict["E_Balance"]["Step_2"]["Heat_to_desorb_N2"][1]
	Work_to_desorb_N2 = results_dict["E_Balance"]["Step_2"]["Work_to_desorb_N2"]
	E_to_heat_adsorbed_N2 = results_dict["E_Balance"]["Step_2"]["E_to_heat_adsorbed_N2"]

	Heat_to_desorb_CO2 = results_dict["E_Balance"]["Step_2"]["Heat_to_desorb_CO2"][1]
	Work_to_desorb_CO2 = results_dict["E_Balance"]["Step_2"]["Work_to_desorb_CO2"]
	E_to_heat_adsorbed_CO2 = results_dict["E_Balance"]["Step_2"]["E_to_heat_adsorbed_CO2"]


	E_to_change_pressure = results_dict["E_Balance"]["Step_2"]["E_to_change_pressure"]

	E_plot = plot(E_to_heat_sorbent[2:end], label = "E sorb", marker =(100, 5, 0.8 ), markerstrokewidth=0)
	E_plot = plot!(Heat_to_desorb_CO2[2:end], label = "Q CO2", marker =(1, 5, 0.8 ), markerstrokewidth=0)
	E_plot = plot!(Work_to_desorb_CO2[2:end], label = "W CO2", marker =(1, 5, 0.8 ), markerstrokewidth=0)
	E_plot = plot!(E_to_heat_adsorbed_CO2[2:end], label = "E CO2", marker =(1, 5, 0.8 ), markerstrokewidth=0)
	
	E_plot = plot!(Heat_to_desorb_N2[2:end], label = "Q N2", marker =(1, 5, 0.8 ), markerstrokewidth=0)
	E_plot = plot!(Work_to_desorb_N2[2:end], label = "W N2", marker =(1, 5, 0.8 ), markerstrokewidth=0)
	E_plot = plot!(E_to_heat_adsorbed_N2[2:end], label = "E N2", marker =(1, 5, 0.8 ), markerstrokewidth=0)

	E_plot = plot!(E_to_change_pressure[2:end], label = "E ΔP", c = 12, marker =(1, 5, 0.8 ), markerstrokewidth=0)

	E_plot = plot!(xlabel = "Path Index")
	E_plot = plot!(ylabel = "Energy (J/kg)")
end

# ╔═╡ 40eb7a37-e1d5-42b8-b066-9771b05304da


# ╔═╡ b5e9619c-c874-470b-8172-1ba8c179aa55


# ╔═╡ 392c29ae-1eb2-47f2-8035-a27a9ef24143


# ╔═╡ d093215c-65a0-4ef8-acf1-3627777365e0


# ╔═╡ a4d14d84-8836-49b7-bbfe-2d56f249742e


# ╔═╡ 88e76da7-93f8-4766-a4e4-ce3623989801
begin
	k_CO2 = results_dict["Refresh_Path"]["Henry_CO2"]
	k_CO2_err = results_dict["Refresh_Path"]["Henry_CO2_err"]
	k_CO2_upper = k_CO2 .+ 1.95 .* k_CO2_err
	k_CO2_lower = k_CO2 .- 1.95 .* k_CO2_err
	
	k_N2 = results_dict["Refresh_Path"]["Henry_N2"]
	k_N2_err = results_dict["Refresh_Path"]["Henry_N2_err"]
	k_N2_upper = k_N2 .+ 1.95 .* k_N2_err
	k_N2_lower = k_N2 .- 1.95 .* k_N2_err
	
	cv_sorb = results_dict["Refresh_Path"]["Specific_heat_sorbent"]
	cv_sorb_err = results_dict["Refresh_Path"]["Specific_heat_sorbent_err"]
	cv_upper = cv_sorb .+ 1.95 .* cv_sorb_err
	cv_lower = cv_sorb .- 1.95 .* cv_sorb_err

	sorb_plot = plot(k_CO2, label = "KH CO2", 
		# marker =(1, 5, 0.8 )
		)
	sorb_plot =	plot!(k_CO2_lower[1], fillrange=k_CO2_upper[1], c=1, alpha = 0.35, label = "KH CO2 CI")
	sorb_plot =plot!(k_N2, label = "KH N2", 
		# marker =(1, 5, 0.8 ),
		c=2)
	sorb_plot =plot!(k_N2_lower[1], fillrange=k_N2_upper[1], c=2, alpha = 0.35, label = "KH N2 CI")
	sorb_plot = plot!(xlabel = "Path Index")
	sorb_plot = plot!(ylabel = "KH (mmol/(kg Pa))")
	sorb_plot = plot!(legend = :bottomright)
	axis2 = twinx()
	sorb_plot = plot!(axis2, cv_sorb, c = 3, label = "Cv Sorb", ylabel = "Specific Heat (J/(kg K))", marker =(1, 5, 0.8 ))
	sorb_plot = plot!(axis2, cv_lower[1], fillrange=cv_upper[1], c = 3, alpha = 0.35, label = "Cv Sorb CI")

	
	# legend()
	# plot!(twinx(), ylabel = "Specific Heat (J/(kg K))")
	
end

# ╔═╡ 9ae7bb9f-d035-410d-ad8f-74439bb58e97
plot(path_plot, sorb_plot, uptake_plot, E_plot)

# ╔═╡ f13a0ec8-ceb8-411d-94dc-17ad405dbac1
begin
	x_s = []
	ξ_s = []

	for name in list_of_materials
		full_file = Base_directory*"/Intrinsic_cycle/"*name*".json"
		results_dict = JSON.parsefile(full_file)

		x = results_dict["Purity_captured_CO2"]
		ξ = results_dict["Intrinsic_capture_efficiency"]

		append!(x_s, x)
		append!(ξ_s, ξ)
	end
end

# ╔═╡ 66b7b9c3-d1ab-466f-9e85-de49e0554d95
begin
	sorted_indexs = sortperm(x_s, rev = true)
	sorted_x_s = x_s[sorted_indexs]
	sorted_ξ_s = ξ_s[sorted_indexs]

	pareto_x_s = [sorted_x_s[1]]
	pareto_ξ_s = [sorted_ξ_s[1]]

	for (i, (ξ, x)) in enumerate(zip(sorted_ξ_s[2:end], sorted_x_s[2:end]))
		if ξ > pareto_ξ_s[end]
			append!(pareto_x_s, x)
			append!(pareto_ξ_s, ξ)
		end
	end
		
end

# ╔═╡ 9d41adfe-2c14-4ed9-a462-df05c25fc98b


# ╔═╡ aacbf2f7-af48-4950-929e-af8449d87450
pareto_x_s

# ╔═╡ 8a583093-6121-484b-8f50-0a4730fe5613
begin
	scatter(ξ_s, x_s)
	plot!(pareto_ξ_s, pareto_x_s, marker = (1,5,0.8))
	plot!(xlabel = "Capture Efficiency (mol/J)")
	plot!(ylabel = "CO2 Purity")
	plot!(legend = nothing)
end

# ╔═╡ 1e359091-4659-4a64-8741-e07bb8ee5795
list_of_material_files[argmax(ξ_s)]

# ╔═╡ 0422aca9-4164-4692-a1f8-4ca99dbf2290


# ╔═╡ Cell order:
# ╠═4faa61c2-b62c-11ee-1952-6b8afe096c06
# ╠═4990f702-0416-4be9-86b6-7afbbffc160c
# ╠═a78c995d-5bba-4361-a27e-be10030f1a58
# ╠═b3ed0739-c814-4663-9853-506d9124b188
# ╠═9517e962-fcdb-4fe3-a077-5799ebaa61b7
# ╠═4d15a7c0-47c6-48b7-93dc-cbd0ad680660
# ╠═544bd32c-4b01-4d73-85c6-ccb2a19bcb14
# ╠═9496c6a4-ba82-4cff-a342-a6c4affad650
# ╠═228faca3-168f-452e-bebe-bd819321733e
# ╠═8bc5cdf8-ae08-4606-a8e2-e6527a3e520b
# ╠═f839cf32-15f2-4a1c-8d70-c988659b5562
# ╠═4bcbc1ae-42d7-4f1f-b2d9-cf4dac7378ae
# ╠═654fcf7f-91b1-4337-b1ca-1406e9274647
# ╠═0f6cac94-d9b9-4802-8042-b97060c6c0f1
# ╠═8853d631-eb8b-4fcf-bfca-1c77d586505c
# ╠═087fe079-6e6f-42ea-92ae-dca11587b311
# ╠═0bb992e1-69fe-4743-a74a-1ded9ef42159
# ╠═aa07906e-809a-4305-864a-a04921498773
# ╠═40eb7a37-e1d5-42b8-b066-9771b05304da
# ╠═b5e9619c-c874-470b-8172-1ba8c179aa55
# ╠═392c29ae-1eb2-47f2-8035-a27a9ef24143
# ╠═d093215c-65a0-4ef8-acf1-3627777365e0
# ╠═a4d14d84-8836-49b7-bbfe-2d56f249742e
# ╠═88e76da7-93f8-4766-a4e4-ce3623989801
# ╠═9ae7bb9f-d035-410d-ad8f-74439bb58e97
# ╠═f13a0ec8-ceb8-411d-94dc-17ad405dbac1
# ╠═66b7b9c3-d1ab-466f-9e85-de49e0554d95
# ╠═9d41adfe-2c14-4ed9-a462-df05c25fc98b
# ╠═aacbf2f7-af48-4950-929e-af8449d87450
# ╠═8a583093-6121-484b-8f50-0a4730fe5613
# ╠═1e359091-4659-4a64-8741-e07bb8ee5795
# ╠═0422aca9-4164-4692-a1f8-4ca99dbf2290
