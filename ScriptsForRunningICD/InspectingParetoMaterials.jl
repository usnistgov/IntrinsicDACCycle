### A Pluto.jl notebook ###
# v0.19.47

using Markdown
using InteractiveUtils

# ╔═╡ 94572302-9606-11ef-1b55-5bfa798d902e
begin 
	using Pkg
	Pkg.activate("/users/asm6/Julia_scripts/IntrinsicDACCycle/")
end

# ╔═╡ 3d123533-32e2-44ce-ab83-124dd7dc90f8
begin
	using JSON
	using DataFrames
end
	

# ╔═╡ fb366d67-8866-4692-8b35-25ba4f119247
begin
	using Plots
end

# ╔═╡ c4254b2d-8c96-4368-bf29-adf019912447
begin
	Pareto_Optimized = JSON.parsefile("/users/asm6/Julia_scripts/Optimized_IDC_pareto.json")
end

# ╔═╡ 9f0deced-8bc7-4da5-8ae9-dd690dfbc0a8
list_of_pareto_names = unique(Pareto_Optimized["Name"])

# ╔═╡ 87fd8c31-6bf4-4f28-b3d3-c49d6da6efbc
begin
	Base_directory = "/users/asm6/DAC_data/"
	Optimized_subdirectory = "Optimized_Intrinsic_Cycle/"  

	name = list_of_pareto_names[2]

	optimized_dict = JSON.parsefile(Base_directory*Optimized_subdirectory*"OptIDC_"*name*".json")


end

# ╔═╡ 420524ec-c8aa-4c5b-aab6-3c1dfe6f661d
alog10(x) = abs(Int(floor(log10(x))))

# ╔═╡ 1034ac65-653c-411e-8217-2f0ebc3b7401
begin
	ICE = reshape(optimized_dict["Intrinsic_capture_efficiency"], (1,:))
	ICE_std = reshape(optimized_dict["Intrinsic_capture_efficiency_std"], (1,:))
	round_ICE_std = round.(ICE_std, sigdigits = 2)
	round_ICE = []
	for (x, sigma_x) in zip(ICE, ICE_std)
		digits = alog10.(sigma_x)
		round_x = round(x, digits = digits)
		push!(round_ICE, round_x)
	end
	round_ICE = reshape(round_ICE, (1,:))
	units = fill("mol/J", size(ICE))
	ICE_symbol = fill("ξ", size(ICE))
	
	Purity = reshape(optimized_dict["Purity_captured_CO2"], (1, :))
	Purity_std = reshape(optimized_dict["Purity_captured_CO2_std"], (1, :))
	round_Purity_std = round.(Purity_std, sigdigits = 2)
	round_Purity = []
	for (x, sigma_x) in zip(Purity, Purity_std)
		digits = alog10.(sigma_x)
		round_x = round(x, digits = digits)
		push!(round_Purity, round_x)
	end
	round_Purity = reshape(round_Purity, (1,:))
	units = fill("mol/J", size(ICE))
	
	ICE_symbol = fill("ξ", size(ICE))
	Purity_symbol = fill("x_end", size(ICE))
	
	
	
	Temperatures = optimized_dict["Refresh_Path"]["Temperatures"]
	Pressures = optimized_dict["Refresh_Path"]["Pressures"]

	Henry_CO2 = optimized_dict["Refresh_Path"]["Henry_CO2"]
	Henry_CO2_err = optimized_dict["Refresh_Path"]["Henry_CO2_err"]

	Henry_N2 = optimized_dict["Refresh_Path"]["Henry_N2"]
	Henry_N2_err = optimized_dict["Refresh_Path"]["Henry_N2_err"]

	Qads_CO2 = optimized_dict["Refresh_Path"]["Heat_of_adsorb_CO2"]
	Qads_CO2_err = optimized_dict["Refresh_Path"]["Heat_of_adsorb_CO2_err"]

	Qads_N2 = optimized_dict["Refresh_Path"]["Heat_of_adsorb_N2"]
	Qads_N2_err = optimized_dict["Refresh_Path"]["Heat_of_adsorb_N2_err"]

	Specific_heat = optimized_dict["Refresh_Path"]["Specific_heat_sorbent"]
	Specific_heat_err = optimized_dict["Refresh_Path"]["Specific_heat_sorbent_err"]

	Moles_CO2 = optimized_dict["Refresh_Path"]["Moles_CO2"]
	Moles_N2 = optimized_dict["Refresh_Path"]["Moles_N2"]

	plots = []

	p1 = plot(Temperatures, Pressures, labels = (round_ICE, round_Purity) )
	xlabel!("Temperature (K)")
	ylabel!("Pressure (Pa)")
	title!(name)
	push!(plots, p1)
	
	# p2 = plot(Temperatures, Henry_CO2)
	# # plot(Temperatures, Henry_CO2 - Henry_CO2_err, fillrange = Henry_CO2 + Henry_CO2_err, fillalpha = 0.35)
	# xlabel!("Temperature (K)")
	# ylabel!("Kh (mmol/(kg Pa))")
	# title!(name)
	# push!(plots, p2)
	p2 = plot(Temperatures, Henry_CO2, color = 1, label = nothing)
	counter = 1
	for (T, Kh, Kh_err) in zip(Temperatures, Henry_CO2, Henry_CO2_err)
		if counter == 1
			label = "CO2"
		else
			label = nothing 
		end
		plot!(T, Kh - 2 .* Kh_err, fillrange = Kh + 2 .* Kh_err, fillalpha = 0.15, color = 1, label = label)
		counter += 1
	end
	xlabel!("Temperature (K)")
	ylabel!("Kh (mmol/(kg Pa))")
	title!(name)
	
	plot!(Temperatures, Henry_N2, color = 2, label = nothing)
	counter = 1
	for (T, Kh, Kh_err) in zip(Temperatures, Henry_N2, Henry_N2_err)
		if counter == 1
			label = "N2"
		else
			label = nothing 
		end
		plot!(T, Kh - 2 .* Kh_err, fillrange = Kh + 2 .* Kh_err, fillalpha = 0.15, color = 2, label = label)
		counter += 1
	end
	push!(plots, p2)

	p3 = plot(Temperatures, Qads_CO2, color = 1, label = "CO2")
	for (T, z, z_err) in zip(Temperatures, Qads_CO2, Qads_CO2_err)
		plot!(T, z[1] - 2 .* z_err[1], fillrange = z[1] + 2 .* z_err[1], fillalpha = 0.15, color = 1, label = "CO2")
	end
	xlabel!("Temperature (K)")
	ylabel!("Qads (J/mol)")
	title!(name)
	
	plot!(Temperatures, Qads_N2, color = 2, label = "N2")
	for (T, z, z_err) in zip(Temperatures, Qads_N2, Qads_N2_err)
		plot!(T, z[1] - 2 .* z_err[1], fillrange = z[1] + 2 .* z_err[1], fillalpha = 0.15, color = 2, label = "N2")
	end
	push!(plots, p3)

	p4 = plot(Temperatures, Specific_heat, color = 1)
	for (T, z, z_err) in zip(Temperatures, Specific_heat, Specific_heat_err)
		plot!(T, z[1] - 2 .* z_err[1], fillrange = z[1] + 2 .* z_err[1], fillalpha = 0.15, color = 1)
	end
	xlabel!("Temperature (K)")
	ylabel!("Specific Heat (J/(kg K))")
	title!(name)
	push!(plots, p4)

	i = 1
	p5 = plot()
	for (T, n_CO2, n_N2) in zip(Temperatures, Moles_CO2, Moles_N2)
		plot!(T, n_CO2, labels = "Path $i CO2")
		plot!(T, n_N2, labels = "Path $i N2")
		i += 1
	end
	xlabel!("Temperature (K)")
	ylabel!("Adsorbed Moles (mol/kg)")
	title!(name)
	push!(plots, p5)
	
	# for (T, P, 
	# 	 Kh_CO2, Kh_CO2_err, 
	# 	 Kh_N2, Kh_N2_err,
	# 	 Q_CO2, Q_CO2_err,
	# 	 Q_N2, Q_N2_err,
	# 	 Cv, Cv_err,
	# 	 n_CO2, n_N2) in zip(Temperatures, Pressures,
	# 	 					 Henry_CO2, Henry_CO2_err,
	# 	 					 Henry_N2, Henry_N2_err,
	# 	 					 Qads_CO2, Qads_CO2_err,
	# 	 					 Qads_N2, Qads_N2_err,
	# 	 					 Specific_heat, Specific_heat_err,
	# 	 					 Moles_CO2, Moles_N2)

	# 	plot(T,P)
	# end
end

# ╔═╡ 0a084f4c-2ed1-4591-8684-c17d117288f0
begin
	plot1 = plot(Temperatures, Pressures, labels = (round_ICE, round_Purity) )
	plot1 = xlabel!("Temperature (K)")
	plot1 = ylabel!("Pressure (Pa)")
	plot1 = title!("a)", titlelocation = :left)

	plot2 = plot(Temperatures, Henry_CO2, color = 1, label = nothing)
	counter1 = 1
	for (T, Kh, Kh_err) in zip(Temperatures, Henry_CO2, Henry_CO2_err)
		if counter1 == 1
			label = "CO2"
		else
			label = nothing 
		end
		plot2 = plot!(T, Kh - 2 .* Kh_err, fillrange = Kh + 2 .* Kh_err, fillalpha = 0.15, color = 1, label = label)
		counter1 += 1
	end
	plot2 = xlabel!("Temperature (K)")
	plot2 = ylabel!("Kh (mmol/(kg Pa))")
	plot2 = title!("b)", titlelocation = :left)

	plot2 = plot!(Temperatures, Henry_N2, color = 2, label = nothing)
	counter1 = 1
	for (T, Kh, Kh_err) in zip(Temperatures, Henry_N2, Henry_N2_err)
		if counter1 == 1
			label = "N2"
		else
			label = nothing 
		end
		plot2 = plot!(T, Kh - 2 .* Kh_err, fillrange = Kh + 2 .* Kh_err, fillalpha = 0.15, color = 2, label = label)
		counter1 += 1
	end

	plot(plot1, plot2, layout = (1, 2))
	plot!(size=(900,400), left_margin = 15*Plots.mm, bottom_margin = 10*Plots.mm)
	savefig("/users/asm6/Julia_scripts/OneSorb_ParetoAndKh.svg")
end


# ╔═╡ 81f209d4-f5cc-4dc5-b87c-98287e9c72fc
optimized_dict["E_Balance"]["Step_2"]

# ╔═╡ 16dd6881-0dd4-4433-aa94-041954d2512b
begin
	plot(Temperatures, Henry_CO2)
	for (T, Kh, Kh_err) in zip(Temperatures, Henry_CO2, Henry_CO2_err)
		plot!(T, Kh - 2 .* Kh_err, fillrange = Kh + 2 .* Kh_err, fillalpha = 0.35)
	end
	xlabel!("Temperature (K)")
	ylabel!("Kh (mmol/(kg Pa))")
	title!(name)
end

# ╔═╡ e1a05db5-0989-4328-8bd7-bad3b73d3bd6
begin
	plot(Temperatures, Specific_heat, color = 1)
	for (T, z, z_err) in zip(Temperatures, Specific_heat, Specific_heat_err)
		plot!(T, z[1] - 2 .* z_err[1], fillrange = z[1] + 2 .* z_err[1], fillalpha = 0.15, color = 1)
	end
	xlabel!("Temperature (K)")
	ylabel!("Specific Heat (J/(kg K))")
	title!(name)
end

# ╔═╡ f4fd753f-77ce-4e75-b428-ffd6a616e922
optimized_dict

# ╔═╡ Cell order:
# ╠═94572302-9606-11ef-1b55-5bfa798d902e
# ╠═3d123533-32e2-44ce-ab83-124dd7dc90f8
# ╠═fb366d67-8866-4692-8b35-25ba4f119247
# ╠═c4254b2d-8c96-4368-bf29-adf019912447
# ╠═9f0deced-8bc7-4da5-8ae9-dd690dfbc0a8
# ╠═87fd8c31-6bf4-4f28-b3d3-c49d6da6efbc
# ╠═420524ec-c8aa-4c5b-aab6-3c1dfe6f661d
# ╠═1034ac65-653c-411e-8217-2f0ebc3b7401
# ╠═0a084f4c-2ed1-4591-8684-c17d117288f0
# ╠═81f209d4-f5cc-4dc5-b87c-98287e9c72fc
# ╠═16dd6881-0dd4-4433-aa94-041954d2512b
# ╠═e1a05db5-0989-4328-8bd7-bad3b73d3bd6
# ╠═f4fd753f-77ce-4e75-b428-ffd6a616e922
