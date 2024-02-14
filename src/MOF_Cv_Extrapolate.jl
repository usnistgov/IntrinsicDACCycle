
"""
This function will read in the JSON file for the MOF that contains the Cv predictions.
    Using all of those predictions (and thier uncertainties) at the specified temperature as training data,
    the function will generate a GP and extrapolate the predictions to the temperatures in Temperatures_test.
"""
function Extrapolate_Cv(directory, name, Temperatures_test)  
    
    # file = directory*"Cv_"*name*"_clean.json"  
    file = directory*"/CV_predictions/"*"Cv_"*name*".json"  
    #Read the Json file
    Cv_dict = JSON.parsefile(file)
    Temperatures = convert(Array{Float64,1},Cv_dict["Temperatures"]) #K
    Cᵥ = convert(Array{Float64,1},Cv_dict["Cv"]) #J kg⁻¹K⁻¹
    σ_Cᵥ = convert(Array{Float64,1},Cv_dict["CV_err"]) #J kg⁻¹K⁻¹

    #Build a GP
    μ_prior = MeanZero() #Prior mean
    K1 = SE(log(10), log(100.0)) #Squared Exponential Kernel
    K2 = Const(log(500)) #Constant Kernel
    K = K1+K2 #Sum the kernels

    logObsNoise = log.(σ_Cᵥ) #Specified oberservation noise

    gp = GP(Temperatures,Cᵥ,μ_prior,K,logObsNoise) #Define GP
    GaussianProcesses.optimize!(gp; noise=false, method = ConjugateGradient()) #Optimize GP

    μ, σ² = predict_f(gp, Temperatures_test, full_cov=true) #Predict GP
    μ = reshape(μ, (:,1))
    σ = reshape(sqrt.(diag(σ²)), (:,1))
    return μ, σ #[J/(kg K)]
end
