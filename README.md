Re-examining the prevalence of chaos in ecological time series, after Rogers et al. (*Nat. Ecol. Evol.*, 2022), while accounting for uncertainty in false positive and false negative rates of chaos detection methods and the presence of nonlinear stochastic dynamics.

For Bayesian analyses of the prevalence of chaos in the GPDD data, use "rstan_code.Rmd".

The script "stochasticity_tests_code.R" provides the code for running stochasticity tests on the GPDD time series. The user-defined functions "embed_udf.R" and "predict_np_udf.R" are associated with this script. The outputs of the delta-epsilon test for presence of determinism are provided in the folder "result_plots". 

The script "simulation_stoch_detect.R" provides the code for running simulations to test the accuracy of methods used to detect chaos and stochasticity in time series. The user-defined functions "nonLinearMaps.R", "npe_heuristic_udf.R" and "predict_np_udf.R" are associated with this script.
