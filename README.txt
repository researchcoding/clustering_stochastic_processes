Project code of clustering weakly stationary stochastic processes


Summary:
The Matlab codes contained in this repository are based on the research (see Citation) of clustering weakly stationary stochastic processes. This technique is applicable to the unsupervised learning with respect to time series that are weakly stationary.


Author: Ran Zhao (All rights reserved.)

Content: 
dist_ts.m: Calculation the distance between two time series by their covariance (or correlation) matrices.

main_simulation_study.m: The main function produces the main conclusions as shown in the simulation studies of cited paper (see citation). The major steps include 1) simulating weakly stationary stochastic processes, 2) constructing offline dataset and performing clustering analysis, and 3) constructing online dataset and performing clustering analysis.

misclassify_rate.m: Calculate the misclassification rate of clustering results. The misclassification rate is defined as the number of mis-clustered observations divided by the total number of observations. 

scale_mean.m: Scaled arbitrate matrix’s rows to specified mean value.

sim_wssp_paths.m: The function simulates one example of weakly stationary stochastic process. 

unsup_wssp_offline_algo.m: This unsupervised clustering script uses offline methodology (Algorithm 1 in cited paper) to perform clustering analysis. 

unsup_wssp_online_algo.m: This unsupervised clustering script uses online methodology (Algorithm 2 in cited paper) to perform clustering analysis.


Citation:
Q. Peng, N. Rao and R. Zhao. Unsupervised Clustering of Weakly Stationary Weakly Ergodic Processes using Similarity Measure of Orthogonality. ArXiv, 2017.
