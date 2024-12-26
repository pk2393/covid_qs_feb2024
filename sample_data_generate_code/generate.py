import pandas as pd
from pgmpy.models import BayesianNetwork
from pgmpy.estimators import MaximumLikelihoodEstimator, HillClimbSearch, BicScore
from pgmpy.inference import VariableElimination

# data
data = pd.read_csv("../learn_data_for_generating_sample_data.csv") # original data: not provided in this repository
data_covariates = data.iloc[:, 0:20]
data_social = data.iloc[:, [1, 2, 19, 20]]

# nodes
all_nodes_covariates = list(data_covariates.columns)
all_nodes_social = list(data_social.columns)

# edges
from itertools import combinations
initial_edges_covariates = list(combinations(all_nodes_covariates, 2)) 
initial_edges_covariates = [
    edge for edge in initial_edges_covariates
    if (edge[0] == "pref" and edge[1] == "pos_feb2024")
    or (edge[1] == "pref" and edge[0] in ["age", "sex"])
    or ("pref" not in edge)
]

initial_edges_social = list(combinations(all_nodes_social, 2)) 

# BN analysis
hc_covariates = HillClimbSearch(data_covariates)
best_model_covariates = hc_covariates.estimate(scoring_method=BicScore(data_covariates), epsilon=1e-4, fixed_edges=initial_edges_covariates)

hc_social = HillClimbSearch(data_social)
best_model_social = hc_social.estimate(scoring_method=BicScore(data_social), epsilon=1e-4, fixed_edges=initial_edges_social)

best_model_bn_covariates = BayesianNetwork(best_model_covariates.edges())
best_model_bn_social = BayesianNetwork(best_model_social.edges())


best_model_bn_covariates.fit(data_covariates, estimator=MaximumLikelihoodEstimator)
best_model_bn_social.fit(data_social, estimator=MaximumLikelihoodEstimator)

# Synthesize sample data
from pgmpy.sampling import BayesianModelSampling
sampler_covariates = BayesianModelSampling(best_model_bn_covariates)
sampler_social = BayesianModelSampling(best_model_bn_social)

seed_idx = 9999
n_generate = 10000
synthetic_data_covariates = sampler_covariates.forward_sample(size=n_generate, seed=seed_idx)
synthetic_data_social = sampler_social.forward_sample(size=n_generate, seed=seed_idx)

# synthetic_data_covariates.to_csv("../synthetic_data_covariates.csv", index=False)
# synthetic_data_social.to_csv("../synthetic_data_census.csv", index=False)
