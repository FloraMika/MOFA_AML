library(MOFA2)
library(data.table)

# (Optional) set up reticulate connection with Python
# reticulate::use_python("/Users/ricard/anaconda3/envs/base_new/bin/python", required = T)

###############
## Load data ##
###############

# Multiple formats are allowed for the input data:

## -- Option 1 -- ##
# nested list of matrices, where the first index refers to the view and the second index refers to the group.
# samples are stored in the rows and features are stored in the columns.
# Missing values must be filled with NAs, including samples missing an entire view

# (...)

## -- Option 2 -- ##
# data.frame with columns ["sample","feature","view","group","value"]
# In this case there is no need to have missing values in the data.frame,
# they will be automatically filled in when creating the corresponding matrices

file = "/data/MOFA_v04_input.txt"
data = fread(file)

#######################
# Create MOFA object ##
#######################

MOFAobject <- create_mofa(data)



####################
## Define options ##
####################

# use DNAm, DSRT, Mutations, Proteomics and RNAseq (default alphabetical order)

# Data options
# - scale_views: if views have different ranges/variances, it is good practice to scale each view to unit variance (default is FALSE)
# - scale_groups: if groups have different ranges/variances, it is good practice to scale each group to unit variance (default is FALSE)
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- T

# Model options
# - likelihoods: likelihood per view (options are "gaussian","poisson","bernoulli"). By default, they are guessed internally
# - num_factors: number of factors. By default K=10
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 15
model_opts$likelihoods <-  c("gaussian","gaussian","bernoulli","gaussian","gaussian")


# (Optional) Set stochastic inference options
# Only recommended with very large sample size (>1e6) and when having access to GPUs
# - batch_size: float value indicating the batch size (as a fraction of the total data set: 0.10, 0.25 or 0.50)
# - learning_rate: learning rate (we recommend values from 0.25 to 0.75)
# - forgetting_rate: forgetting rate (we recommend values from 0.1 to 0.5)
# - start_stochastic: first iteration to apply stochastic inference (recommended > 5)
# stochastic_opts <- get_default_stochastic_options(MOFAobject)

##############################################
## Prepare MOFA object and  Train the model ##
##############################################

set.seed(12345)
seeds <- floor(10^7*runif(5))

mofa.list <- list()
for(i in 1:5){
  
  # Training options
  # - maxiter: number of iterations
  # - convergence_mode: "fast", "medium", "slow". For exploration, the fast mode is good enough.
  # - startELBO: initial iteration to compute the ELBO (the objective function used to assess convergence)
  # - freqELBO: frequency of computations of the ELBO (the objective function used to assess convergence)
  # - drop_factor_threshold: minimum variance explained criteria to drop factors while training. Default is -1 (no dropping of factors)
  # - gpu_mode: use GPU mode? This needs cupy installed and a functional GPU, see https://cupy.chainer.org/
  # - verbose: verbose mode?
  # - seed: random seed
  train_opts <- get_default_training_options(MOFAobject)
  train_opts$convergence_mode <- "slow"
  train_opts$drop_factor_threshold <- 0.02
  train_opts$maxiter <- 10000
  train_opts$seed <- seeds[i]
  
  
  mofa.toTrain <- prepare_mofa(MOFAobject,
                               data_options = data_opts,
                               model_options = model_opts,
                               training_options = train_opts
  )
  
  outfile <- paste0("/data/MOFA_v04_train_5_views_",i,".hdf5")
  mofa.list[[i]] <- run_mofa(mofa.toTrain, outfile, use_basilisk = T)
}

save(mofa.list, file="/data/MOFA_v04_train_list.RData")