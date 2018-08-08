# Bootstrap resampling
# Determine whether mutational profiles of two groups are significantly different (qualitatively)
# Francis Blokzijl

library(MutationalPatterns)
library(ggplot2)
library(reshape2)

# Read mutational profiles
mut_mat_all = read.table("~/surfdrive/ERCC1/Results/SNV/mutational_profile.txt", header = T)

# Load global variables from MutationalPatterns to get triplets for row.names
rownames(mut_mat_all) = TRIPLETS_96

# Subset mut matrix for WT and MUT liver samples
mut_mat_WT_liver = mut_mat_all[,1:3]
mut_mat_MUT_liver = mut_mat_all[,4:6]

# Normalize mutation count matrices for mutation load
mut_mat_WT_liver_norm = prop.table(as.matrix(mut_mat_WT_liver),2)
mut_mat_MUT_liver_norm = prop.table(as.matrix(mut_mat_MUT_liver),2)

# FUNCTION

generate_replicas = function(m, n_rep)
{
  replicas_list = list()
  # for each sample i
  for(i in 1:ncol(m))
  {
    profile = m[,i]
    replicas = matrix(nrow=96,ncol=n_rep)
    for(j in 1:n_rep)
    {
      # total number of mutations in sample
      n = sum(profile)
      # probability vectorof 96 trinucleotide changes
      p = profile/n
      # randomly take sample with replacement
      x = sample(TRIPLETS_96, size = n, prob = p, replace = T)
      # count the triplet appearances
      res = table(factor(x,lev=TRIPLETS_96))
      replicas[,j] = res
    }
    replicas_list[i] = list(replicas)
  }
  # combine replicas in a single matrix
  replicas_combined = do.call(cbind, replicas_list)
  return(replicas_combined)
}

# FUNCTION 
# Calculate centroid of matrix

centroid = function(m)
{
  m_centroid = rowSums(m)/ncol(m)
  return(m_centroid)
}

# FUNCTION

# Randomly select n_samples from the replicas_matrix
# Calculate distance of centroid with centroid of original matrix (m_centroid)
# Perform n_perm permutations
# return vector with distances

profile_distance_permutations = function(replicas_matrix, n_samples, original_matrix, n_perm, dist_measure)
{
  dist_vector = c()
  for(y in 1:n_perm)
  {
    # Randomy select n_samples from the combined replicas matrix
    # Total number of replicates
    n_rep = ncol(replicas_matrix)
    # Sample n_samples without replacement
    select = sample(1:n_rep, n_samples, replace = F)
    # Get profiles of these 3 samples
    replicas_selected = replicas_matrix[,select]
    
    # Calculate distance of replicas centroid with centroid original matrix
    if(dist_measure == "euclidean normalized")
    {
      # Normalize profiles of replicas for mutation load
      replicas_selected_norm = prop.table(replicas_selected, 2)
      # Normalize profiles of original mutation count matrix
      original_matrix_norm = prop.table(as.matrix(original_matrix), 2)
      # Calculate centroids
      c1 = centroid(replicas_selected_norm)
      c2 = centroid(original_matrix_norm)
      # Calculate Euclidean distance
      res = dist(rbind(c1,c2))
    }
    if(dist_measure == "cosine")
    {
      m_centroid = centroid(original_matrix)
      # Calculate cosine similarity with centroid of original samples
      res = 1 - cos_sim(centroid(replicas_selected), m_centroid)
    }
    dist_vector = c(dist_vector, res)
  }
  return(dist_vector)
}

# FUNCTION
# Plot bootstrap distribution
plot_density = function(values, title, ylab, pline)
{
  df = as.data.frame(values)
  N = length(values)
  ggplot(df, aes(x=values)) +
    geom_density(alpha=.3,fill="#FF6666") +
    ggtitle(title) +
    ylab(ylab) +
    theme_bw() +
    geom_vline(aes(xintercept=pline), colour="#990000", linetype="dashed")
}

plot_density = function(values, title, xlab)
{
  # get value for P value = 0.01
  dist_001 = quantile(values, 0.99)
  df = as.data.frame(values)
  ggplot(df, aes(x=values)) +
    geom_density(alpha=.3,fill="#FF6666") +
    ggtitle(title) +
    xlab(xlab) +
    theme_bw() +
    geom_vline(aes(xintercept=dist_001), colour="#990000", linetype="dashed")
}

plot_density2 = function(distances_WT, distances_MUT, title, xlab, line)
{
  df = data.frame(WT = distances_WT, MUT = distances_MUT)
  df = melt(df)
  # get values for P value = 0.01
  dist_001_WT = quantile(distances_WT, 0.99)
  dist_001_MUT = quantile(distances_MUT, 0.99)
  ggplot(df, aes(x=value, group=variable, fill=variable)) +
    geom_density(alpha=.3) +
    ggtitle(title) +
    xlab(xlab) +
    theme_bw() +
    geom_vline(aes(xintercept=dist_001_WT), colour="red", linetype="dashed") +
    geom_vline(aes(xintercept=dist_001_MUT), colour="lightblue", linetype="dashed") +
    geom_vline(aes(xintercept=line), colour="black")
}

# ------ CALCULATE ACTUAL DISTANCES BETWEEN WT & MUT SAMPLES ------

# Euclidean normalized
euc_dist_norm_WT_MUT = dist(rbind(centroid(mut_mat_WT_liver_norm), centroid(mut_mat_MUT_liver_norm))) # 0.05547004

# Cosine distance
cos_dist_WT_MUT = 1 - cos_sim(centroid(mut_mat_WT_liver), centroid(mut_mat_MUT_liver)) # 0.07493675

# ---- GENERATE BOOTSTRAP DISTRIBUTION FOR WT SAMPLES ------

# generate replicas matrix for WT
mut_mat_WT_liver_rep = generate_replicas(mut_mat_WT_liver, 7000)
# calculate distances with
euc_dist_norm_WT_perm = profile_distance_permutations(mut_mat_WT_liver_rep, original_matrix = mut_mat_WT_liver, 
                                                 n_samples = 3, n_perm = 10000, dist_measure = "euclidean normalized")
cos_dist_WT_perm = profile_distance_permutations(mut_mat_WT_liver_rep, original_matrix = mut_mat_WT_liver, 
                                                 n_samples = 3, n_perm = 10000, dist_measure = "cosine")

# ---- GENERATE BOOTSTRAP DISTRIBUTION FOR MUT SAMPLES ------

# generate replicas matrix for MUT
mut_mat_MUT_liver_rep = generate_replicas(mut_mat_MUT_liver, 7000)
# calculate distances with replicas
euc_dist_norm_MUT_perm = profile_distance_permutations(mut_mat_MUT_liver_rep, original_matrix = mut_mat_MUT_liver, 
                                                      n_samples = 3, n_perm = 10000, dist_measure = "euclidean normalized")
cos_dist_MUT_perm = profile_distance_permutations(mut_mat_MUT_liver_rep, original_matrix = mut_mat_MUT_liver, 
                                                 n_samples = 3, n_perm = 10000, dist_measure = "cosine")

# Plot distribution
plot_density2(euc_dist_norm_WT_perm, euc_dist_norm_MUT_perm, "Bootstrap distributions", "Euclidean distance", line = as.vector(euc_dist_norm_WT_MUT))
plot_density2(cos_dist_WT_perm, cos_dist_MUT_perm, "Bootstrap distributions", "Cosine distance", line = as.vector(cos_dist_WT_MUT))

# ----- FIT TO SIGNATURES -----

# FUNCTION
# Calculate distances between original signature contribution vector
# And contribution vector of 3 (combined) replicas

sig_contr_distance_permutations = function(replicas_matrix, n_samples, signatures, original_matrix, n_perm, dist_measure)
{
  # Get signature contributions for centroid of original matrix
  c1 = centroid(original_matrix)
  fit1 = fit_to_signatures(as.matrix(c1), signatures)
  original_contribution = as.vector(prop.table(fit1$contribution,2))
  
  # replicas
  dist_vector = c()
  for(y in 1:n_perm)
  {
    # Randomy select n_samples from the combined replicas matrix
    # Total number of replicates
    n_rep = ncol(replicas_matrix)
    # Sample n_samples without replacement
    select = sample(1:n_rep, n_samples, replace = F)
    # Get profiles of these samples
    replicas_selected = replicas_matrix[,select]
    # Calculate centroid or selected replicas
    c2 = centroid(replicas_selected)
    # Find signature contribution
    fit2 = fit_to_signatures(as.matrix(c2), signatures)
    # Normalize
    replicas_contribution = as.vector(prop.table(fit2$contribution,2))
 
    
    # Calculate distance of replicas centroid with centroid original matrix
    if(dist_measure == "euclidean normalized")
    {
      # Calculate Euclidean distance
      res = dist(rbind(original_contribution, replicas_contribution))
    }
    if(dist_measure == "cosine")
    {
      # Calculate cosine distance
      res = 1 - cos_sim(original_contribution, replicas_contribution)
    }
    dist_vector = c(dist_vector, res)
  }
  return(dist_vector)
} 

# ----- Read COSMIC mutational signatures ------
sp_url <- paste("http://cancer.sanger.ac.uk/cancergenome/assets/", "signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
# Match the order of the mutation types to MutationalPatterns standard
new_order = match(TRIPLETS_96, cancer_signatures$Somatic.Mutation.Type)
# Reorder cancer signatures dataframe
cancer_signatures = cancer_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
# Keep only 96 contributions of the signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])

# Calculate signature contributions of centroids
sig_contr_30_WT_liver = fit_to_signatures(as.matrix(centroid(mut_mat_WT_liver)), cancer_signatures)$contribution
sig_contr_30_MUT_liver = fit_to_signatures(as.matrix(centroid(mut_mat_MUT_liver)), cancer_signatures)$contribution
# Normalize signature contributions
sig_contr_30_WT_liver_norm = prop.table(sig_contr_30_WT_liver,2)
sig_contr_30_MUT_liver_norm = prop.table(sig_contr_30_MUT_liver,2)

# Plot signature contribution of WT liver & MUT liver
plot_contribution(cbind(sig_contr_30_WT_liver, sig_contr_30_MUT_liver))

# Calculate euclidean distance between signature contributions of WT liver & MUT liver
euc_dist_contr_WT_MUT = dist(rbind(t(sig_contr_30_WT_liver_norm), t(sig_contr_30_MUT_liver_norm))) # 0.4516852
# Cosine
cos_dist_contr_WT_MUT = cos_sim(as.vector(sig_contr_30_WT_liver_norm), as.vector(sig_contr_30_MUT_liver_norm)) # 0.6922431

# ------ Generate bootstrap distribution WT -------

# Euclidean

euc_dist_contr_WT_perm = sig_contr_distance_permutations(replicas_matrix = mut_mat_WT_liver_rep,
                                       signatures = cancer_signatures,
                                       n_samples = 3,
                                       original_matrix = mut_mat_WT_liver,
                                       n_perm = 1000,
                                       dist_measure = "euclidean normalized")

euc_dist_contr_MUT_perm = sig_contr_distance_permutations(replicas_matrix = mut_mat_MUT_liver_rep,
                                                         signatures = cancer_signatures,
                                                         n_samples = 3,
                                                         original_matrix = mut_mat_MUT_liver,
                                                         n_perm = 1000,
                                                         dist_measure = "euclidean normalized")

plot_density2(euc_dist_contr_WT_perm, euc_dist_contr_MUT_perm, 
              "Bootstrap distributions", 
              "Euclidean distance \n 30 signature contributions", 
              line = as.vector(euc_dist_contr_WT_MUT))

# Cosine

cos_dist_contr_WT_perm = sig_contr_distance_permutations(replicas_matrix = mut_mat_WT_liver_rep,
                                                         signatures = cancer_signatures,
                                                         n_samples = 3,
                                                         original_matrix = mut_mat_WT_liver,
                                                         n_perm = 1000,
                                                         dist_measure = "cosine")

cos_dist_contr_MUT_perm = sig_contr_distance_permutations(replicas_matrix = mut_mat_MUT_liver_rep,
                                                          signatures = cancer_signatures,
                                                          n_samples = 3,
                                                          original_matrix = mut_mat_MUT_liver,
                                                          n_perm = 1000,
                                                          dist_measure = "cosine")

plot_density2(cos_dist_contr_WT_perm, cos_dist_contr_MUT_perm, 
              "Bootstrap distributions", 
              "Cosine distance \n 30 signature contributions", 
              line = as.vector(cos_dist_contr_WT_MUT))



# ----- RECONSTRUCT WITH SUBSET OF COSMIC SIGNATURES -----

# Signature 1 3 8 9 10 11 14 18 29 30

select = c(1,3,8,9,10,11,14,18,29,30)
COSMIC_subset = cancer_signatures[,select]

# repeat analyses with subset of signatures

# Calculate signature contributions of centroids
sig_contr_10_WT_liver = fit_to_signatures(as.matrix(centroid(mut_mat_WT_liver)), COSMIC_subset)$contribution
sig_contr_10_MUT_liver = fit_to_signatures(as.matrix(centroid(mut_mat_MUT_liver)), COSMIC_subset)$contribution
# Normalize signature contributions
sig_contr_10_WT_liver_norm = prop.table(sig_contr_10_WT_liver,2)
sig_contr_10_MUT_liver_norm = prop.table(sig_contr_10_MUT_liver,2)

# Plot signature contribution of WT liver & MUT liver
plot_contribution(cbind(sig_contr_10_WT_liver, sig_contr_10_MUT_liver))

# Calculate euclidean distance between signature contributions of WT liver & MUT liver
euc_dist_contr_10_WT_MUT = dist(rbind(t(sig_contr_10_WT_liver_norm), t(sig_contr_10_MUT_liver_norm))) # 0.4516852
# Cosine
cos_dist_contr_10_WT_MUT = cos_sim(as.vector(sig_contr_10_WT_liver_norm), as.vector(sig_contr_10_MUT_liver_norm)) # 0.6922431

# ------ Generate bootstrap distribution WT -------

# Euclidean

euc_dist_contr_10_WT_perm = sig_contr_distance_permutations(replicas_matrix = mut_mat_WT_liver_rep,
                                                         signatures = COSMIC_subset,
                                                         n_samples = 3,
                                                         original_matrix = mut_mat_WT_liver,
                                                         n_perm = 1000,
                                                         dist_measure = "euclidean normalized")

euc_dist_contr_10_MUT_perm = sig_contr_distance_permutations(replicas_matrix = mut_mat_MUT_liver_rep,
                                                          signatures = COSMIC_subset,
                                                          n_samples = 3,
                                                          original_matrix = mut_mat_MUT_liver,
                                                          n_perm = 1000,
                                                          dist_measure = "euclidean normalized")

plot_density2(euc_dist_contr_10_WT_perm, euc_dist_contr_10_MUT_perm, 
              "Bootstrap distributions", 
              "Euclidean distance \n 10 signature contributions", 
              line = as.vector(euc_dist_contr_10_WT_MUT))

# Cosine

cos_dist_contr_10_WT_perm = sig_contr_distance_permutations(replicas_matrix = mut_mat_WT_liver_rep,
                                                         signatures = COSMIC_subset,
                                                         n_samples = 3,
                                                         original_matrix = mut_mat_WT_liver,
                                                         n_perm = 1000,
                                                         dist_measure = "cosine")

cos_dist_contr_10_MUT_perm = sig_contr_distance_permutations(replicas_matrix = mut_mat_MUT_liver_rep,
                                                          signatures = COSMIC_subset,
                                                          n_samples = 3,
                                                          original_matrix = mut_mat_MUT_liver,
                                                          n_perm = 1000,
                                                          dist_measure = "cosine")

plot_density2(cos_dist_contr_10_WT_perm, cos_dist_contr_10_MUT_perm, 
              "Bootstrap distributions", 
              "Cosine distance \n 10 signature contributions", 
              line = as.vector(cos_dist_contr_10_WT_MUT))




# Now specifically test Signature 8 contribution

sig8_contr_permutations = function(replicas_matrix, n_samples, signatures, original_matrix, n_perm, dist_measure)
{
  # Get signature contributions for centroid of original matrix
  c1 = centroid(original_matrix)
  fit1 = fit_to_signatures(as.matrix(c1), signatures)
  original_contribution = as.vector(prop.table(fit1$contribution,2))
  
  # replicas
  dist_vector = c()
  for(y in 1:n_perm)
  {
    # Randomy select n_samples from the combined replicas matrix
    # Total number of replicates
    n_rep = ncol(replicas_matrix)
    # Sample n_samples without replacement
    select = sample(1:n_rep, n_samples, replace = F)
    # Get profiles of these samples
    replicas_selected = replicas_matrix[,select]
    # Calculate centroid or selected replicas
    c2 = centroid(replicas_selected)
    # Find signature contribution
    fit2 = fit_to_signatures(as.matrix(c2), signatures)
    # Normalize
    replicas_contribution = as.vector(prop.table(fit2$contribution,2))
    
    # Calculate difference

    res = replicas_contribution[8] - original_contribution[8]

    dist_vector = c(dist_vector, res)
  }
  return(dist_vector)
} 

diff_sig8_contr_WT_MUT = 

diff_sig8_contr_WT_perm = sig8_contr_permutations(replicas_matrix = mut_mat_WT_liver_rep,
                                                         signatures = cancer_signatures,
                                                         n_samples = 3,
                                                         original_matrix = mut_mat_WT_liver,
                                                         n_perm = 1000)

diff_sig8_contr_MUT_perm = sig8_contr_permutations(replicas_matrix = mut_mat_MUT_liver_rep,
                                                      signatures = cancer_signatures,
                                                      n_samples = 3,
                                                      original_matrix = mut_mat_MUT_liver,
                                                      n_perm = 1000)



plot_density2(diff_sig8_contr_WT_perm, diff_sig8_contr_MUT_perm, 
              "Bootstrap distributions", 
              "Difference in Signature 8 contribution",
              line = 0.5)

