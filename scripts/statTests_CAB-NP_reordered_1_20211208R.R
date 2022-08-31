# Settings
# Number of permutations
P_num. <- 99999
weighting_methods. <- 'unweighted'
sample_size_multi <- function(x) {
  1
}

# Turn on compiler
library(compiler)
library(mixdist)
enableJIT(3)

# Modified sample function
sample_mod <- function(x) {
  if (length(x) == 1) {
    return(x)
  } else {
    return(sample(x))
  }
}

# Reading in data
file_path <-
  '/Users/stephanie/OneDrive/szmeta/CAB-NP_reordered_input_cole12_n428_max'
base_MDD. <-
  read.csv(
    file.path(file_path, 'tested_CAB-NP_reordered_MDD-HC.csv'),
    check.names = FALSE,
    row.names = 1
  )
base_BD. <-
  read.csv(
    file.path(file_path, 'tested_CAB-NP_reordered_BD-HC.csv'),
    check.names = FALSE,
    row.names = 1
  )
base_SZ. <-
  read.csv(
    file.path(file_path, 'tested_CAB-NP_reordered_SZ-HC.csv'),
    check.names = FALSE,
    row.names = 1
  )

hyper_MDD. <-
  read.csv(
    file.path(file_path, 'effects_CAB-NP_reordered_MDD-HC_hyper.csv'),
    check.names = FALSE,
    row.names = 1
  )
hyper_BD. <-
  read.csv(
    file.path(file_path, 'effects_CAB-NP_reordered_BD-HC_hyper.csv'),
    check.names = FALSE,
    row.names = 1
  )
hyper_SZ. <-
  read.csv(
    file.path(file_path, 'effects_CAB-NP_reordered_SZ-HC_hyper.csv'),
    check.names = FALSE,
    row.names = 1
  )

hypo_MDD. <-
  read.csv(
    file.path(file_path, 'effects_CAB-NP_reordered_MDD-HC_hypo.csv'),
    check.names = FALSE,
    row.names = 1
  )
hypo_BD. <-
  read.csv(
    file.path(file_path, 'effects_CAB-NP_reordered_BD-HC_hypo.csv'),
    check.names = FALSE,
    row.names = 1
  )
hypo_SZ. <-
  read.csv(
    file.path(file_path, 'effects_CAB-NP_reordered_SZ-HC_hypo.csv'),
    check.names = FALSE,
    row.names = 1
  )

n_MDD. <-
  read.csv(
    file.path(file_path, 'n_MDD.csv')
  )
n_BD. <-
  read.csv(
    file.path(file_path, 'n_BD.csv')
  )
n_SZ. <-
  read.csv(
    file.path(file_path, 'n_SZ.csv')
  )

for (ii in 1:nrow(base_MDD.)) {
  base_MDD.[ii, ] <- base_MDD.[ii, ] * sample_size_multi(n_MDD.[ii, 2])
  hyper_MDD.[ii, ] <- hyper_MDD.[ii, ] * sample_size_multi(n_MDD.[ii, 2])
  hypo_MDD.[ii, ] <- hypo_MDD.[ii, ] * sample_size_multi(n_MDD.[ii, 2])
}

for (ii in 1:nrow(base_BD.)) {
  base_BD.[ii, ] <- base_BD.[ii, ] * sample_size_multi(n_BD.[ii, 2])
  hyper_BD.[ii, ] <- hyper_BD.[ii, ] * sample_size_multi(n_BD.[ii, 2])
  hypo_BD.[ii, ] <- hypo_BD.[ii, ] * sample_size_multi(n_BD.[ii, 2])
}

for (ii in 1:nrow(base_SZ.)) {
  base_SZ.[ii, ] <- base_SZ.[ii, ] * sample_size_multi(n_SZ.[ii, 2])
  hyper_SZ.[ii, ] <- hyper_SZ.[ii, ] * sample_size_multi(n_SZ.[ii, 2])
  hypo_SZ.[ii, ] <- hypo_SZ.[ii, ] * sample_size_multi(n_SZ.[ii, 2])
}

# Do the hyper permutations

hyper_MDD.mat. <- matrix(NA, P_num., ncol(base_MDD.))
for (jj in 1:P_num.) {
  samp. <- 1:ncol(base_MDD.)
  base_MDD.perm. <- base_MDD.[, samp.]
  hyper_MDD.perm. <- hyper_MDD.[, samp.]
  for (ii in 1:nrow(hyper_MDD.)) {
    hyper_MDD.perm.[ii, which(base_MDD.perm.[ii, ] > 0)] <-
      sample_mod(hyper_MDD.perm.[ii, which(base_MDD.perm.[ii, ] > 0)])
  }
  hyper_MDD.mat.[jj, ] <-
    colSums(hyper_MDD.perm.) / colSums(base_MDD.perm.)
  print(jj)
}
hyper.p_val. <- c()
for (ii in 1:ncol(hyper_MDD.)) {
  hyper.p_val.[ii] <-
    mean(c(hyper_MDD.mat.[, ii] > (sum(hyper_MDD.[, ii]) / sum(base_MDD.[, ii])), 1))
}
MDD_hyper.p_val. <- hyper.p_val.

hyper_BD.mat. <- matrix(NA, P_num., ncol(base_BD.))
for (jj in 1:P_num.) {
  samp. <- 1:ncol(base_BD.)
  base_BD.perm. <- base_BD.[, samp.]
  hyper_BD.perm. <- hyper_BD.[, samp.]
  for (ii in 1:nrow(hyper_BD.)) {
    hyper_BD.perm.[ii, which(base_BD.perm.[ii, ] > 0)] <-
      sample_mod(hyper_BD.perm.[ii, which(base_BD.perm.[ii, ] > 0)])
  }
  hyper_BD.mat.[jj, ] <-
    colSums(hyper_BD.perm.) / colSums(base_BD.perm.)
  print(jj)
}
hyper.p_val. <- c()
for (ii in 1:ncol(hyper_BD.)) {
  hyper.p_val.[ii] <-
    mean(c(hyper_BD.mat.[, ii] > (sum(hyper_BD.[, ii]) / sum(base_BD.[, ii])), 1))
}
BD_hyper.p_val. <- hyper.p_val.

hyper_SZ.mat. <- matrix(NA, P_num., ncol(base_SZ.))
for (jj in 1:P_num.) {
  samp. <- 1:ncol(base_SZ.)
  base_SZ.perm. <- base_SZ.[, samp.]
  hyper_SZ.perm. <- hyper_SZ.[, samp.]
  for (ii in 1:nrow(hyper_SZ.)) {
    hyper_SZ.perm.[ii, which(base_SZ.perm.[ii, ] > 0)] <-
      sample_mod(hyper_SZ.perm.[ii, which(base_SZ.perm.[ii, ] > 0)])
  }
  hyper_SZ.mat.[jj, ] <-
    colSums(hyper_SZ.perm.) / colSums(base_SZ.perm.)
  print(jj)
}
hyper.p_val. <- c()
for (ii in 1:ncol(hyper_SZ.)) {
  hyper.p_val.[ii] <-
    mean(c(hyper_SZ.mat.[, ii] > (sum(hyper_SZ.[, ii]) / sum(base_SZ.[, ii])), 1))
}
SZ_hyper.p_val. <- hyper.p_val.

# Do the hypo permutations

hypo_MDD.mat. <- matrix(NA, P_num., ncol(base_MDD.))
for (jj in 1:P_num.) {
  samp. <- 1:ncol(base_MDD.)
  base_MDD.perm. <- base_MDD.[, samp.]
  hypo_MDD.perm. <- hypo_MDD.[, samp.]
  for (ii in 1:nrow(hypo_MDD.)) {
    hypo_MDD.perm.[ii, which(base_MDD.perm.[ii, ] > 0)] <-
      sample_mod(hypo_MDD.perm.[ii, which(base_MDD.perm.[ii, ] > 0)])
  }
  hypo_MDD.mat.[jj, ] <-
    colSums(hypo_MDD.perm.) / colSums(base_MDD.perm.)
  print(jj)
}
hypo.p_val. <- c()
for (ii in 1:ncol(hypo_MDD.)) {
  hypo.p_val.[ii] <-
    mean(c(hypo_MDD.mat.[, ii] > (sum(hypo_MDD.[, ii]) / sum(base_MDD.[, ii])), 1))
}
MDD_hypo.p_val. <- hypo.p_val.

hypo_BD.mat. <- matrix(NA, P_num., ncol(base_BD.))
for (jj in 1:P_num.) {
  samp. <- 1:ncol(base_BD.)
  base_BD.perm. <- base_BD.[, samp.]
  hypo_BD.perm. <- hypo_BD.[, samp.]
  for (ii in 1:nrow(hypo_BD.)) {
    hypo_BD.perm.[ii, which(base_BD.perm.[ii, ] > 0)] <-
      sample_mod(hypo_BD.perm.[ii, which(base_BD.perm.[ii, ] > 0)])
  }
  hypo_BD.mat.[jj, ] <-
    colSums(hypo_BD.perm.) / colSums(base_BD.perm.)
  print(jj)
}
hypo.p_val. <- c()
for (ii in 1:ncol(hypo_BD.)) {
  hypo.p_val.[ii] <-
    mean(c(hypo_BD.mat.[, ii] > (sum(hypo_BD.[, ii]) / sum(base_BD.[, ii])), 1))
}
BD_hypo.p_val. <- hypo.p_val.

hypo_SZ.mat. <- matrix(NA, P_num., ncol(base_SZ.))
for (jj in 1:P_num.) {
  samp. <- 1:ncol(base_SZ.)
  base_SZ.perm. <- base_SZ.[, samp.]
  hypo_SZ.perm. <- hypo_SZ.[, samp.]
  for (ii in 1:nrow(hypo_SZ.)) {
    hypo_SZ.perm.[ii, which(base_SZ.perm.[ii, ] > 0)] <-
      sample_mod(hypo_SZ.perm.[ii, which(base_SZ.perm.[ii, ] > 0)])
  }
  hypo_SZ.mat.[jj, ] <-
    colSums(hypo_SZ.perm.) / colSums(base_SZ.perm.)
  print(jj)
}
hypo.p_val. <- c()
for (ii in 1:ncol(hypo_SZ.)) {
  hypo.p_val.[ii] <-
    mean(c(hypo_SZ.mat.[, ii] > (sum(hypo_SZ.[, ii]) / sum(base_SZ.[, ii])), 1))
}
SZ_hypo.p_val. <- hypo.p_val.

# Calculate P-values hyper between diseases

# MDD/BD
store. <- matrix(NA, P_num., ncol(base_MDD.))
base_both. <- rbind(base_MDD.,base_BD.)
effect_both. <- rbind(hyper_MDD.,hyper_BD.)
row1. <- nrow(base_MDD.)
row2. <- nrow(base_BD.)
for (jj in 1:P_num.) {
  samp. <- 1:(row1.+row2.)
  base_both.perm. <- base_both.[samp.,]
  effect_both.perm. <- effect_both.[samp.,]
  for (ii in 1:ncol(store.)) {
    effect_both.perm.[which(base_both.perm.[,ii] > 0),ii] <-
      sample_mod(effect_both.perm.[which(base_both.perm.[,ii] > 0),ii])
  }
  store.[jj,] <- abs(colSums(effect_both.perm.[1:row1.,])/colSums(base_both.perm.[1:row1.,]) - 
    colSums(effect_both.perm.[(1+row1.):(row1.+row2.),])/colSums(base_both.perm.[(1+row1.):(row1.+row2.),]))
  print(jj)
}
test_stats. <- abs(colSums(effect_both.[1:row1.,])/colSums(base_both.[1:row1.,]) - 
                     colSums(effect_both.[(1+row1.):(row1.+row2.),])/colSums(base_both.[(1+row1.):(row1.+row2.),]))
p_val. <- c()
for (ii in 1:ncol(base_MDD.)) {
  p_val.[ii] <- mean(c(store.[,ii],Inf)>=test_stats.[ii])
}
hyper_diffGroup_MDD_BD.p_val. <- p_val.

# BD/SZ
store. <- matrix(NA, P_num., ncol(base_BD.))
base_both. <- rbind(base_BD.,base_SZ.)
effect_both. <- rbind(hyper_BD.,hyper_SZ.)
row1. <- nrow(base_BD.)
row2. <- nrow(base_SZ.)
for (jj in 1:P_num.) {
  samp. <- 1:(row1.+row2.)
  base_both.perm. <- base_both.[samp.,]
  effect_both.perm. <- effect_both.[samp.,]
  for (ii in 1:ncol(store.)) {
    effect_both.perm.[which(base_both.perm.[,ii] > 0),ii] <-
      sample_mod(effect_both.perm.[which(base_both.perm.[,ii] > 0),ii])
  }
  store.[jj,] <- abs(colSums(effect_both.perm.[1:row1.,])/colSums(base_both.perm.[1:row1.,]) - 
                       colSums(effect_both.perm.[(1+row1.):(row1.+row2.),])/colSums(base_both.perm.[(1+row1.):(row1.+row2.),]))
  print(jj)
}
test_stats. <- abs(colSums(effect_both.[1:row1.,])/colSums(base_both.[1:row1.,]) - 
                     colSums(effect_both.[(1+row1.):(row1.+row2.),])/colSums(base_both.[(1+row1.):(row1.+row2.),]))
p_val. <- c()
for (ii in 1:ncol(base_BD.)) {
  p_val.[ii] <- mean(c(store.[,ii],Inf)>=test_stats.[ii])
}
hyper_diffGroup_BD_SZ.p_val. <- p_val.

# MDD/SZ
store. <- matrix(NA, P_num., ncol(base_MDD.))
base_both. <- rbind(base_MDD.,base_SZ.)
effect_both. <- rbind(hyper_MDD.,hyper_SZ.)
row1. <- nrow(base_MDD.)
row2. <- nrow(base_SZ.)
for (jj in 1:P_num.) {
  samp. <- 1:(row1.+row2.)
  base_both.perm. <- base_both.[samp.,]
  effect_both.perm. <- effect_both.[samp.,]
  for (ii in 1:ncol(store.)) {
    effect_both.perm.[which(base_both.perm.[,ii] > 0),ii] <-
      sample_mod(effect_both.perm.[which(base_both.perm.[,ii] > 0),ii])
  }
  store.[jj,] <- abs(colSums(effect_both.perm.[1:row1.,])/colSums(base_both.perm.[1:row1.,]) - 
                       colSums(effect_both.perm.[(1+row1.):(row1.+row2.),])/colSums(base_both.perm.[(1+row1.):(row1.+row2.),]))
  print(jj)
}
test_stats. <- abs(colSums(effect_both.[1:row1.,])/colSums(base_both.[1:row1.,]) - 
                     colSums(effect_both.[(1+row1.):(row1.+row2.),])/colSums(base_both.[(1+row1.):(row1.+row2.),]))
p_val. <- c()
for (ii in 1:ncol(base_MDD.)) {
  p_val.[ii] <- mean(c(store.[,ii],Inf)>=test_stats.[ii])
}
hyper_diffGroup_MDD_SZ.p_val. <- p_val.


# Calculate P-values hypo between diseases

# MDD/BD
store. <- matrix(NA, P_num., ncol(base_MDD.))
base_both. <- rbind(base_MDD.,base_BD.)
effect_both. <- rbind(hypo_MDD.,hypo_BD.)
row1. <- nrow(base_MDD.)
row2. <- nrow(base_BD.)
for (jj in 1:P_num.) {
  samp. <- 1:(row1.+row2.)
  base_both.perm. <- base_both.[samp.,]
  effect_both.perm. <- effect_both.[samp.,]
  for (ii in 1:ncol(store.)) {
    effect_both.perm.[which(base_both.perm.[,ii] > 0),ii] <-
      sample_mod(effect_both.perm.[which(base_both.perm.[,ii] > 0),ii])
  }
  store.[jj,] <- abs(colSums(effect_both.perm.[1:row1.,])/colSums(base_both.perm.[1:row1.,]) - 
                       colSums(effect_both.perm.[(1+row1.):(row1.+row2.),])/colSums(base_both.perm.[(1+row1.):(row1.+row2.),]))
  print(jj)
}
test_stats. <- abs(colSums(effect_both.[1:row1.,])/colSums(base_both.[1:row1.,]) - 
                     colSums(effect_both.[(1+row1.):(row1.+row2.),])/colSums(base_both.[(1+row1.):(row1.+row2.),]))
p_val. <- c()
for (ii in 1:ncol(base_MDD.)) {
  p_val.[ii] <- mean(c(store.[,ii],Inf)>=test_stats.[ii])
}
hypo_diffGroup_MDD_BD.p_val. <- p_val.

# BD/SZ
store. <- matrix(NA, P_num., ncol(base_BD.))
base_both. <- rbind(base_BD.,base_SZ.)
effect_both. <- rbind(hypo_BD.,hypo_SZ.)
row1. <- nrow(base_BD.)
row2. <- nrow(base_SZ.)
for (jj in 1:P_num.) {
  samp. <- 1:(row1.+row2.)
  base_both.perm. <- base_both.[samp.,]
  effect_both.perm. <- effect_both.[samp.,]
  for (ii in 1:ncol(store.)) {
    effect_both.perm.[which(base_both.perm.[,ii] > 0),ii] <-
      sample_mod(effect_both.perm.[which(base_both.perm.[,ii] > 0),ii])
  }
  store.[jj,] <- abs(colSums(effect_both.perm.[1:row1.,])/colSums(base_both.perm.[1:row1.,]) - 
                       colSums(effect_both.perm.[(1+row1.):(row1.+row2.),])/colSums(base_both.perm.[(1+row1.):(row1.+row2.),]))
  print(jj)
}
test_stats. <- abs(colSums(effect_both.[1:row1.,])/colSums(base_both.[1:row1.,]) - 
                     colSums(effect_both.[(1+row1.):(row1.+row2.),])/colSums(base_both.[(1+row1.):(row1.+row2.),]))
p_val. <- c()
for (ii in 1:ncol(base_BD.)) {
  p_val.[ii] <- mean(c(store.[,ii],Inf)>=test_stats.[ii])
}
hypo_diffGroup_BD_SZ.p_val. <- p_val.

# MDD/SZ
store. <- matrix(NA, P_num., ncol(base_MDD.))
base_both. <- rbind(base_MDD.,base_SZ.)
effect_both. <- rbind(hypo_MDD.,hypo_SZ.)
row1. <- nrow(base_MDD.)
row2. <- nrow(base_SZ.)
for (jj in 1:P_num.) {
  samp. <- 1:(row1.+row2.)
  base_both.perm. <- base_both.[samp.,]
  effect_both.perm. <- effect_both.[samp.,]
  for (ii in 1:ncol(store.)) {
    effect_both.perm.[which(base_both.perm.[,ii] > 0),ii] <-
      sample_mod(effect_both.perm.[which(base_both.perm.[,ii] > 0),ii])
  }
  store.[jj,] <- abs(colSums(effect_both.perm.[1:row1.,])/colSums(base_both.perm.[1:row1.,]) - 
                       colSums(effect_both.perm.[(1+row1.):(row1.+row2.),])/colSums(base_both.perm.[(1+row1.):(row1.+row2.),]))
  print(jj)
}
test_stats. <- abs(colSums(effect_both.[1:row1.,])/colSums(base_both.[1:row1.,]) - 
                     colSums(effect_both.[(1+row1.):(row1.+row2.),])/colSums(base_both.[(1+row1.):(row1.+row2.),]))
p_val. <- c()
for (ii in 1:ncol(base_MDD.)) {
  p_val.[ii] <- mean(c(store.[,ii],Inf)>=test_stats.[ii])
}
hypo_diffGroup_MDD_SZ.p_val. <- p_val.


# Calculations of P-values for differences between hyper and hypo

# MDD
store. <- matrix(NA, P_num., ncol(base_MDD.))
base_both. <- rbind(base_MDD.,base_MDD.)
effect_both. <- rbind(hyper_MDD.,hypo_MDD.)
row1. <- nrow(base_MDD.)
row2. <- nrow(base_MDD.)
for (jj in 1:P_num.) {
  samp. <- 1:(row1.+row2.)
  base_both.perm. <- base_both.[samp.,]
  effect_both.perm. <- effect_both.[samp.,]
  for (ii in 1:ncol(store.)) {
    effect_both.perm.[which(base_both.perm.[,ii] > 0),ii] <-
      sample_mod(effect_both.perm.[which(base_both.perm.[,ii] > 0),ii])
  }
  store.[jj,] <- abs(colSums(effect_both.perm.[1:row1.,])/colSums(base_both.perm.[1:row1.,]) - 
                       colSums(effect_both.perm.[(1+row1.):(row1.+row2.),])/colSums(base_both.perm.[(1+row1.):(row1.+row2.),]))
  print(jj)
}
test_stats. <- abs(colSums(effect_both.[1:row1.,])/colSums(base_both.[1:row1.,]) - 
                     colSums(effect_both.[(1+row1.):(row1.+row2.),])/colSums(base_both.[(1+row1.):(row1.+row2.),]))
p_val. <- c()
for (ii in 1:ncol(base_MDD.)) {
  p_val.[ii] <- mean(c(store.[,ii],Inf)>=test_stats.[ii])
}
MDD_diffHypoHyper.pval. <- p_val.

# BD
store. <- matrix(NA, P_num., ncol(base_BD.))
base_both. <- rbind(base_BD.,base_BD.)
effect_both. <- rbind(hyper_BD.,hypo_BD.)
row1. <- nrow(base_BD.)
row2. <- nrow(base_BD.)
for (jj in 1:P_num.) {
  samp. <- 1:(row1.+row2.)
  base_both.perm. <- base_both.[samp.,]
  effect_both.perm. <- effect_both.[samp.,]
  for (ii in 1:ncol(store.)) {
    effect_both.perm.[which(base_both.perm.[,ii] > 0),ii] <-
      sample_mod(effect_both.perm.[which(base_both.perm.[,ii] > 0),ii])
  }
  store.[jj,] <- abs(colSums(effect_both.perm.[1:row1.,])/colSums(base_both.perm.[1:row1.,]) - 
                       colSums(effect_both.perm.[(1+row1.):(row1.+row2.),])/colSums(base_both.perm.[(1+row1.):(row1.+row2.),]))
  print(jj)
}
test_stats. <- abs(colSums(effect_both.[1:row1.,])/colSums(base_both.[1:row1.,]) - 
                     colSums(effect_both.[(1+row1.):(row1.+row2.),])/colSums(base_both.[(1+row1.):(row1.+row2.),]))
p_val. <- c()
for (ii in 1:ncol(base_BD.)) {
  p_val.[ii] <- mean(c(store.[,ii],Inf)>=test_stats.[ii])
}
BD_diffHypoHyper.pval. <- p_val.

# SZ
store. <- matrix(NA, P_num., ncol(base_SZ.))
base_both. <- rbind(base_SZ.,base_SZ.)
effect_both. <- rbind(hyper_SZ.,hypo_SZ.)
row1. <- nrow(base_SZ.)
row2. <- nrow(base_SZ.)
for (jj in 1:P_num.) {
  samp. <- 1:(row1.+row2.)
  base_both.perm. <- base_both.[samp.,]
  effect_both.perm. <- effect_both.[samp.,]
  for (ii in 1:ncol(store.)) {
    effect_both.perm.[which(base_both.perm.[,ii] > 0),ii] <-
      sample_mod(effect_both.perm.[which(base_both.perm.[,ii] > 0),ii])
  }
  store.[jj,] <- abs(colSums(effect_both.perm.[1:row1.,])/colSums(base_both.perm.[1:row1.,]) - 
                       colSums(effect_both.perm.[(1+row1.):(row1.+row2.),])/colSums(base_both.perm.[(1+row1.):(row1.+row2.),]))
  print(jj)
}
test_stats. <- abs(colSums(effect_both.[1:row1.,])/colSums(base_both.[1:row1.,]) - 
                     colSums(effect_both.[(1+row1.):(row1.+row2.),])/colSums(base_both.[(1+row1.):(row1.+row2.),]))
p_val. <- c()
for (ii in 1:ncol(base_SZ.)) {
  p_val.[ii] <- mean(c(store.[,ii],Inf)>=test_stats.[ii])
}
SZ_diffHypoHyper.pval. <- p_val.


# Write out file

all_out. <- rbind(
  c(MDD_hyper.p_val.),
  c(MDD_hypo.p_val.),
  c(BD_hyper.p_val.),
  c(BD_hypo.p_val.),
  c(SZ_hyper.p_val.),
  c(SZ_hypo.p_val.),
  c(MDD_diffHypoHyper.pval.),
  c(BD_diffHypoHyper.pval.),
  c(SZ_diffHypoHyper.pval.),
  c(hyper_diffGroup_MDD_BD.p_val.),
  c(hyper_diffGroup_BD_SZ.p_val.),
  c(hyper_diffGroup_MDD_SZ.p_val.),
  c(hypo_diffGroup_MDD_BD.p_val.),
  c(hypo_diffGroup_BD_SZ.p_val.),
  c(hypo_diffGroup_MDD_SZ.p_val.)
)
colnames(all_out.) <- c(names(base_BD.))
file_name. <-
  paste('uncorrected_', weighting_methods., '.csv', sep = '')
write.csv(
  all_out.,
  file.path(file_path, file_name.),
  row.names = c(
    'MDD_hyper',
    'MDD_hypo',
    'BD_hyper',
    'BD_hypo',
    'SZ_hyper',
    'SZ_hypo',
    'MDD_diffHypoHyper',
    'BD_diffHypoHyper',
    'SZ_diffHypoHyper',
    'MDD_BD_diffGroup_hyper',
    'BD_SZ_diffGroup_hyper',
    'MDD_SZ_diffGroup_hyper',
    'MDD_BD_diffGroup_hypo',
    'BD_SZ_diffGroup_hypo',
    'MDD_SZ_diffGroup_hypo'
  )
)

# Get all p_values
uncorrected_pvalues <- as.matrix(all_out.)
all_p. <- as.numeric(uncorrected_pvalues)

# Conduct FDR on all p_values
library(mixdist)
z.scores. <- qnorm(1 - all_p.)
grouped.data. <- mixgroup(qnorm(1 - all_p.), breaks = 'FD')
grouped.data.[1, 2] <- grouped.data.[1, 2]
grouped.data.[nrow(grouped.data.), 2] <-
  sum(z.scores. == Inf) + grouped.data.[nrow(grouped.data.), 2]

para. <- mixparam(mu = c(0, 2), sigma = c(1, 1))
mix. <-
  mix(
    grouped.data.,
    para.,
    emsteps = 1000,
    iterlim = 1000,
    steptol = 1e-10
  )

# Calculate tau
small. <- which.min(mix.$parameters$mu)
tau0. <-
  mix.$parameters$pi[small.] * dnorm(z.scores., mix.$parameters$mu[small.], mix.$parameters$sigma[small.]) /
  (
    mix.$parameters$pi[small.] * dnorm(z.scores., mix.$parameters$mu[small.], mix.$parameters$sigma[small.]) +
      mix.$parameters$pi[3 - small.] * dnorm(z.scores., mix.$parameters$mu[3 -
                                                                             small.], mix.$parameters$sigma[3 - small.])
  )
tau0.[is.na(tau0.)] <- 1
tau0.sided. <- tau0. 
tau0.sided.[which(z.scores. < mix.$parameters$mu[small.])] <- (tau0.[which(z.scores. < mix.$parameters$mu[small.])] + 1)

# Calculate the FDR at cutoff values corresponding to posterior probabilities
cc <- pmin(sort(tau0.sided.),1)
FDRc <- c()
for (ii in 1:length(cc)) {
  FDRc[ii] <- sum(tau0.sided. * (tau0.sided. <= cc[ii])) /
    sum(tau0.sided. <= cc[ii])
}
p.adj.fdr. <- FDRc
p.adj.fdr.[order(tau0.sided.)] <- FDRc

# Write out fdr csv
corrected_qvalues <- uncorrected_pvalues
corrected_qvalues[,] <- p.adj.fdr.

signif_bool_qvalues <- (corrected_qvalues <= 0.05)
signif_qvalues <- signif_bool_qvalues * corrected_qvalues
signif_pvalues <- signif_bool_qvalues * uncorrected_pvalues

if (all(signif_pvalues <= signif_qvalues)) {
  file_name. <-
    paste('fdr', weighting_methods., '.csv', sep = '_')
  write.csv(
    corrected_qvalues,
    file.path(file_path, file_name.),
    row.names = c(
      'MDD_hyper',
      'MDD_hypo',
      'BD_hyper',
      'BD_hypo',
      'SZ_hyper',
      'SZ_hypo',
      'MDD_diffHypoHyper',
      'BD_diffHypoHyper',
      'SZ_diffHypoHyper',
      'MDD_BD_diffGroup_hyper',
      'BD_SZ_diffGroup_hyper',
      'MDD_SZ_diffGroup_hyper',
      'MDD_BD_diffGroup_hypo',
      'BD_SZ_diffGroup_hypo',
      'MDD_SZ_diffGroup_hypo'
    )
  )
}
