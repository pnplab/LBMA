

# Get all p_values
folder_path <- '/Users/stephanie/OneDrive/szmeta/CAB-NP_reordered_input_cole12_n428_max'
all_out2. <- read.csv(file.path(folder_path,'99999p_n428_cole12_uncorrected_unweighted.csv'), row.names = 1)
uncorrected_pvalues <- as.matrix(all_out2.)
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
cc <- pmin(sort(tau0.sided.), 1)
FDRc <- c()
for (ii in 1:length(cc)) {
  FDRc[ii] <- sum(pmin(tau0.sided.,1) * (pmin(tau0.sided.,1) <= cc[ii])) /
    sum(pmin(tau0.sided.,1) <= cc[ii])
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
    paste('fdr', '.csv', sep = '_')
  write.csv(
    corrected_qvalues,
    file.path(folder_path, file_name.),
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
