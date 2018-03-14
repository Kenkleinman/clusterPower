# parameters common to at least two outcomes
alphatext <- HTML('&alpha; (alpha)')
alphatooltip <- 'Type I error rate. Should be between 0 and 1, preferably close to 0 (e.g. 0.05).'
alphavalidmsg <- 'Type I error rate should be between 0 and 1.'

powertext <- 'Power (power)'
powertooltip <- 'Power of the test. Should be between 0 and 1, preferably close to 1 (e.g. 0.80 or 0.90).'
powervalidmsg <- 'Power should be between 0 and 1.'

mtext <- 'Clusters per arm (m)'
mtooltip <- 'The number of clusters per arm.'

ntext <- 'Cluster size (n)'
ntooltip <- 'The mean sample size per cluster.'

p1text <- 'Proportion 1 (p1)'
p1tooltip <- 'The expected proportion in the treatment group.'

p2text <- 'Proportion 2 (p2)'
p2tooltip <- 'The proportion in the control group.'

p1inctext <- 'p1 > p2'
p1inctooltip <- 'Select to indicate that the treatment group proportion is greater than the control group proportion. This option is needed only when the target quantity is either "p1" or "p2". If both "p1" and "p2" are given this option has no effect.'

icctext <- 'ICC (icc)'
icctooltip <- 'Intracluster correlation coefficient.'

cvtext <- 'Cluster size CV (cv)'
cvtooltip <- 'Coefficient of variation of the cluster sizes. When this equals 0, all clusters have the same size.'

# -----------------------------------------------------------------------------
# specific to 2mean:
varttext <- 'Total variance (vart)'
varttooltip <- 'The total variation of the outcome.'

dtext <- 'Difference (d)'
dtooltip <- 'Expected difference in condition means.'

methodtext <- 'Unequal Cluster Size Adjustment'
methodtooltip <- 'Method for calculating the variance inflation and design effect due to unequal cluster sizes. When CV = 0, "method" has no effect.'

# -----------------------------------------------------------------------------
# specific to 2meanD:
dDtext <- 'Difference in difference (d)'
dDtooltip <- 'Expected net difference in outcome statistics (means, proportions, etc.).'

rho_ctext <- 'Within cluster correlation (rho_c)'
rho_ctooltip <- 'The correlation between baseline and post-test outcomes at the cluster level. Used in both cohort and cross-sectional designs. A value of "0" is a conservative estimate.'

rho_stext <- 'Within subject correlation (rho_s)'
rho_stooltip <- 'The correlation between baseline and post-test outcomes at the subject level. For a purely cross-sectional design, this value should be 0.'

#
#
rho_mtext <- 'Matching correlation (rho_m)'
rho_mtooltip <- 'The correlation in outcome used between matched clusters.'



# -----------------------------------------------------------------------------
# specific to 2prop:
pooledtext <- 'Pooled'
pooledtooltip <- 'Select to indicate if pooled variance is desired.'


# -----------------------------------------------------------------------------
# specific to 2propD:
ptext <- "Expected mean proportion (p)"
ptooltip <- "The expected mean proportion at the post-test, averaged across treatment and control arms."

covdftext <- "Covariate degrees of freedom (covdf)"
covdftooltip <- "The degrees of freedom used by group-level covariates. A value of 0 means no regression adjustment using covariates."

pvar_ctext <- "Cluster-level proportion of variance explained (pvar_c)"
pvar_ctooltip <- "The expected cluster-level proportion of variance to be explained by regression adjustment for covariates."

pvar_stext <- "Subeject-level proportion of variance explained (pvar_s)"
pvar_stooltip <- "The expected subject-level proportion of variance to be explained by regression adjustment for covariates."

# -----------------------------------------------------------------------------
# specific to 2propM:
cvmtext <- "Within-pair Outcome CV (cvm)"
cvmtooltip <- "The coefficient of variation in the outcome within matched clusters."

# specific to 2rate:
r1text <- 'Rate 1 (r1)'
r1tooltip <- 'The expected rate in the treatment group.'

r2text <- 'Rate 2 (r2)'
r2tooltip <- 'The expected rate in the control group.'

r1inctext <- 'r1 > r2'
r1inctooltip <- 'Select to indicate that the treatment group rate is greater than the control group rate. This option is needed only when the target quantity is either "r1" or "r2". If both "r1" and "r2" are given this option has no effect.'

pytext <- 'Person-years per cluster (py)'
pytooltip <- 'Person years per cluster.'

cvbtext <- 'Between-cluster CV (cvb)'
cvbtooltip <- 'The coefficient of variation of the person years per cluster. Analogous to ICC for two continuous outcomes.'

# static/button text
defaulttext <- 'Defaults'
clearalltext <- 'Clear All'
calctext <- 'Calculate'
dltext <- 'Download'
powercheck <- 'Power must be between 0 and 1.'
alphacheck <- 'Type I error rate, &alpha;, must be between 0 and 1.'
credittext <- 'App created by Jon Moyer and Ken Kleinman; support from NIGMS grant R01GM121370.'

# graphs
ylab <- 'Y'
xlab <- 'X'
grouplab <- 'Group'
colorlab <- 'Color by Group'
rowlab <- 'Facet Row'
collab <- 'Facet Column'
heightlab <- 'Plot Height'
psizelab <- 'Point Size'
lsizelab <- 'Line Size'

# labels for the various quantities
graph_labels <- c(
  `alpha` = "Level of Significance",
  `power` = "Power",
  `d` = "Difference",
  `m` = "Mean Clusters Per Arm",
  `n` = "Mean Cluster Size",
  `icc` = "ICC",
  `varw` = "Within-Cluster Variance",
  `cv` = "Cluster Size CV",
  `rho_c` = "Within-Cluster Correlation",
  `rho_s` = "Within-Subject Correlation",
  `rho_m` = "Matching Correlation",
  `p1` = "Treatment Proportion",
  `p2` = "Control Propotion",
  `cvm` = "Within-Pair Outcome CV",
  `r1` = "Treatment Rate",
  `r2` = "Control Rate",
  `py` = "Person-Years Per Cluster",
  `cvb` = "Between-Cluster CV") 
