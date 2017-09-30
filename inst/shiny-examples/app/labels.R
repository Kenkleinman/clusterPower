# parameters common to at least two outcomes
alphatext <- '&alpha; (alpha)'
alphatooltip <- 'Type I error rate. Should be between 0 and 1, preferably close to 0 (e.g. 0.05).'
alphavalidmsg <- 'Type I error rate should be between 0 and 1.'

powertext <- 'Power (power)'
powertooltip <- 'Power of the test. Should be between 0 and 1, and preferably close to 1 (e.g. 0.80 or 0.90).'
powervalidmsg <- 'Power should be between 0 and 1.'

mtext <- 'Clusters per arm (m)'
mtooltip <- 'The number of clusters per arm.'

ntext <- 'Cluster size (n)'
ntooltip <- 'The mean sample size per cluster.'

icctext <- 'ICC (icc)'
icctooltip <- 'Intracluster correlation coefficient.'

cvtext <- 'Cluster size CV (cv)'
cvtooltip <- 'Coefficient of variation of the cluster sizes. When this equals 0, all clusters have the same size.'

# specific to 2mean:
varwtext <- 'Within variance (varw)'
varwtooltip <- 'Within cluster variance. Assumed to be the same for all clusters.'

dtext <- 'Difference (d)'
dtooltip <- 'Expected difference in condition means.'

methodtext <- 'Unequal Cluster Size Adjustment'
methodtooltip <- 'Method for calculating the variance inflation and design effect due to unequal cluster sizes. When CV = 0, "method" has no effect.'

# specific to 2prop:
p1text <- 'Proportion 1 (p1)'
p1tooltip <- 'The expected proportion in the treatment group.'

p2text <- 'Proportion 2 (p2)'
p2tooltip <- 'The proportion in the control group.'

p1inctext <- 'p1 > p2'
p1inctooltip <- 'Select to indicate that the treatment group proportion is greater than the control group proportion. This selection only matters if the target values are "p1" or "p2".'

pooledtext <- 'Pooled'
pooledtooltip <- 'Select to indicate if pooled variance is desired.'

# specific to 2rate:
r1text <- 'Rate 1 (r1)'
r1tooltip <- 'The expected rate in the treatment group.'

r2text <- 'Rate 2 (r2)'
r2tooltip <- 'The expected rate in the control group.'

pytext <- 'Person-years per cluster (py)'
pytooltip <- 'Person years per cluster.'

cvbtext <- 'Btwn-cluster CV (cvb)'
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