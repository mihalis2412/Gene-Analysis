# if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install(version = "3.14")
# BiocManager::install("leukemiasEset")
# BiocManager::install("qvalue")

## Exercise 1 ##
## Question 1 & 3 ##
library("leukemiasEset")
data(leukemiasEset)
x <- exprs(leukemiasEset)
str(x)
summary(x)
length(rownames(x)) 
# Number of genes profiled

types <- leukemiasEset$LeukemiaType
summary(types) 
# How many patients in each type of leukemia group

new_mat <- x[,leukemiasEset$LeukemiaType=='AML' | leukemiasEset$LeukemiaType=='NoL']
dim(new_mat) 
# Check that we have the right data set to work with

summary(new_mat) 
# Summary statistics

which(is.na(new_mat)) 
# No missing values

# Plot the histogram of the whole expression matrix and check the data distribution
hist(x, breaks = 100)
# We observe it's right skewed!


# Deeper analysis of the whole expression matrix using the histogram
quantile_Data <- quantile(x)
med <- quantile_Data[3]
q25 <- quantile_Data[2]
q75 <- quantile_Data[4]
me <- mean(x)

hist(x, breaks = 100)
abline(v = med, col = "red")
abline(v = q25, col = "green")
abline(v = q75, col = "blue")
abline(v = me, col = "yellow")
legend("topright", col = c("red", "green", "blue", "yellow"), legend = c("median", 
"25 quantile", "75 quantile", "mean"), pch = 20)



# Another way to check the data distribution of our data set 
qqnorm(new_mat)



# Creating a small multiple chart
df1 <- melt(as.data.frame(new_mat))
ggplot(data = df1, aes(x = value)) + 
  stat_density() + 
  facet_wrap(~variable, scales = "free")
# We observe that our data are nowhere near being normal
# Right skewed data so the log transformation should be more appropriate



# Log transformation
df2 <- melt(as.data.frame(log(new_mat)))
ggplot(data = df2, aes(x = value)) + 
  stat_density() + 
  facet_wrap(~variable, scales = "free")
# We observe that the picture now looks closer to normality than before!



# Useful summary statistics
apply(new_mat,2,range)
#  It’s not often used because it’s very sensitive to outliers

apply(new_mat,2,IQR)
# It’s pretty robust to outliers. It’s used a lot in combination with the median.

apply(new_mat,2,median)

apply(new_mat,2,mean)




# Checking the correlations
library(Hmisc)
mat <- as.matrix(new_mat)
newmat <- mat[,1:5]
rcorr(newmat)
# Really strong correlations observed


# Examine data distributions for all individual patients using colored boxplots
boxplot(x, horizontal = T, las = 2, col = leukemiasEset$LeukemiaType)
legend("topright", col = unique(as.numeric(leukemiasEset$LeukemiaType)), legend = unique(leukemiasEset$LeukemiaType), 
pch = 20)


# Boxplots of the 2 groups of people (those with AML and the control group)
boxplot(new_mat, las = 2, col = rep(c(1,2), each = 12), main = "Boxplots of the 2 groups")
legend('topleft',legend = c('AML','NOL'),fill = c('black','red'))



# Multiple means testing for the 2 groups (people with AML and healthy ones)
n <- dim(new_mat)[1]
pval <- numeric(n)
diff <- numeric(n)
std.err <- numeric(n)
for (i in 1:n) {
  pval[i] <- t.test(new_mat[i,1:12],new_mat[i,13:24],alternative = c('two.sided'), var.equal = T)$p.value
  diff[i] <- t.test(new_mat[i,1:12],new_mat[i,13:24],alternative = c('two.sided'), var.equal = T)$estimate
  std.err[i] <-  t.test(new_mat[i,1:12],new_mat[i,13:24],alternative = c('two.sided'), var.equal = T)$stderr
}
hist(pval)

estimated_var <- sqrt(std.err)

# Create a data frame that has the estimated means and their corresponding p_values
dat <- data.frame(diff,pval)
important_dat <- dat[dat$pval < 0.05,]



# Plot the differences of the means of the 2 groups with their corresponding p_values
ggplot(dat, aes(diff, pval)) +                       
  geom_point(aes(color = pval), size = 1) +
  scale_color_gradient(low = "yellow", high = "darkblue") 
# We observe that the darker the color the less the p_value hence when the p_values are for example less than 0.05 we reject H0 => there's a difference in the means of the 2 groups 




# Estimating means for each group
n <- dim(new_mat)[1]
m1 <- numeric(n)
m2 <- numeric(n)
for (i in 1:n) {
  m1[i] <- mean(new_mat[i,1:12])
  m2[i] <- mean(new_mat[i,13:24])  
}

# Plotting the means for each group
plot(m1,m2, col=c('red','green'), main='Scatterplot of the estimated means of the 2 groups', xlab='Means of patients with AML',ylab='Means of healthy people')
legend('topleft',legend = c('AML','NOL'),fill = c('red','green'))



# Plotting the difference of the means of the 2 groups
plot(diff, col=c('red'),main='Scatterplot of the estimated difference of the means of the 2 groups', 
xlab='Number of genes',ylab='Estimated difference of the means of the 2 groups')

# Heatmap of the initial data matrix of the 2 groups
heatmap(new_mat)




## Question 2  ##
# PCA
set.pr <- princomp(scale(new_mat)) 
summary(set.pr)
# The 1st PC explains the 92.2% of the total variability
# The first 2 PCs explain the 94.4% of the total variability
# The first 3 PCs explain the 95.5% of the total variability
# The first 4 PCs explain the 96% of the total variability
# We observe that the standard deviations with value more than 1, is the first one ONLY


# Overall plot of all the Leukemia types using the first 2 PCs
plot(set.pr$scores[, 1:2], col = leukemiasEset$LeukemiaType, pch = 20)
legend("topleft", legend = unique(leukemiasEset$LeukemiaType), pch = 15, col = unique(as.numeric(leukemiasEset$LeukemiaType)))


# Plot of the 2 Leukemia types we're interested in, using the first 2 PCs
factorlvl <- gl(2, 12, labels = c("AML", "NOL"))
plot(set.pr$scores[, 1:2], col = factorlvl  , pch = 20)
legend("topleft", legend = unique(factorlvl), pch = 15, col = unique(factorlvl))
# We observe that there's a lot of overlap between the 2 groups





screeplot(set.pr,type="lines")
# The screeplot indicates to use the first 2 PCs as the greatest angle is located between the first and the second component

biplot(set.pr, choices=c(1,2) ,cex=.6)
# Biplot offers useful insights about the variables!



# Second way to perform PCA
library(factoextra)

# Compute PCA
res.pca <- prcomp(new_mat, scale = TRUE)
summary(res.pca)

# Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.
fviz_eig(res.pca)




# Graph of variables. Positive correlated variables point to the same side of the plot. Negative correlated variables point to opposite sides of the graph.
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
fviz_pca_var(res.pca,
             repel = TRUE     # Avoid text overlapping
)









# PCA on the transposed matrix of the new data
dim(t(new_mat))
# Check that we got the correct data

t_set.pr <- prcomp(scale(t(new_mat))) 
summary(t_set.pr)
# We observe that we need about 11 PCs to have 81.6% of the total variability explained


# Plot of the 2 Leukemia types we're interested in, using the first 2 PCs
plot(t_set.pr$x[,1:2] , col = factorlvl  , pch = 20)
legend("topleft", legend = unique(factorlvl), pch = 15, col = unique(factorlvl))
# We observe that with the transposed data the partitioning is much better (there's not as much overlap as before!)


# Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.
fviz_eig(t_set.pr)


biplot(t_set.pr, choices=c(1,2) ,cex=.6)
# Biplot of the transposed data now




# Graph of variables. Positive correlated variables point to the same side of the plot. Negative correlated variables point to opposite sides of the graph.
fviz_pca_var(t_set.pr,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(t_set.pr,
             repel = TRUE     # Avoid text overlapping
)













## Question 4 ##
# Once we’ve examined the distribution of the p-values the function qvalue can be used to calculate the q-values:
qobj <- qvalue(pval)
names(qobj)

# Summarizing results
summary(qobj)

# One very important statistic that is obtained with the software is an estimate of the overall proportion of true null hypotheses, ð0:
pi0 <- qobj$pi0
pi0
hist(qobj)


# The q-value is the minimum FDR incurred when calling a test significant. The q-values can be extracted from the qvalue object by:
qvalues <- qobj$qvalues



# Visualizing results
# The hist and plot functions can be used to visualize the results from qvalue. The function plot allows one to view several useful plots:
# The estimated ð0 versus the tuning parameter ë
# The q-values versus the p-values
# The number of significant tests versus each q-value cut-off
# The number of expected false positives versus the number of significant tests
# Applying plot to the qvalue object, we get:
plot(qobj)






## Question 5 ##
# BiocManager::install("qvalue")
library(qvalue)
library(VennDiagram)
m <- length(pval)
qStar <- 0.01

# Return a list of significant genes, while controlling FWER at q = 0.01
fwer_genes <- m * pval < qStar
table(fwer_genes)
# 92 genes are differentially expressed at q = 0.01 (based on Bonferroni procedure)

# Return a list of significant genes, while controlling FDR at q = 0.05
fdr_genes <- qvalue(p = pval, lambda = 0, fdr.level = qStar)$significant
table(fdr_genes)
# 953 genes are differentially expressed at q = 0.01 (based on BH procedure)

# Return a list of significant genes, while controlling FDR at q = 0.05
qvalue_genes <- qvalue(p = pval, fdr.level = qStar)$significant
table(qvalue_genes)
# 953 genes are differentially expressed at q = 0.01 (based on pFDR procedure)


# Venn them
n1 <- sum(fwer_genes)
n2 <- sum(fdr_genes)
n3 <- sum(qvalue_genes)

n12 <- sum(fwer_genes * fdr_genes)
n13 <- sum(fwer_genes * qvalue_genes)
n23 <- sum(fdr_genes * qvalue_genes)

n123 <- sum(fwer_genes * fdr_genes * qvalue_genes )

pdf(file = 'venn_type1errors.pdf', width = 6, height = 6)
grid.newpage()
draw.triple.venn(area1 = n1, area2 = n2, area3 = n3, n12 = n12, n23 = n23, n13 = n13, alpha = rep(0.7, 3),
                 n123 = n123, category = c("FWER", "FDR", "qValue"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"), cat.dist = c( -0.01, -0.02, -0.02), cat.cex = 2, cat.pos = 180, cex = 2)
dev.off()









## Question 6 ##
# Plotting the colored genes when controlling the FDR
plot(diff, estimated_var,  col=c('red','green'), main='Scatterplot of the differentially expressed genes', xlab='Difference of the means',
ylab='Estimated Variances')
legend('topleft',legend = c('Differentially Expressed','Not Differentially Expressed'),fill = c('red','green'))


plot(diff , col=as.factor(fdr_genes[fdr_genes==T]), main='Scatterplot of the differentially expressed genes', xlab='Number of genes',
ylab='Estimated Diferrence in the means of the 2 groups')

range(diff)



library(gplots)
fdr <- p.adjust(p=pval,method = 'fdr')
heat <- new_mat[fdr < 0.01,]
factorlvl <- gl(2, 12, labels = c("AML", "NOL"))
heatmap.2(heat, trace = "none", ColSideColors = as.character(as.numeric(factorlvl)), dendrogram = c('column'))
legend('bottomleft',legend = c('AML','NOL'),fill = c('red','black'))




plot(set.pr$scores[, 1:2], col = as.factor(fdr_genes), pch = 20,
main='Scatterplot of the differentially expressed genes using the first 2 PCs', xlab='PC1',
ylab='PC2')
legend("topleft", legend=c('AML','NOL') ,fill = c('red','black'))
# We observe that there's a lot of overlap between the 2 groups but the situation is clearly better than before











# Exercise 2
# Question 1
# Multiple testing simulation study
# Simulate the explanatory variables from independent normal distributions:
n <- 500
p <- 100
b <- numeric(p)
m <- 10000
ground.truth <- numeric(m)
pvals <- numeric(m)
for (j in 1:m) {
  x <- matrix(rnorm(n*p),nrow = n, ncol = p)
  # Generate the p regression coefficients b1, . . . , bp as follows:
  b <- numeric(p)
  if( runif(1) < 0.3){ b[1] <- rnorm(1) }
  # Generate the values of the response variable from a typical normal linear model, that is:
  y <- x%*%b + rnorm(n)
  ground.truth[j] <- ifelse(b[1]==0, 1, 2)
# This means that âi = 0 for all i > 2, while the first coefficient (b1) is zero with probability
# 0.7, while it is different than zero with probability 0.3.
  df <- data.frame(x)
  f <- lm(y~., data = df)
  pvals[j] <- pf(summary(f)$fstatistic[1], summary(f)$fstatistic[2], summary(f)$fstatistic[3], lower.tail=FALSE)
}
hist(pvals,30)
table(ground.truth) # 1 refers to true null hypothesis, while 2 refers to non-true null hypothesis

qstar <- 0.05
# Type I error rate target value

sortedP <- pvals[order(pvals)]
# Necessary for the Bonferroni procedure

bonfSelected <- sortedP < qstar/m 
table(bonfSelected)
# Based on Bonferroni procedure we will reject 2020 hypotheses and we will fail to reject 7980 hypotheses

# Now let's create a confusion matrix regarding the Bonferroni procedure
tt <- table(ground.truth[order(pvals)],bonfSelected)
tt
tt[2,2]/sum(tt[2,])
# Estimated power

# Second way of controlling FWER with Bonferroni procedure using the p.adjust command
adjustedP_Bonf <- p.adjust(p=pvals,method = 'bonf')
bonfSelected2 <- adjustedP_Bonf < qstar
table(bonfSelected2)
# As expected it's the same output as before

# Now let's create a confusion matrix regarding the Bonferroni procedure
tt2 <- table(ground.truth,bonfSelected2)
tt2
tt2[2,2]/sum(tt2[2,])
# Estimated power



# Controlling FWER with Holm procedure
adjustedP_holm <- p.adjust(p=pvals,method = 'holm')
holmSelected <- adjustedP_holm < qstar
table(holmSelected)


# Now let's create a confusion matrix regarding the Holm procedure
tt.h <- table(ground.truth,holmSelected)
tt.h
tt.h[2,2]/sum(tt.h[2,])
# Estimated power




# Controlling FWER with Hochberg procedure
adjustedP_hochberg <- p.adjust(p=pvals,method = 'hochberg')
hochbergSelected <- adjustedP_hochberg < qstar
table(hochbergSelected)


# Now let's create a confusion matrix regarding the Hochberg procedure
tt.hb <- table(ground.truth,hochbergSelected)
tt.hb
tt.hb[2,2]/sum(tt.hb[2,])
# Estimated power




# Controlling FWER with Hommel procedure
adjustedP_hommel <- p.adjust(p=pvals,method = 'hommel')
hommelSelected <- adjustedP_hommel < qstar
table(hommelSelected)


# Now let's create a confusion matrix regarding the Hommel procedure
tt.hml <- table(ground.truth,hommelSelected)
tt.hml
tt.hml[2,2]/sum(tt.hml[2,])
# Estimated power





# Controlling FDR with BH procedure
adjustedP_BH <- p.adjust(p=pvals,method = 'BH')
BHSelected <- adjustedP_BH < qstar
table(BHSelected)


# Now let's create a confusion matrix regarding the BH procedure
tt.BH <- table(ground.truth,BHSelected)
tt.BH
tt.BH[2,2]/sum(tt.BH[2,])
# Estimated power





# Controlling FDR with BY procedure
adjustedP_BY <- p.adjust(p=pvals,method = 'BY')
BYSelected <- adjustedP_BY < qstar
table(BYSelected)


# Now let's create a confusion matrix regarding the BY procedure
tt.BY <- table(ground.truth,BYSelected)
tt.BY
tt.BY[2,2]/sum(tt.BY[2,])
# Estimated power




# Controlling pFDR
pFDRSelected <- qvalue(p = pvals, fdr.level = qstar)$significant
table(pFDRSelected)


# Now let's create a confusion matrix regarding the pFDR procedure
tt.pFDR <- table(ground.truth,pFDRSelected)
tt.pFDR
tt.pFDR[2,2]/sum(tt.pFDR[2,])
# Estimated power










## Exercise 2 ##
##  Question 2 ##
alpha <- seq(0.001,0.999,0.01)
bonf <- numeric(length(alpha))
holm <- numeric(length(alpha))
hb <- numeric(length(alpha))
hml <- numeric(length(alpha))
bh <- numeric(length(alpha))
by <- numeric(length(alpha))
q <- numeric(length(alpha))


for (i in 1:length(alpha)) {
  adjustedP_Bonf <- p.adjust(p=pvals,method = 'bonf')
  bonfSelected <- adjustedP_Bonf < alpha[i]
  table.bonf <- table(ground.truth, bonfSelected)
  bonf[i] <- table.bonf[2,2]/sum(table.bonf[2,])
}



for (i in 1:length(alpha)) {
  adjustedP_holm <- p.adjust(p = pvals, method = 'holm')
  holmSelected <- adjustedP_holm < alpha[i]
  table.holm <- table(ground.truth, holmSelected)
  holm[i] <- table.holm[2,2]/sum(table.holm[2,])
}



for (i in 1:length(alpha)) {
  adjustedP_hb <- p.adjust(p = pvals, method = 'hochberg')
  hbSelected <- adjustedP_hb < alpha[i]
  table.hb <- table(ground.truth, hbSelected)
  hb[i] <- table.hb[2,2]/sum(table.hb[2,])
}



for (i in 1:length(alpha)) {
  adjustedP_hml<- p.adjust(p = pvals, method = 'hommel')
  hmlSelected <- adjustedP_hml < alpha[i]
  table.hml <- table(ground.truth, hmlSelected)
  hml[i] <- table.hml[2,2]/sum(table.hml[2,])
}



for (i in 1:length(alpha)) {
  adjustedP_BH <- p.adjust(p = pvals, method = 'BH')
  BHSelected <- adjustedP_BH < alpha[i]
  table.bh <- table(ground.truth, BHSelected)
  bh[i] <- table.bh[2,2]/sum(table.bh[2,])
}



for (i in 1:length(alpha)) {
  adjustedP_BY <- p.adjust(p = pvals, method = 'BY')
  BYSelected <- adjustedP_BY < alpha[i]
  table.by <- table(ground.truth, BYSelected)
  by[i] <- table.by[2,2]/sum(table.by[2,])
}




for (i in 1:length(alpha)) {
  q.value <- qvalue(p = pvals, fdr.level = alpha[i])$significant
  if (sum(qvalue(p = pvals, fdr.level = alpha[i])$significant)==10000) {
    q[i] <- 1} else {
      table.q <- table(ground.truth, q.value)
      q[i] <- table.q[2,2]/sum(table.q[2,])
    }
}





plot(alpha, bonf, type = "l", ylab = "Power", ylim=c(0,1), col='red')
lines(alpha, holm, type="l", col='blue')
lines(alpha, hb, type="l", col='green')
lines(alpha, hml, type="l", col='yellow')
lines(alpha, bh, type="l", col='orange')
lines(alpha, by, type="l", col='purple')
lines(alpha, q, type="l", col='cyan')
legend("bottomright", fill = c('red','blue','green','yellow','orange','purple','cyan'), legend = c("Bonferroni", "Holm", 
"Hochberg", "Hommel", "Benjamini & Çochberg", "Benjamini & Yekuteli", "q-Value"))




















