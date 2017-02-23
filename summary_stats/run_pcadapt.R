library(pcadapt)
library(qvalue)
library(tools)

args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("At least two arguments must be supplied (input file and ploidy)", call.=FALSE)}

###Warning! Following the inspection of scree plot after the first run, 
#the script then requires manual changing of the final_k parameter 
#which is the k value best explaining the population structure.

#Input VCF file
input_file <- args[1]
#Choose ploidy (1 or 2)
ploidy_choice <- args[2]

###Script modelled after tutorial: https://cran.r-project.org/web/packages/pcadapt/vignettes/pcadapt.html

filename <- read.pcadapt(input_file,type="vcf",ploidy=ploidy_choice)
#Check how many individuals present in the file. 
#K used below has to be a maximum number_of_individuals - 1
input_pca <- paste(file_path_sans_ext(input_file), ".pcadapt", sep="")
pca_table <- read.delim(input_pca, header=FALSE, sep=" ")
k_value <- ncol(pca_table) - 1
x <- pcadapt(filename,K=k_value,ploidy=ploidy_choice)

#Choosing the number K of Principal Components
#The ‘scree plot’ displays in decreasing order the percentage of variance explained by each PC.
#The ideal pattern in a scree plot is a steep curve followed by a bend and a straight line. 
#The eigenvalues that correspond to random variation lie on a straight line
#whereas the ones that correspond to population structure lie on a steep curve.
#We recommend to keep PCs that correspond to eigenvalues to the left of the straight line. 
out_scree <- paste(file_path_sans_ext(input_file), "_scree.pdf", sep="")
pdf(out_scree, width=11, height=8)
plot(x,option="screeplot")
dev.off()

#Computing the test statistic based on PCA
#!!!Change the final_k value below:
final_k <- k_value
x <- pcadapt(filename,K=final_k)

summary(x)

#A Manhattan plot displays −log10−log10 of the p-values.
out_man <- paste(file_path_sans_ext(input_file), "_manhattan.pdf", sep="")
pdf(out_man, width=11, height=8)
plot(x,option="manhattan")
dev.off()

#The user is also given the possibility to check the distribution of the p-values using a Q-Q plot
out_qq <- paste(file_path_sans_ext(input_file), "_qq.pdf", sep="")
pdf(out_qq, width=11, height=8)
plot(x,option="qqplot",threshold=0.1)
dev.off()

#A histogram of p-values
out_hist <- paste(file_path_sans_ext(input_file), "_pval.pdf", sep="")
pdf(out_hist, width=11, height=8)
hist(x$pvalues,xlab="p-values",main=NULL,breaks=50)
dev.off()

#Plot a histogram of the test statistic Dj
out_dj <- paste(file_path_sans_ext(input_file), "_dj.pdf", sep="")
pdf(out_dj, width=11, height=8)
plot(x,option="stat.distribution")
dev.off()

####Choosing a cutoff for outlier detection
#For a given αα (real valued number between 0 and 1), SNPs with q-values less than α
#will be considered as outliers with an expected false discovery rate bounded by α.
#The false discovery rate is defined as the percentage of false discoveries among the list of candidate SNPs. 
#Here is an example of how to provide a list of candidate SNPs for the geno3pops data, 
#for an expected false discovery rate lower than 10%:
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliers <- which(qval<alpha)

#It may be interesting to know which principal components are actually 
#the most correlated with the oulier SNPs. The function get.pc allows to achieve that:
snp_pc <- get.pc(x,outliers)