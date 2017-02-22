library(pcadapt)
library(qvalue)

args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Two arguments must be supplied (input file)", call.=FALSE)}

input_file <- argv[1]
#Choose ploidy (1 or 2)
ploidy_choice <- argv[2]

filename <- read.pcadapt(input_file,type="vcfR",ploidy=ploidy_choice)
x <- pcadapt(filename,K=20,ploidy=ploidy_choice)

plot(x,option="screeplot")

# With names
poplist.names <- c(rep("POP1",50),rep("POP2",50),rep("POP3",50))
plot(x,option="scores",pop=poplist.names)

plot(x,option="scores",i=3,j=4,pop=poplist.names)

x <- pcadapt(filename,K=2)

summary(x)

plot(x,option="manhattan")
plot(x,option="qqplot",threshold=0.1)

hist(x$pvalues,xlab="p-values",main=NULL,breaks=50)
plot(x,option="stat.distribution")

qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliers <- which(qval<alpha)

snp_pc <- get.pc(x,outliers)