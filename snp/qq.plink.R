#!/usr/bin/Rscript

# Usage qq.plink.R gwas.txt "A plot title"
# Note: Required that columns name with P values is labeled as 'P'.
library(tools)
args <- commandArgs(trailingOnly = TRUE)

plink.f = args[1]
title = args[2]
qq.f = paste(plink.f, ".qq.pdf", sep="")

message("Input file = ", plink.f)
message("Title = ", title)
message("Output file = ", qq.f)

message("Reading data ...")

d <- tryCatch({ readRDS(plink.f) },
              warning = function(war) { read.table(plink.f, header=TRUE) },
              error = function(err) { read.table(plink.f, header=TRUE) },
              finally = {
                  message("")
              }
              )

do <- d

obs <- sort(d$P)
obs <- obs[!is.na(obs)]
obs <- obs[is.finite(-log10(obs))]

lobs <- -(log10(obs))

exp <- c(1:length(obs)) 
lexp <- -(log10(exp / (length(exp)+1)))

chisq <- qchisq(1-obs,1)

lambda <- round(median(chisq)/qchisq(0.5,1),2)

message("Drawing Plot")
pdf(file=qq.f, height = 10, width = 10)
par(mar=c(5,5,4,2))
plot(0, 0, type="n",	
	xlim=c(0,max(lexp)), ylim=c(0,max(lobs)), 
	xlab=expression("Expected (-log"[10]*"(P))"),
	ylab=expression("Observed (-log"[10]*"(P))"),
	#las=1, xaxs="i", yaxs="i", bty="l",
	bty='n', 
	main=title)
legend("topleft", c(as.expression(bquote(lambda == .(lambda))),
		  as.expression(bquote(N == .(length(lobs))))),	
		  , bty='n')
abline(0,1,col="red",lwd=3, lty=3)
points(lexp, lobs, cex=1, pch=21, bg="#0000ff33", col="#0000ff")
dev.off()
message("Lambda = ", lambda)
