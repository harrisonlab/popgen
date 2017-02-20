#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file)", call.=FALSE)}

input <- args[1]

source('/home/sobczm/bin/bayescan2.1/R_functions/plot_R.r')