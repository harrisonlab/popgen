## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence 
## interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95) {
  library(doBy)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # Collapse the data
  formula <- as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "))
  datac <- summaryBy(formula, data=data, FUN=c(length2,mean,sd), na.rm=na.rm)
  
  # Rename columns
  names(datac)[ names(datac) == paste(measurevar, ".mean",    sep="") ] <- measurevar
  names(datac)[ names(datac) == paste(measurevar, ".sd",      sep="") ] <- "sd"
  names(datac)[ names(datac) == paste(measurevar, ".length2", sep="") ] <- "N"
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#Read in table
a <- Sys.glob("*.txt")
results <- read.delim(a)

#Rename the samples
library(plyr)
results <- rename(results, c("PRO1467_S1"="Emily_0_1", "PRO1467_S2"="Emily_0_2", "PRO1467_S3"="Emily_0_3", 
                  "PRO1467_S10"="Fenella_0_1","PRO1467_S11"="Fenella_0_2", "PRO1467_S12"="Fenella_0_3",
                  "PRO1467_S4"="Emily_12_1","PRO1467_S5"="Emily_12_2", "PRO1467_S6"="Emily_12_3",
                  "PRO1467_S14"="Fenella_12_1", "PRO1467_S15"="Fenella_12_2", "PRO1467_S22"="Fenella_12_3",
                  "PRO1467_S7"="Emily_48_1", "PRO1467_S20"="Emily_48_2", "PRO1467_S9"="Emily_48_3",
                  "PRO1467_S17"="Fenella_48_1", "PRO1467_S18"="Fenella_48_2", "PRO1467_S28"="Fenella_48_3"))
#drop R.Tpe and QTL columns
results <- results[ -c(2:3) ]
library(reshape2)
results_long <- melt(results,
                 id.vars = "Gene",
                 measure.vars = c("Emily_0_1", "Emily_0_2", "Emily_0_3", 
                                   "Fenella_0_1", "Fenella_0_2", "Fenella_0_3",
                                  "Emily_12_1", "Emily_12_2", "Emily_12_3",
                                  "Fenella_12_1", "Fenella_12_2", "Fenella_12_3",
                                  "Emily_48_1", "Emily_48_2", "Emily_48_3",
                                  "Fenella_48_1", "Fenella_48_2", "Fenella_48_3"),
                 variable.name = "sample")
head(results_long)

Timepoint <- factor(rep(NA, length(results$sample) ), 
                 levels=c("0h", "12h", "48h") )   
results_long$Timepoint[ results_long$sample %in% c("Emily_0_1", "Emily_0_2", "Emily_0_3", 
                                         "Fenella_0_1", "Fenella_0_2", "Fenella_0_3")] <- "0h"
results_long$Timepoint[ results_long$sample %in% c("Emily_12_1", "Emily_12_2", "Emily_12_3", 
                                                   "Fenella_12_1", "Fenella_12_2", "Fenella_12_3")] <- "12h"
results_long$Timepoint[ results_long$sample %in% c("Emily_48_1", "Emily_48_2", "Emily_48_3", 
                                                   "Fenella_48_1", "Fenella_48_2", "Fenella_48_3")] <- "48h"
Cultivar <- factor(rep(NA, length(results$sample) ), 
                    levels=c("Fenella", "Emily") )   
results_long$Cultivar[ results_long$sample %in% c("Emily_0_1", "Emily_0_2", "Emily_0_3", 
                                                  "Emily_12_1", "Emily_12_2", "Emily_12_3",
                                                  "Emily_48_1", "Emily_48_2", "Emily_48_3")] <- "Emily"
results_long$Cultivar[ results_long$sample %in% c("Fenella_0_1", "Fenella_0_2", "Fenella_0_3", 
                                                  "Fenella_12_1", "Fenella_12_2", "Fenella_12_3",
                                                  "Fenella_48_1", "Fenella_48_2", "Fenella_48_3")] <- "Fenella"



library(ggplot2)
results_long_sum <- summarySE(results_long, measurevar="value", groupvars=c("Timepoint","Cultivar", "Gene"))
for(i in seq_along(levels(results_long_sum$Gene))) 
{ 
  val <- levels(results_long_sum$Gene)[i] 
  dfc <- subset(results_long_sum, Gene == val) 
  outfile <- paste(val, ".pdf", sep="")
  barchart <- ggplot(dfc, aes(x=Timepoint, y=value, fill=Cultivar)) + ylab("FPKM") +
    geom_bar(position=position_dodge(), stat="identity") +
    scale_fill_manual(values=c("cornflowerblue", "lightcoral")) +
    geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.2, position=position_dodge(.9))                  
    ggsave(outfile, barchart, width = 6, height = 4)
  
}

library(pheatmap)
results2 <- results[,-1]
rownames(results2) <- results[,1]
res_matrix <- data.matrix(results2, rownames.force = "True")
res_matrix[res_matrix == 0] <- 1
res_matrix_log <- log(res_matrix)
pheatmap(res_matrix_log, cluster_col=FALSE)