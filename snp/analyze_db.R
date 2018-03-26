library("ggplot2")
input_file <- "qc_table"
stats <- read.csv(input_file, sep="\t", quote='', stringsAsFactors=TRUE,header=TRUE)
overview <- summary(stats)
####pipeline1_number_istraw35_genotypes
## Number of genotypes : Number of markers
## 0 (marker has a defined marker id and Affx code on the chip but no genotypes available) 
## NA (marker DOES NOT have a defined markerd id and Affx code on the chip - missing from the chip entirely)
# 0    :11964                                                                          
# 379  :20613                                                                 
# 758  :1683 <- counted twice, as two Affx codes available for the marker 
# # NA   :69649

####pipeline2_number_istraw35_genotypes
#0  :34260  
#NA :69649 

####pipeline1_number_istraw90_genotypes
#0  :103909

####pipeline2_number_istraw90_genotypes
#Plot it               
histogram <- ggplot(data=stats, aes(stats$pipeline2_number_istraw90_genotypes)) + geom_histogram(col="black", position='dodge') + scale_fill_manual(values=c("black")) + xlab("Number of genotypes") + ylab("Number of markers") + ggtitle("Istraw90 - Pipeline 2")
ggsave("pipeline2_istraw90_markers.png", histogram)
#Conclusion
#Pipeline 1 - istraw35 samples
#Pipeline 2 - istraw90 samples

