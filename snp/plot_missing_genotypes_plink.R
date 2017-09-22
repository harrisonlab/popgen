library("ggplot2")
input_file <- "clean-inds-GWA-data.lmiss"
stats <- read.csv(input_file, sep="\t", quote='', stringsAsFactors=TRUE,header=TRUE)
basic <- ggplot(data=stats, aes(stats$F_MISS)) + geom_histogram(col="black", position='dodge') + scale_fill_manual(values=c("black")) + xlab("Fraction missing genotypes") + ylab("Number of SNPs") 
grids <- basic + theme(panel.grid.major = element_line(NA), panel.grid.minor = element_line(NA)) + scale_y_continuous(expand=c(0,0)) 
histogram_length <- grids + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title=element_text(size=12), legend.text=element_text(size=12)) + scale_x_continuous(breaks = pretty(stats$F_MISS, n = 15)) 
ggsave("clean-inds-GWA-data.pdf", histogram_length, dpi=300, height=5.7, width=9)
