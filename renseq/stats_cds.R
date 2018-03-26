
library(tools)
library("ggplot2")
args = commandArgs(trailingOnly=TRUE)
input_file1 <- args[1]
input_file2 <- args[2]
output <- args[3]
stats1 <- read.csv(input_file1, sep=" ", quote='', stringsAsFactors=TRUE,header=FALSE)
stats2 <- read.csv(input_file2, sep=" ", quote='', stringsAsFactors=TRUE,header=FALSE)
names(stats1)[1:6] <- c("query","hit","aln_length","qlen", "hlen", "coverage")
names(stats2)[1:6] <- c("query","hit","aln_length","qlen", "hlen", "coverage")

a <- data.frame(Cultivar = "Hapil", value = stats1$coverage, length=stats1$aln_length, qlen=stats1$qlen, hlen=stats1$hlen)
b <- data.frame(Cultivar = "RedGauntlet", value = stats2$coverage, length=stats2$aln_length, qlen=stats2$qlen, hlen=stats2$hlen)

# Combine into one long data frame
plot.data <- rbind(a, b)

##Coverage histogram
basic <- ggplot(data=plot.data, aes(plot.data$value,fill=Cultivar)) + geom_histogram(col="black", position='dodge') + scale_fill_manual(values=c("black", "red2")) + xlab("Hit coverage") + ylab("Number of RBB pairs") 
grids <- basic + theme(panel.grid.major = element_line(NA), panel.grid.minor = element_line(NA)) + scale_y_continuous(expand=c(0,0)) 
coverage_histogram <- grids + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title=element_text(size=12), legend.text=element_text(size=12)) + scale_x_continuous(breaks = pretty(plot.data$value, n = 10)) 
ggsave(paste(output, "_coverage_histogram.pdf", sep=""), coverage_histogram, dpi=300, height=5.7, width=8)

##CDS length histogram
basic <- ggplot(data=plot.data, aes(plot.data$length,fill=Cultivar)) + geom_histogram(col="black", position='dodge') + scale_fill_manual(values=c("black", "red2")) + xlab("Query CDS length") + ylab("Frequency") 
grids <- basic + theme(panel.grid.major = element_line(NA), panel.grid.minor = element_line(NA)) + scale_y_continuous(expand=c(0,0)) 
cds_length <- grids + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title=element_text(size=12), legend.text=element_text(size=12)) + scale_x_continuous(breaks = pretty(plot.data$length, n = 10)) 
ggsave(paste(output,"_cds_length.pdf", sep=""), cds_length, dpi=300, height=5.7, width=8)

# Dot plot (CDS length versus Query length)
basic <- ggplot(data=plot.data, aes(x=plot.data$value,y=plot.data$qlen,colour=Cultivar)) + geom_point() + scale_colour_manual(values=c("black", "red2")) + xlab("Hit coverage") + ylab("Query CDS length") 
grids <- basic + theme(panel.grid.major = element_line(NA), panel.grid.minor = element_line(NA)) + scale_y_continuous(expand=c(0,0), breaks = pretty(plot.data$qlen, n = 10)) 
dotplot1 <- grids + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title=element_text(size=12), legend.text=element_text(size=12)) + scale_x_continuous(breaks = pretty(plot.data$value, n = 10)) 
ggsave(paste(output,"_cds_coverage_query_length.pdf", sep=""), dotplot1, dpi=300, height=8, width=12)

# Dot plot (CDS length versus Hit length)
basic <- ggplot(data=plot.data, aes(x=plot.data$value,y=plot.data$hlen,colour=Cultivar)) + geom_point() + scale_colour_manual(values=c("black", "red2")) + xlab("Hit coverage") + ylab("Hit CDS length") 
grids <- basic + theme(panel.grid.major = element_line(NA), panel.grid.minor = element_line(NA)) + scale_y_continuous(expand=c(0,0), breaks = pretty(plot.data$qlen, n = 10)) 
dotplot2 <- grids + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title=element_text(size=12), legend.text=element_text(size=12)) + scale_x_continuous(breaks = pretty(plot.data$value, n = 10)) 
ggsave(paste(output, "_cds_coverage_hit_length.pdf", sep=""), dotplot2, dpi=300, height=8, width=12)
