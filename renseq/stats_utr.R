library(tools)
library("ggplot2")
args = commandArgs(trailingOnly=TRUE)
input_file1 <- args[1]
input_file2 <- args[2]
output <- args[3]
utr <- args[4]
stats1 <- read.csv(input_file1, sep="\t", quote='', stringsAsFactors=TRUE,header=FALSE)
stats2 <- read.csv(input_file2, sep="\t", quote='', stringsAsFactors=TRUE,header=FALSE)
names(stats1)[1:5] <- c("query","hit","utr_length","qlen", "hlen")
names(stats2)[1:5] <- c("query","hit","utr_length","qlen", "hlen")

a <- data.frame(Cultivar = "Hapil", length=stats1$utr_length, qlen=stats1$qlen, hlen=stats1$hlen)
b <- data.frame(Cultivar = "RedGauntlet", length=stats2$utr_length, qlen=stats2$qlen, hlen=stats2$hlen)

# Combine into one long data frame
plot.data <- rbind(a, b)

##UTR length
basic <- ggplot(data=plot.data, aes(plot.data$length,fill=Cultivar)) + geom_histogram(col="black", position='dodge') + scale_fill_manual(values=c("black", "red2")) + xlab(paste("Query ", utr, "' UTR length", sep="")) + ylab("Frequency") 
grids <- basic + theme(panel.grid.major = element_line(NA), panel.grid.minor = element_line(NA)) + scale_y_continuous(expand=c(0,0)) 
utr_length <- grids + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title=element_text(size=12), legend.text=element_text(size=12)) + scale_x_continuous(breaks = pretty(plot.data$length, n = 10)) 
ggsave(paste(output, utr, "utr_length.pdf", sep=""), utr_length, dpi=300, height=5.7, width=8)

# Dot plot (UTR length versus Query length)
basic <- ggplot(data=plot.data, aes(x=plot.data$length,y=plot.data$qlen,colour=Cultivar)) + geom_point() + scale_colour_manual(values=c("black", "red2")) + xlab("Hit coverage") + ylab(paste("Query ", utr,  "' UTR length")) 
grids <- basic + theme(panel.grid.major = element_line(NA), panel.grid.minor = element_line(NA)) + scale_y_continuous(expand=c(0,0), breaks = pretty(plot.data$qlen, n = 10)) 
dotplot1 <- grids + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title=element_text(size=12), legend.text=element_text(size=12)) + scale_x_continuous(breaks = pretty(plot.data$length, n = 10)) 
ggsave(paste(output, utr, "utr_length_query_length.pdf", sep=""), dotplot1, dpi=300, height=8, width=12)

# Dot plot (UTR length versus Hit length)
basic <- ggplot(data=plot.data, aes(x=plot.data$length,y=plot.data$hlen,colour=Cultivar)) + geom_point() + scale_colour_manual(values=c("black", "red2")) + xlab("Hit coverage") + ylab(paste("Hit ", utr,  "' UTR length"))
grids <- basic + theme(panel.grid.major = element_line(NA), panel.grid.minor = element_line(NA)) + scale_y_continuous(expand=c(0,0), breaks = pretty(plot.data$hlen, n = 10)) 
dotplot2 <- grids + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title=element_text(size=12), legend.text=element_text(size=12)) + scale_x_continuous(breaks = pretty(plot.data$length, n = 10)) 
ggsave(paste(output, utr, "utr_length_hit_length.pdf", sep=""), dotplot2, dpi=300, height=8, width=12)