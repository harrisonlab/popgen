library("ggplot2")
input_file <- "all_everbearer_kinship.genome"
stats <- read.csv(input_file, sep="\t", quote='', stringsAsFactors=TRUE,header=TRUE)
##Plot Pi hat histogram
basic <- ggplot(data=stats, aes(stats$PI_HAT)) + geom_histogram(col="black", position='dodge') + scale_fill_manual(values=c("black")) + xlab("PI_HAT") + ylab("Frequency") 
grids <- basic + theme(panel.grid.major = element_line(NA), panel.grid.minor = element_line(NA)) + scale_y_continuous(expand=c(0,0)) 
histogram_length <- grids + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title=element_text(size=12), legend.text=element_text(size=12)) + scale_x_continuous(breaks = pretty(stats$PI_HAT, n = 10)) 
ggsave("Pihat.pdf", histogram_length, dpi=300, height=5.7, width=8)

##Plot DST histogram
basic <- ggplot(data=stats, aes(stats$DST)) + geom_histogram(col="black", position='dodge') + scale_fill_manual(values=c("black")) + xlab("DST") + ylab("Frequency") 
grids <- basic + theme(panel.grid.major = element_line(NA), panel.grid.minor = element_line(NA)) + scale_y_continuous(expand=c(0,0)) 
histogram_length <- grids + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title=element_text(size=12), legend.text=element_text(size=12)) + scale_x_continuous(breaks = pretty(stats$DST, n = 10)) 
ggsave("DST.pdf", histogram_length, dpi=300, height=5.7, width=8)

#Plot DST versus PIHAT scatterplot
basic <- ggplot(data=stats, aes(x=stats$DST, y=stats$PI_HAT)) + geom_point(col="black", position='dodge') + scale_fill_manual(values=c("black")) + xlab("DST") + ylab("PI_HAT") 
grids <- basic + theme(panel.grid.major = element_line(NA), panel.grid.minor = element_line(NA)) + scale_y_continuous(expand=c(0,0)) 
histogram_length <- grids + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title=element_text(size=12), legend.text=element_text(size=12)) + scale_x_continuous(breaks = pretty(stats$DST, n = 10)) 
ggsave("DST_ver_Pihat.pdf", histogram_length, dpi=300, height=5.7, width=8)

#sample table
input_file <- "sample"
sample_table <- read.csv(input_file, sep="\t", quote='', stringsAsFactors=TRUE, header=TRUE)
sample_table <- sample_table[,1:2]

merge1 <- merge(stats, sample_table, by.x="IID1", by.y="id")
colnames(merge1)[colnames(merge1)=="clone_id"] <- "IID1_clone_id"
merge2 <- merge(merge1, sample_table, by.x="IID2", by.y="id")
colnames(merge2)[colnames(merge2)=="clone_id"] <- "IID2_clone_id"
write.table(merge2, file="all_everbearer_samples_kinship.genome", sep="\t", quote=FALSE, row.names=FALSE)