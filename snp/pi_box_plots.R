library("ggplot2")
all_sites <- read.csv("nucleotide_diversity.all.sites.pi", sep="\t", quote='', stringsAsFactors=TRUE,header=TRUE)
syn_sites <- read.csv("nucleotide_diversity.syn.sites.pi", sep="\t", quote='', stringsAsFactors=TRUE,header=TRUE)
nonsyn_sites <- read.csv("nucleotide_diversity.nonsyn.sites.pi", sep="\t", quote='', stringsAsFactors=TRUE,header=TRUE)
fourfold_sites <- read.csv("nucleotide_diversity.syn4fd.sites.pi", sep="\t", quote='', stringsAsFactors=TRUE,header=TRUE)

a <- data.frame(group = "all sites", value = all_sites$PI)
b <- data.frame(group = "syn sites", value = syn_sites$PI)
c <- data.frame(group = "4-fold syn sites", value = fourfold_sites$PI)
d <- data.frame(group = "nonsyn sites", value = nonsyn_sites$PI)

# Combine into one long data frame
plot.data <- rbind(a, b, c, d)
fig <- ggplot(plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot() + ylab("Pi") + theme(axis.title.x=element_blank())
ggsave("pi_boxplot.pdf", fig, dpi=300, height=5.7, width=8)