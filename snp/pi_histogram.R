library("ggplot2")
all_sites_istraw90 <- read.csv("nucleotide_diversity.all.sites.pi.istraw90", sep="\t", quote='', stringsAsFactors=TRUE,header=TRUE)
all_sites_all <- read.csv("nucleotide_diversity.all.sites.pi", sep="\t", quote='', stringsAsFactors=TRUE,header=TRUE)

a <- data.frame(group = "all istraw90 sites", value = all_sites_istraw90$PI)
b <- data.frame(group = "istraw35 and istraw90 shared", value = all_sites_all$PI)

# Combine into one long data frame
plot.data <- rbind(a, b)
fig <- ggplot(data=plot.data, aes(plot.data$value, fill=plot.data$group)) + geom_histogram(col="black") + xlab("Pi") + labs(fill="Marker subset")
ggsave("pi_histogram.pdf", fig, dpi=300, height=5.7, width=8)