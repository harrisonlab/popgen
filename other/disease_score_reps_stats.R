library("psych")
filename <- "EF_reps.txt"
my_data <- read.csv(filename, sep="\t", header=TRUE)
my_data2 <- my_data[,1:2]
my_data2$Genotype <- as.factor(my_data2$Genotype)
summary <- describeBy(my_data2, group=my_data2$Genotype, mat=TRUE)
write.table(summary, file="stats.txt", sep="\t")

my_DF <- data.frame(Genotype=character(), Mean=character(), Range=character()) 
for (i in levels(my_data2$Genotype)) 
{
  test <- my_data2[my_data2$Genotype == i,]
  Disease_Score <- sample(test$Disease_Score, z, replace=TRUE)
  my_mean <- mean(Disease_Score)
  range_val = max(Disease_Score) - min(Disease_Score)
  my_DF <- rbind(my_DF, data.frame(Genotype=i, Mean=my_mean, Range=range_val))
}
write.table(my_DF, file="original_data.txt", quote=FALSE)

for (z in 5:10)
{
my_file <- paste("my_sim_n=", z, ".txt", sep="")
my_DF <- data.frame(Genotype=character(), Mean=character(), Range=character()) 
n = 0
while (n < 100) 
{
for (i in levels(my_data2$Genotype)) 
{
test <- my_data2[my_data2$Genotype == i,]
Disease_Score <- sample(test$Disease_Score, z, replace=TRUE)
my_mean <- mean(Disease_Score)
range_val = max(Disease_Score) - min(Disease_Score)
my_DF <- rbind(my_DF, data.frame(Genotype=i, Mean=my_mean, Range=range_val))
}
n = n + 1
}
write.table(my_DF, file=my_file, quote=FALSE)
}

#Plot#
library(ggplot2)
for (k in 5:9)
{
filename <- paste("residuals_", k, ".txt", sep="")
my_data <- read.table(filename, header=TRUE)
my_data$Genotype <- as.factor(my_data$Genotype)
my_title <- paste("Reps = ", k, sep="")
mean_plot <- ggplot(my_data, aes(x=Genotype,y=Mean)) + geom_boxplot(fill="lightblue", coef = 6) + ggtitle(my_title)
range_plot <- ggplot(my_data, aes(x=Genotype,y=Range)) + geom_boxplot(fill="pink", coef = 6) + ggtitle(my_title)
my_mean <- paste("r_", k, "_mean.pdf", sep="")
ggsave(my_mean, mean_plot, width = 50, height = 6, limitsize=FALSE)
my_range <- paste("r_", k, "_range.pdf", sep="")
ggsave(my_range, range_plot, width = 50, height = 6, limitsize=FALSE)
}