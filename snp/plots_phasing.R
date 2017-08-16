library(ggplot2)
library(tools)
#aggregate data from all 28 linkage groups -> default
args = commandArgs(trailingOnly=TRUE)
my_comparison <- args[1]
dir_original <- args[2]
dir_comparison <- args[3]

setwd(dir_original)
df <- NULL
for(hg in 1:7)
{
  for(sg in c('A','B','C','D'))
  {
    fname <- paste0(hg,sg,"_out")
    tmp <- read.csv(fname,sep=" ",header=F)
    colnames(tmp) <- c("sample","same","different","switches","inconsistent")
    tmp$lg <- paste0(hg,sg)
    tmp$type <- "default"
    df <- rbind(df,tmp)
    }
}

setwd(dir_comparison)
for(hg in 1:7)
{
  for(sg in c('A','B','C','D'))
  {
      fname <- paste0(hg,sg,"_out")
      tmp <- read.csv(fname,sep=" ",header=F)
      colnames(tmp) <- c("sample","same","different","switches","inconsistent")
      tmp$lg <- paste0(hg,sg)
      tmp$type <- "tested"
      df <- rbind(df,tmp)
  }
}

df$lg <- as.factor(df$lg)
df$type <- as.factor(df$type)
df$sample <- as.factor(df$sample)

my_plot <- ggplot(df,aes(x=switches,fill=type)) +
  geom_histogram(binwidth=1,position='identity',alpha=0.5) +
  xlab("haplotyping errors (switches) per progeny") +
  ggtitle(paste("shapeit haplotype imputation errors on rgxha map: default versus ", my_comparison))
outfile <- paste(my_comparison, ".pdf", sep="")
ggsave(outfile, my_plot, width = 8, height = 5)
