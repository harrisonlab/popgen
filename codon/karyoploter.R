library("karyoploteR")

for (contig in 1:34)
{
cname <- paste("contig_", contig, "_pilon", sep="")
custom.genome <- toGRanges(paste("chrom/",cname, sep=""))
genome_df <- read.delim(paste("chrom/",cname, sep=""), header=TRUE)
all.genes <- read.delim(paste("all_genes/",cname, sep=""), header=TRUE)
segmental <- read.delim(paste("segmental/",cname, sep=""), header=TRUE)
tandem <- read.delim(paste("tandem/",cname, sep=""), header=TRUE)

pdf(paste(cname, ".pdf", sep=""), width=6.95, height=4.27)
kp <- plotKaryotype(genome = custom.genome, plot.type=2, labels.plotter=NULL)
kpDataBackground(kp, data.panel = 1, r0=0, r1=0.8, color="lavenderblush")
kpDataBackground(kp, data.panel = 2, r0=0.05, r1=0.7, color="aliceblue")
#Vary x scale according to the contig size. 
#For small contigs, smaller than 1 million bp (contigs 15 and above), minor tick every 10,000 bp
#For big contigs, minor tick every 100,000 bp
if (contig > 14) {
  minor_tick=10000
  major_tick=1000000
}
else {
  minor_tick=100000
  major_tick=1000000
}

kpAddBaseNumbers(kp, col="gray50", tick.dist = major_tick, tick.len = 10, cex=0.6,
                 minor.tick.dist = minor_tick, minor.tick.len = 5)
#Custom x axis - eliminiated due to addition of kpAddBaseNumbers function
#l <- seq(0, genome_df$end, by = genome_df$end/10)
#lab <- pretty(l, n=10)
#axis(1, at=lab/genome_df$end, labels=format(lab, scientific=FALSE), cex.axis=0.7)


title(main=paste("contig", contig), col.main="gray65")

#generate counts of genes in given range
interval = 10000
number_of_intervals <- genome_df$end / interval
curr_interval = 0
start_r <- c()
end_r <- c()
obs_all <- c()
obs_tan <- c()
obs_seg <- c()

all_intervals <- 0:number_of_intervals
  for (i in all_intervals) {
  curr_interval = i
  min_bp = curr_interval * interval
  start_r <- c(start_r, min_bp)
  curr_interval = i + 1
  max_bp = curr_interval * interval
  if (max_bp > genome_df$end) {max_bp = genome_df$end}
  end_r <- c(end_r, max_bp)
  no_all <- sum(all.genes$end > min_bp & all.genes$end < max_bp)
  no_seg <- sum(segmental$end > min_bp & segmental$end < max_bp)
  no_tan <- sum(tandem$end > min_bp & tandem$end < max_bp)
  obs_all <- c(obs_all, no_all)
  obs_seg <- c(obs_seg, no_seg)
  obs_tan <- c(obs_tan, no_tan)
  }
#Calculate the remainder:
if (max_bp < genome_df$end) {
min_bp = max_bp
max_bp = genome_df$end
start_r <- c(start_r, min_bp)
end_r <- c(end_r, max_bp)
no_all <- sum(all.genes$end > min_bp & all.genes$end < max_bp)
no_seg <- sum(segmental$end > min_bp & segmental$end < max_bp)
no_tan <- sum(tandem$end > min_bp & tandem$end < max_bp)
obs_all <- c(obs_all, no_all)
obs_seg <- c(obs_seg, no_seg)
obs_tan <- c(obs_tan, no_tan)
}
#Create data frames
all_genes_input <- data.frame(start_r, end_r, obs_all)
cc <- rep(cname,length(obs_all))
all_genes_input <- cbind(cc, all_genes_input)

tandem_input <- data.frame(start_r, end_r, obs_tan)
cc <- rep(cname,length(obs_tan))
tandem_input <- cbind(cc, tandem_input)

segmental_input <- data.frame(start_r, end_r, obs_seg)
cc <- rep(cname,length(obs_seg))
segmental_input <- cbind(cc, segmental_input)

#plot all genes

#plot tandem duplications

#plot segmental duplications

#kpAxis(kp, ymin=0, ymax=max(all_genes_input$obs_all), r0=0, r1=0.8, col="gray50", cex=0.5, numticks=max(all_genes_input$obs_all)+1)
#Fixed at 10 genes
ymaximum=10
kpAxis(kp, ymin=0, ymax=ymaximum, r0=0, r1=0.8, col="gray50", cex=0.5, numticks=11)


kpLines(kp, chr=all_genes_input$cc, x=all_genes_input$start_r, y=all_genes_input$obs_all,
         col="#AA88FF", ymin=0, ymax=ymaximum, r0=0, r1=0.8)

kpLines(kp, chr=tandem_input$cc, x=tandem_input$start_r, y=tandem_input$obs_tan,
        col="black",  ymin=0, ymax=ymaximum, r0=0, r1=0.8)

kpLines(kp, chr=segmental_input$cc, x=segmental_input$start_r, y=segmental_input$obs_seg,
        col="brown", ymin=0, ymax=ymaximum, r0=0, r1=0.8)

#plot secreted effector P genes
effectorp <- read.delim(paste("effectorp_secreted/",cname, sep=""), header=TRUE)
y <- rep(1,length(effectorp$chr))
input_effectorP <- cbind(effectorp, y)

kpPoints(kp, chr=input_effectorP$chr, x=input_effectorP$start, y=input_effectorP$y,
         col="burlywood2", ymin=0, ymax=1, r0=0, r1=0.15, data.panel = 2, pch=".", cex=3.5)

#plot secreted CAZY
scazy <- read.delim(paste("cazy_secreted/",cname, sep=""), header=TRUE)
y <- rep(1,length(scazy$chr))
input_scazy <- cbind(scazy, y)

kpPoints(kp, chr=input_scazy$chr, x=input_scazy$start, y=input_scazy$y,
         col="forestgreen", ymin=0, ymax=1, r0=0.15, r1=0.3, data.panel = 2, pch=".", cex=3.5)

#plot secreted genes within 2 kbp of mimps
smimp <- read.delim(paste("mimp_secreted/",cname, sep=""), header=TRUE)
y <- rep(1,length(smimp$chr))
input_smimp <- cbind(smimp, y)

kpPoints(kp, chr=input_smimp$chr, x=input_smimp$start, y=input_smimp$y,
         col="slateblue3", ymin=0, ymax=1, r0=0.3, r1=0.45, data.panel = 2, pch=".", cex=3.5)

#plot secreted  metabolite gene clusters
smet <- read.delim(paste("metabolite_clusters/",cname, sep=""), header=TRUE)
y <- rep(1,length(smet$chr))
input_smet <- cbind(smet, y)

kpPoints(kp, chr=input_smet$chr, x=input_smet$start, y=input_smet$y,
         col="violetred", ymin=0, ymax=1, r0=0.45, r1=0.6, data.panel = 2, pch=".", cex=3.5)

#plot SIX genes
six <- read.delim(paste("six/",cname, sep=""), header=TRUE)
y <- rep(1,length(six$chr))
input_six <- cbind(six, y)
kpPoints(kp, chr=input_six$chr, x=input_six$start, y=input_six$y,
         col="black", ymin=0, ymax=1, r0=0.6, r1=0.66, data.panel = 2, pch=".", cex=5.5)

dev.off()
}
