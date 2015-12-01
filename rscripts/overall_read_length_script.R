#!/usr/bin/env/Rscript
#
# This script generates an read-length plots for each file
#
#
#
# @date 17 November 2015.
# @author James Boocock
#
library(optparse)
library(gridExtra)
library(Biostrings)
library(reshape2)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Rsamtools)
library(scales)

option_list = list(
    make_option(c("-r","--reference_fasta")),
    make_option(c("-c","--contamination_mapping")),
    make_option(c("-o","--overall_coverage_file")),
    make_option(c("-d","--directory"))
    )
parser = OptionParser(usage = "%prog [options] coverage_files", option_list=option_list)
arguments = parse_args(parser, positional_arguments = 1)
opt = arguments$options
bam_files =arguments$args

contamination_mapping = opt$contamination_mapping
reference_fasta = opt$reference_fasta
output_directory = opt$output_directory
overall_coverage_file = opt$overall_coverage_file 

m_d = data.frame()
map_q = data.frame()
length_c_t = data.frame()
r_names = c()
r_names = c()
MAP_FILTER=20
for(i in seq(1,length(bams_files),by=1)){
# Skip MS10129http://127.0.0.1:21580/graphics/plot_zoom_png?width=1028&height=859
bam = scanBam(bam_files[[i]])

seq_length= cut(width(bam[[1]]$seq),breaks=seq(from=25,to=150,by=10))
mapq_length = (bam[[1]]$mapq)
m = split(mapq_length ,seq_length)
m2 = split(cbind(bam[[1]]$seq,bam[[1]]$mapq,bam[[1]]$pos), seq_length)

m_percent = lapply(m,function(x){sum(!is.na(x) & x > MAP_FILTER)/length(x) })
l = apply(cbind(bam[[1]]$pos, width(bam[[1]]$seq)),1, function(x){
  if(!is.na(x[1])){
    get_region_from_reference(dog_output,x[1],x[1]+x[2])
  }else{
    NA
  }
  
})

mapq = width(bam[[1]]$seq[which(bam[[1]]$mapq > MAP_FILTER)])
lapply(m2, function(x){  })
# First step of the analysis.
name = strsplit(basename(bam_files[[i]]),split = ".",fixed =T)[[1]][1]
mapq = cbind.data.frame(name,(mapq))
m_new = c(name,as.numeric(unlist(m_percent)))
m_d = rbind(m_d, m_new)
map_q = rbind.data.frame(map_q, mapq)
}
}
m_d[,2:ncol(m_d)]  = apply(m_d[,2:ncol(m_d)],2,as.numeric)
colnames(m_d) = c("name",as.numeric(seq(from=25,to=139,by=10)))
m_d = m_d %>% gather(read_length, sample, -name)
m_d[,2] = as.numeric(as.character(m_d[,2]))
colnames(map_q) = c("name", "read_length")
p = map_q %>% ggplot(aes(read_length)) + geom_density(fill='grey') + facet_wrap(~name,ncol=4,nrow=11) + theme_bw() +  theme(legend.position="none") + 
scale_x_continuous(breaks=pretty_breaks(10))  + xlab("Read Length") +
ylab("Mapped read density") + ggtitle("Merged Read Length Distributions")
