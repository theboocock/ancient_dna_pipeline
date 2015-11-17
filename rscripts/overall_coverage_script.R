#!/usr/bin/env/Rscript
#
# This script generates an overall coverage plot and individual coverag# plots for each file, currently redundant slightly with the overall script
#
#
# @date 17 November 2015.
# @author James Boocock
#
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
coverage_files = arguments$args

contamination_mapping = opt$contamination_mapping
reference_fasta = opt$reference_fasta
output_directory = opt$output_directory
overall_coverage_file = opt$overall_coverage_file 
if (!is.null(reference_fasta)){
    referenceGenome = readDNAStringSet(reference_fasta, format="fasta")
    if(!is.null(contamination_mapping)){
        one_to_use =  grep(contamination_mapping,names(referenceGenome),fixed=T)
        xlim_coord = c(1,length(referenceGenome[[one_to_use]]))
    }
}
cov_df = data.frame()
cov_top = data.frame()
for (cov_file in coverage_files){
    cov_tmp=rep(NA,total_regions)
    cov = (read.table(cov_file,header=F,sep='\t'))
    for(j in 1:total_regions){
      if (j %in% cov[,2]){
        cov_tmp[j] = cov[which(cov[,2]==j),4]
      } else {
        cov_tmp[j] = 0
      }
    }
     name = strsplit(basename(cov_files[[i]]),split = ".",fixed =T)[[1]][1]
    t = (cut(cov_tmp,breaks = c(seq(0,19,by=1),10000),right = T))
    tab = table(t)/total_regions
    ggplot_cov = data.frame(read_depth=c(seq(1,19,by=1),">=20"),tabs=tab)
    ggplot_cov = cbind.data.frame(name,ggplot_cov)
    cov_df = rbind.data.frame(cov_df, ggplot_cov)
    cov_top = rbind.data.frame(cov_top, cov_tmp)
} 
 cov_df[,2] = factor(cov_df[,2],levels=c(seq(1,19,by=1),">=20"))
  #cov_df %>% ggplot(aes(x=name)) + geom_bar(stat='identity',aes(y=tabs.Freq)) + coord_flip() + 
  #  xlab("Sample ID")   + ylab("Coverage Proportion")
  wide_cov = dcast(cov_df, name ~ read_depth)
  sorted_wide_cov= wide_cov[(order(apply(wide_cov[,c(2:21)],1,sum))),]
  sorted_wide_cov[,1] = factor(sorted_wide_cov[,1])
  #sorted_wide_cov[,1] = reorder(sorted_wide_cov[,1],do.call(order,wide_cov[,c(2:21)]))
  p = as.character((sorted_wide_cov[,1]))
  w_cov = melt(sorted_wide_cov, id=("name"))
  w_cov[,1] = factor(w_cov[,1],levels=p)
  w_cov = cbind.data.frame(w_cov)
  names(w_cov) = c("name","read_depth","value")
  png(overall_coverage_file)
  w_cov%>% ggplot(aes(x=name, fill=read_depth)) + geom_bar(stat='identity',aes(y=value)) + 
  coord_flip() + xlab("Sample") + ylab("Proportion of bases covered by > 2 reads") + ggtitle("Mitochondrial Coverage") + scale_fill_discrete(name="Read Depth") +  theme_bw()
  dev.off()
  cov_tran = t(cov_top)
  colnames(cov_tran)= unique(cov_df$name)
  calculate_mean_sd = cov_tran
  meanien = apply(calculate_mean_sd,2,mean)
  sds = apply(calculate_mean_sd,2,sd)
  cov_tran = cov_tran %>% melt()
  r = rep(seq(1:total_regions),length(unique(w_cov[,1])))
  cov_tran = cbind(r,cov_tran)
  cov_tran = cov_tran[,-2]
  names(cov_tran) = c("pos","name","depth")
  column_names = unique(cov_tran$name)
    
for (i in 1:length(coverage_file)){
    # Extract  sample name 

    temp_map_q = cov_tran[which(cov_tran$name %in% column_names[i]),]
    png(paste0(directory,"/",column_names[i],"_cov.png"))
    p = temp_map_q %>% ggplot(aes(x=pos,y=depth)) + geom_area() + theme_bw() + 
    theme(legend.position="none") + scale_x_continuous(breaks=pretty_breaks(10))  + 
    xlab("Genome Position") + ylab("Read Depth") + ggtitle(column_names[i])
    dev.off()
} 
  





