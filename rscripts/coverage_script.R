#!/usr/bin/env Rscript
# Generate coverage plots
#
# $1 rscript Folder
# $2 results folder 
# $3 x low
# $4 x high
# $5 Coverage output folder

library(getopt)
spec = matrix(c(
    'coverage_file', 'c', 1, "character",
    'reference_fasta', 'r', 2, "character",
    'sample_name','s',2,"character",
    'output_prefix','o' , 2, "character",
    'directory_output', 'd' ,2,"character",
    'total_coverage_file', 't', 2, 'character',
    'contamination_mapping','C', 2, 'character',
    'help','h',0,'logical'
    ),byrow=TRUE, ncol=4)
                
opt=getopt(spec)
if (!is.null(opt$help)){
    cat(getopt(spec, usage=TRUE));
    q(status=1);
}
if(!is.null(opt$total_coverage_file)){
    if(!file.exists(opt$total_coverage_file)){
        file.create(opt$total_coverage_file)
    }
    opt$total_coverage_file=normalizePath(opt$total_coverage_file)
}
opt$coverage_file <- normalizePath(opt$coverage_file) 
if(!is.null(opt$directory_output)){
    setwd(opt$directory_output)
}

#TODO improve so we don't use a hardcoded path for the script
source('~/Programming/OpenSource/MyGitHub/ancient_dna_pipeline/rscripts/coverage_plot.R')
coverage_plot(opt$coverage_file,opt$reference_fasta, opt$sample_name, opt$output_prefix,opt$total_coverage_file,opt$contamination_mapping)
