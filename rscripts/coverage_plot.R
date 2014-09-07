#http://learnr.wordpress.com/2009/04/29/ggplot2-labelling-data-series-and-adding-a-data-table/


coverage_plot=function(cov,reference_fasta=NULL, sample_name=NULL,output_prefix=NULL,total_coverage_file=NULL){
  coverage=read.table(cov,header=F,sep='\t')
  if (!is.null(reference_fasta)){
    require(Biostrings)
    referenceGenome = readDNAStringSet(reference_fasta, format="fasta")
    xlim_coord=c(1, length(referenceGenome[[1]])) 
    }else{
    #Just default to the human coordinates if nothing is specified
    xlim_coord=c(1,16569)
  }
  require(ggplot2)
  gg=ggplot()
  cov_data=data.frame(x=coverage[,2],y=(coverage[,4]))
  gg = gg + geom_area(data=cov_data,aes(x=x,y=y),colour='blue', stat="identity")
  gg = gg + coord_cartesian(xlim=xlim_coord)
  if(!is.null(sample_name)){
    if(!is.null(total_coverage_file)){
        #Covered by atleast 10 reads, What percentage are coverered by atleast
        #10 reads.
        reads_greater_than_10=coverage[coverage[,2] >=10,]
        read_gr_10 = round((length(reads_greater_than_10[,1])/xlim_coord[2])*100,digits=2)
        print(read_gr_10)
        cat(paste(sample_name,round((length(coverage[,1])/xlim_coord[2])* 100,digits=2), read_gr_10 ,sep='\t'),'\n',file=total_coverage_file,append=T)
    }
    gg = gg+ggtitle(paste("Coverage plot for",sample_name, "\n Cov = ", round((length(coverage[,1])/xlim_coord[2])* 100,digits=1),'%'))
      
  } else {
    gg = gg+ggtitle(paste("Coverage plot for",cov,"\n Cov = ",round((length(coverage[,1])/xlim_coord[2])* 100,digits=1),'%'))
  }
  gg = gg + xlab("Position (Bp)") + ylab("Coverage")
  if(is.null(output_prefix)){
    output_prefix=cov
  }
  print(output_prefix)
  png(filename=paste0(basename(output_prefix),'.png'))
  plot(gg)
  dev.off()
}

