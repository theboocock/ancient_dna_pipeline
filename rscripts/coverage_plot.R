#http://learnr.wordpress.com/2009/04/29/ggplot2-labelling-data-series-and-adding-a-data-table/


coverage_plot=function(cov,reference_fasta=NULL, sample_name=NULL,output_prefix=NULL,total_coverage_file=NULL,contamination_mapping = NULL){
  coverage=read.table(cov,header=F,sep='\t')
  if (!is.null(reference_fasta)){
    require(Biostrings)
    referenceGenome = readDNAStringSet(reference_fasta, format="fasta")
    if(!is.null(contamination_mapping)){
        one_to_use =  grep(contamination_mapping,names(referenceGenome),fixed=T)
        print(one_to_use)
        xlim_coord = c(1,length(referenceGenome[[one_to_use]]))
    }else{
    # Must just be the first one because no mapping was specified
    xlim_coord=c(1, length(referenceGenome[[1]])) 
    }
    }else{
    #Just default to the human coordinates if nothing is specified
    xlim_coord=c(1,16569)
  }
  require(ggplot2)
  coverage_new = data.frame(pos=c(1:xlim_coord[2]),cov=rep(0,xlim_coord[2]))
  for (i in 1:nrow(coverage)){
      coverage_new[coverage[i,2],2] = coverage[i,4] 
  }
  gg=ggplot()
  cov_data=data.frame(x=coverage_new[,1],y=(coverage_new[,2]))
  gg = gg + geom_area(data=cov_data,aes(x=x,y=y),colour='blue', stat="identity")
  gg = gg + coord_cartesian(xlim=xlim_coord)
  if(!is.null(sample_name)){
    if(!is.null(total_coverage_file)){
        #Covered by atleast 10 reads, What percentage are coverered by atleast
        #10 reads.
        reads_greater_than_10=coverage[coverage[,4] >=10,]
        read_gr_10 = round((length(reads_greater_than_10[,1])/xlim_coord[2])*100,digits=2)
        mean_dat = sum(coverage[,4])/xlim_coord[2]
        diff = xlim_coord[2] - length(coverage[,2])
        add_zeros = rep(0,diff)
        dat = c(coverage[,4],add_zeros)
        std_dev = (dat - mean_dat)^2
        std_dev = sqrt(sum(std_dev)/(xlim_coord[2]-1)) 
        cat(paste(sample_name,round((length(coverage[,1])/xlim_coord[2])* 100,digits=2), mean_dat, std_dev ,sep='\t'),'\n',file=total_coverage_file,append=T)
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
  print(length(coverage[,1]))
  print(xlim_coord)
  png(filename=paste0(basename(output_prefix),'.png'))
  plot(gg)
  dev.off()
}


