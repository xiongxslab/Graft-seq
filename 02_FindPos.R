suppressMessages(library(Rsamtools))
suppressMessages(library(ShortRead))
suppressMessages(library(stringr))
suppressMessages(library(BSgenome))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(GenomicRanges))
suppressMessages(library(data.table))
suppressMessages(library(doParallel))
suppressMessages(library(dplyr))
######################################################
# 1.identify the position of the first base of reads #
######################################################
anl_rm <- function(){
  dir <- getwd()
  # bam file
  bamlist <- list.files(dir, pattern = "\\.sort.bam$")
  # screen bamlist ---------
  for(bn in bamlist){ # Process all bam files in the current path
    print(bn)
    b1 <- scanBam(file.path(dir,bn))
    b1 <- data.frame(b1[[1]], stringsAsFactors = F) 
    b1$rname <- as.character(b1$rname)
    if(length(which(b1$flag==4)) > 0){ # remove non-map
      tb1 <- b1[-which(b1$flag==4),]
    }else{
      tb1 <- b1
    }
    tb1[which(tb1$strand=='-'),5] <- tb1[which(tb1$strand=='-'),5] + tb1[which(tb1$strand=='-'),6] -1
    
    x <- paste(tb1$rname,tb1$pos,tb1$strand,sep=":")
    count <- as.data.frame(table(x))
    count <- cbind(as.data.frame(str_split_fixed(count$x, ':', 3)),count$Freq)
    colnames(count)<-c("chr","pos","strand","freq")
    count$strand <- ifelse(count$strand == '+','-','+')
    write.table(count[order(count$freq,decreasing = T),],paste0(unlist(strsplit(bn,'[.]'))[1],'.pos.coverage.txt'),sep="\t",col.names = F,row.names = F,quote=F)
  }
}
anl_rm()



##############################################
#### 2.identify the position of nearest A ####
##############################################
search_m6A <- function(x){
  # substract base(upstream 5nt,downstream 10nt)
  if(ori.info.df$down.pos[x] < length(getSeq(hg38,ori.info.df$chr[x]))){
    # Un-exceeded the chromosome length regions
    seq_info <- getSeq(hg38, ori.info.df$chr[x], ori.info.df$up.pos[x], 
                       ori.info.df$down.pos[x], strand=ori.info.df$strand[x])
  } else{
    # Exceeded the chromosome length regions
    seq_info <- getSeq(hg38, ori.info.df$chr[x], ori.info.df$up.pos[x], 
                       length(getSeq(hg38,ori.info.df$chr[x])), strand=ori.info.df$strand[x]) 
  }
  # Current position base
  base <- as.character(getSeq(hg38, ori.info.df$chr[x], ori.info.df$pos[x], 
                              ori.info.df$pos[x], strand=ori.info.df$strand[x]))
 
  # Determine whether there is A base in the the extracted regions
  if(unlist(gregexpr("A", seq_info))[1] == '-1'){
    min.distance = NA
    second.distance = NA
  } else {
    # Prioritize A at the 3' end within 5 nt
    a.distance = gregexpr("A", seq_info)[[1]] - 6
    a3.rank = which(a.distance>0 & a.distance <6)
    if(length(a3.rank)>0){
      rest.rank <- setdiff(order(abs(gregexpr("A", seq_info)[[1]] - 6)),a3.rank)
      rank <- c(a3.rank,rest.rank)
    } else{
      rank<-order(abs(gregexpr("A", seq_info)[[1]] - 6))
    }
    min.distance <- (gregexpr("A", seq_info)[[1]] - 6)[rank[1]]
    second.distance <- (gregexpr("A", seq_info)[[1]] - 6)[rank[2]]
  }
  result<-c(base,as.character(seq_info),min.distance,second.distance)
}
cl <- makeCluster(10) # Adjust the number of threads according to the situation
registerDoParallel(cl)
files<-list.files(pattern = '.pos.coverage.txt$')
for(file in files){
  print(file)
  ori.info.df<-fread(file)
  colnames(ori.info.df)<-c("chr","pos","strand","count")
  # filter chromosome
  common_chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                          "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                          "chr20", "chr21", "chr22", "chrM", "chrX", "chrY")
  ori.info.df <- subset(ori.info.df, chr %in% common_chromosomes)
  
  hg38 <- BSgenome.Hsapiens.UCSC.hg38
  ori.info.df<-ori.info.df %>% filter(count > 50)
  
  # substract region
  ori.info.df$up.pos='NA'
  ori.info.df$down.pos='NA'
  ori.info.df$up.pos[which(ori.info.df$strand=='+')]<-ori.info.df$pos[which(ori.info.df$strand=='+')]-5
  ori.info.df$down.pos[which(ori.info.df$strand=='+')]<-ori.info.df$pos[which(ori.info.df$strand=='+')]+10
  ori.info.df$up.pos[which(ori.info.df$strand=='-')]<-ori.info.df$pos[which(ori.info.df$strand=='-')]-10
  ori.info.df$down.pos[which(ori.info.df$strand=='-')]<-ori.info.df$pos[which(ori.info.df$strand=='-')]+5
  ori.info.df$up.pos<-as.numeric(ori.info.df$up.pos)
  ori.info.df$up.pos[which(ori.info.df$up.pos<0)] <- 1
  ori.info.df$down.pos<-as.numeric(ori.info.df$down.pos)
  
  info <- foreach(x=1:nrow(ori.info.df),.combine = 'rbind',.packages="BSgenome") %dopar% search_m6A(x)
  # info <- foreach(x=1000:1500,.combine = 'rbind',.packages="BSgenome") %dopar% search_m6A(x)
  result<-cbind(ori.info.df,info)
  colnames(result)<-c("chr","pos","strand","count","up.pos","down.pos","base","-5+10nt.seq","min.distance","second.distance")
  # result$min.pos<-result$pos + as.numeric(result$min.distance)
  result$min.pos<-ifelse(result$strand=='-',result$pos-as.numeric(result$min.distance),result$pos+as.numeric(result$min.distance))
  
  # filter no A base whithin regions
  result <- result[-which(is.na(result$min.distance)),]
  write.table(result,paste0('./1.3nearest.A/',unlist(strsplit(file,'[.]'))[1],'.pos.coverage.m6A.txt'),
              sep='\t',quote=F,row.names=F,col.names=T)
}
stopCluster(cl)
