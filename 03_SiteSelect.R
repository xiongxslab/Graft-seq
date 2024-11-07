suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
#######################
# 1.replicated filter #
#######################
setwd("/data/") # The folder needs to have 4 files corresponding to rep1/rep2 of Exp/Ctrl
files=list.files(pattern = '*.pos.coverage.m6A.txt')
for(i in 1:2){
  #### rep1 ####
  con1 <- fread(files[4*i-3]) %>% 
    select(chr,min.pos,strand,count) %>%
    mutate(id = paste(chr,min.pos,strand,sep='|'))
  print(paste("total con1 m6A sites: ",nrow(con1)))
  print(paste("uniq con1 m6A sites: ",length(unique(con1$id))))
  con1 <- con1 %>% group_by(id) %>% summarise(count = sum(count)) # For repeated positions, add their counts
  con1 <- cbind(data.frame(do.call('rbind', strsplit(as.character(con1$id), '\\|'))),con1$count,con1$id)
  colnames(con1) <- c("chr", "pos", "strand", "count", "id")
  
  jump1<-fread(files[4*i-1]) %>% 
    select(chr,min.pos,strand,count) %>%
    mutate(id = paste(chr,min.pos,strand,sep='|'))
  print(paste("total jump1 m6A sites: ",nrow(jump1)))
  print(paste("uniq jump1 m6A sites: ",length(unique(jump1$id))))
  jump1 <- jump1 %>% group_by(id) %>% summarise(count = sum(count))
  jump1 <- cbind(data.frame(do.call('rbind', strsplit(as.character(jump1$id), '\\|'))),jump1$count,jump1$id)
  colnames(jump1) <- c("chr", "pos", "strand", "count", "id")
  
  # calculate value (if Exp yes，Ctrl no，value=Exp count; if both yes, value=Exp/Ctrl)
  df_merged<-left_join(jump1,con1,by='id',suffix=c('.jump','.con')) %>%
    mutate(value = case_when(is.na(count.con) ~ count.jump,
                             TRUE ~ count.jump/count.con)) %>%
    select(-id,-chr.con,-pos.con,-strand.con,-count.con)
  colnames(df_merged) <- c("chr", "pos", "strand", "count", "value")
  df_merged <- df_merged[order(df_merged$value,decreasing = T),]
  
  print(paste("uniq rep1 m6A sites: ",nrow(df_merged)))
  name=paste0(unlist(strsplit(files[4*i-3],'-'))[2],'-',unlist(strsplit(files[4*i-3],'-'))[3],"-rep1.m6A.filter.txt")
  write.table(df_merged,file=name,sep="\t",
              row.names = F,col.names = T,quote=F)
  
  #### rep2 ####
  con2 <- fread(files[4*i-2]) %>% 
    select(chr,min.pos,strand,count) %>%
    mutate(id = paste(chr,min.pos,strand,sep='|'))
  print(paste("total con2 m6A sites: ",nrow(con2)))
  print(paste("uniq con2 m6A sites: ",length(unique(con2$id))))
  con2 <- con2 %>% group_by(id) %>% summarise(count = sum(count))
  con2 <- cbind(data.frame(do.call('rbind', strsplit(as.character(con2$id), '\\|'))),con2$count,con2$id)
  colnames(con2) <- c("chr", "pos", "strand", "count", "id")
  
  jump2<-fread(files[4*i]) %>% 
    select(chr,min.pos,strand,count) %>%
    mutate(id = paste(chr,min.pos,strand,sep='|'))
  print(paste("total jump2 m6A sites: ",nrow(jump2)))
  print(paste("uniq jump2 m6A sites: ",length(unique(jump2$id))))
  jump2 <- jump2 %>% group_by(id) %>% summarise(count = sum(count))
  jump2 <- cbind(data.frame(do.call('rbind', strsplit(as.character(jump2$id), '\\|'))),jump2$count,jump2$id)
  colnames(jump2) <- c("chr", "pos", "strand", "count", "id")
  
  df_merged<-left_join(jump2,con2,by='id',suffix=c('.jump','.con')) %>%
    mutate(value = case_when(is.na(count.con) ~ count.jump,   # a表有，b表没有的行，value = jump count
                             TRUE ~ count.jump/count.con)) %>%
    select(-id,-chr.con,-pos.con,-strand.con,-count.con)
  colnames(df_merged) <- c("chr", "pos", "strand", "count", "value")
  df_merged <- df_merged[order(df_merged$value,decreasing = T),]
  
  print(paste("uniq rep2 m6A sites: ",nrow(df_merged)))
  name=paste0(unlist(strsplit(files[4*i-3],'-'))[2],'-',unlist(strsplit(files[4*i-3],'-'))[3],"-rep2.m6A.filter.txt")
  write.table(df_merged,file=name,sep="\t",
              row.names = F,col.names = T,quote=F)
}

#######################
# 2.two rep intersect #
#######################
files=list.files(pattern = '*.m6A.filter.txt')
# cutoff
c1 = 0 # count_rep1
v1 = 0 # value_rep1
c2 = 0 # count_rep2
v2 = 0 # value_rep2
for(i in 1:1){
  rep1<-fread(files[2*i-1]) %>% mutate(label = paste(chr,pos,strand,sep='|'))
  rep2<-fread(files[2*i]) %>% mutate(label = paste(chr,pos,strand,sep='|')) %>% select(label,count,value)
  full.rep <- full_join(rep1,rep2,by="label",suffix = c("_rep1", "_rep2")) %>%
    select(chr,pos,strand,count_rep1,value_rep1,label,count_rep2,value_rep2)
  name1=paste0(unlist(strsplit(files[2*i-1],'-'))[1],'-',unlist(strsplit(files[2*i-1],'-'))[2],".full.m6A.txt")
  print(paste0('full m6A sites: ',nrow(full.rep)))
  write.table(full.rep,file=name1,sep="\t",
              row.names = F,col.names = T,quote=F)
  
  inter.rep<-inner_join(rep1,rep2,by="label",suffix = c("_rep1", "_rep2")) %>% 
    select(chr,pos,strand,count_rep1,value_rep1,count_rep2,value_rep2) %>%
    filter(count_rep1>c1,value_rep1>v1,count_rep2>c2,value_rep2>v2)
  name=paste0(unlist(strsplit(files[2*i-1],'-'))[1],'-',unlist(strsplit(files[2*i-1],'-'))[2],".m6A.txt")
  print(paste0('m6A sites: ',nrow(inter.rep)))
  write.table(inter.rep,file=name,sep="\t",
              row.names = F,col.names = T,quote=F)
}
