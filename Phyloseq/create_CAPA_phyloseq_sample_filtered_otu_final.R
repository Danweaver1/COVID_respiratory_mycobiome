#### Takes in individual count files, filters data,
##      and creates an otu matrix for use by Phyloseq
#uses raw R species counts.csv files as input (counts in long format: <Counts>,<Species>,<Sample>)
#NB raw files do not contain header
#converts to matrix (samples filtered by total read counts cutoff, set by 'val')
#filters samples based on highest total sample counts NEC per batch
#     samples with counts less than the highest corresponding NEC are removed
#converts passing sample data to otu table for phyloseq
#prints otu table as csv in upper directory

#load libraries
library("stringr")
library("dplyr")
library(tidyr)

## Requires: 
# count data as '*counts.csv' files described above
# key of extraction groups passing QC 'Extraction_group_key_batches_passing.csv'
#         (stored in directory above)
# 

full_file_list=(list.files(pattern="counts.csv"))
#filter samples by total read counts
file_list=full_file_list[file.size(full_file_list) > 0]

files_to_keep <- list()
for (i in seq_along(file_list)) {
  filename =file_list[[i]]
  samplename=str_replace(filename,".csv","")
  #read counts into R
  #no header
  counts=as.data.frame(read.csv(filename, header=F, sep=","))
  #add header
  colnames(counts) <- c("Counts","Species","Sample")
  val=500
  if (sum(counts$Counts) > val) {
    files_to_keep <- append(files_to_keep, list(filename))
  }
}


tbl <- lapply(files_to_keep, function(x){
  read.table(x, sep=",", as.is=TRUE, header=F)
})

tbl2 <- do.call(rbind, tbl)

#add header
colnames(tbl2) <- c("Counts","Species","Sample")

#calculate total counts per sample
counts_per_sample <- tbl2 %>% dplyr::group_by(Sample) %>% 
  dplyr::mutate(total_sample_counts = sum(Counts)) %>% select(total_sample_counts,Sample) %>% unique()

#assign samples to extraction batch
Xn_group_key=as.data.frame(read.csv("../Extraction_group_key_batches_passing.csv", header=F, sep=","))
colnames(Xn_group_key) <- c("Sample","Extraction_group")

#collate extraction group key with data
ps_groups <- merge(Xn_group_key, counts_per_sample, by = "Sample",all=T)
ps_groups$total_sample_counts[is.na(ps_groups$total_sample_counts)] <- 0
#subset to only NECs
ps_groups_NECs <-  ps_groups[(grep("^NEC",ps_groups$Sample)),]

#get highest NEC counts per extraction batch  
count_cutoffs <- ps_groups_NECs %>% dplyr::group_by(Extraction_group) %>%
  dplyr::mutate(count_cutoff = max(total_sample_counts)) %>% select(Extraction_group,count_cutoff) %>% unique()
count_cutoffs$Extraction_group[is.na(count_cutoffs$Extraction_group)] <- 0
#subset to only samples
ps_groups_samples.tmp <-  ps_groups[-(grep("^NEC",ps_groups$Sample)),]
ps_groups_samples <-  ps_groups_samples.tmp[-(grep("^NTC",ps_groups_samples.tmp$Sample)),]
ps_groups_samples$Extraction_group[is.na(ps_groups_samples$Extraction_group)] <- 0

#subset samples using their corresponding NEC count cutoff
ps_groups_samples$contam <- NA

for (i in seq_along(unique(ps_groups_samples$Sample))) {
  samp <- ps_groups_samples$Sample
  print(samp[i])
  group <- ps_groups_samples$Extraction_group
  print(group[i])
  cutoff <- count_cutoffs$count_cutoff[which(count_cutoffs$Extraction_group==group[i])]
  print(cutoff)
  if (ps_groups_samples$total_sample_counts[i] > cutoff) {
    ps_groups_samples$contam[i] <- "KEEP"
    print(ps_groups_samples$total_sample_counts[i])
    print(cutoff)
  } else {
    ps_groups_samples$contam[i] <- "REMOVE"
    print(ps_groups_samples$total_sample_counts[i])
    print(cutoff)
  }
}

#reduce info to only those passing NEC per batch filter
NEC_pb_filtered_samples <-  ps_groups_samples[(-grep("REMOVE",ps_groups_samples$contam)),] 
#write filtered file set to file
write.table(NEC_pb_filtered_samples[,-c(3:4)],file="sp_NEC_filtered_sample_list.csv",sep=",",quote=F,col.names = F,row.names = F)
#subset full count data to filtered samples
new_tbl <- subset(tbl2, Sample %in% NEC_pb_filtered_samples$Sample )

write.table(new_tbl,file="full_counts_sp_NEC_filtered.csv",sep=",",quote=F,col.names = T,row.names = F)



#convert long table to wide format (matrix) - species as rows

mat <- new_tbl %>% tidyr::spread(Sample, Counts)
#convert NA to zeros
mat[is.na(mat)] <- 0

#remove dummy database row 
#mat <- mat[,-which(colnames(mat)=="database")]


write.csv(mat,"../out",row.names = F)
