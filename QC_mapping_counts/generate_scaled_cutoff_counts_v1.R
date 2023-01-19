#reads in raw R species count files and, using run group key, perform run to run scaling 
# a new set of scaled count files are created (can be used for phyloseq analyses)
# the new scaled count files are then subject to species two species filter steps (see below)

#Requires: ITS_run_group_key_batches_passing.csv - NB this contains only the samples which
#     will be processed ie. batches passing initial batch QC

library("stringr")
library("dplyr")
library("ggplot2")

dir.create("../raw_counts_scaled")
#### create run to run scaled counts ####
#read in raw data
full_file_list=list.files(pattern=".Rspecies_counts.csv")
file_list=full_file_list[file.size(full_file_list) > 0]

#combine into table
tbl <- lapply(file_list, function(x){
  read.table(x, sep=",", as.is=TRUE, header=FALSE)
})
new_tbl <- do.call(rbind, tbl)
colnames(new_tbl) <- c("Counts","Species","Sample")

#convert read counts to read pair counts 
#halve the counts (so that 1 = 1 read pair)
new_tbl$Counts <- new_tbl$Counts/2
#assign samples to run group 
ITS_group_key=as.data.frame(read.csv("../../Key_files/ITS_run_group_key_batches_passing.csv", header=F, sep=","))
colnames(ITS_group_key) <- c("Sample","Run_group")

#collate run group key with data
new_tbl <- merge(ITS_group_key, new_tbl, by = "Sample",all=T)

#converts NaNs  counts to 0 
new_tbl$Counts[is.na(new_tbl$Counts)] <- 0
#converts NaNs  counts to 0 
new_tbl$Species[is.na(new_tbl$Species)] <- "Empty"

### prepare table for run to run scaling 
## Calculate average sample counts scaling factor - average sample counts (over full dataset) ###
#get total counts of dataset
dataset_counts <- new_tbl %>% 
  dplyr::mutate(total_counts = sum(Counts)) %>% select(total_counts) %>% unique()

#get number of samples in dataset
dataset_samples <- new_tbl %>% select(Sample,Run_group) %>% unique() %>% tally()

#get average sample counts (for scaling all groups to)
av_counts <- round(mean(dataset_counts$total_counts/dataset_samples$n))


## use av acounts to define new total counts per group 
#get number of samples per group
samples_per_group <- new_tbl %>% select(Sample,Run_group) %>% unique() %>%
  dplyr::group_by(Run_group) %>% tally()

samples_per_group$scaled.group.counts <- samples_per_group$n * av_counts

#merge newly scaled group counts with all data 
#create merge of sample counts, group counts and and run group data
group_merge.tmp <- merge(samples_per_group,new_tbl, by = "Run_group", all=T)

#calculate existing total group counts
group_data.tmp2 <- group_merge.tmp %>% 
  dplyr::group_by(Run_group) %>% 
  dplyr::mutate(group.counts = sum(Counts)) 

#calculate total counts per sample
group_data <- group_data.tmp2 %>% 
  dplyr::group_by(Sample) %>% 
  dplyr::mutate(sample_counts = sum(Counts))

### Perform run to run scaling 
##scale sample counts - ratio of each sample relative to total group counts, multiplied by newly scaled group counts 
#scaled.group.counts * (sample_counts/group.counts)
for (i in seq_along(group_data$Sample)) {
  print(group_data$Sample[i])
  group_data$scaled_sample_counts[i] <- round(group_data$scaled.group.counts[i] * 
                                                     (group_data$sample_counts[i]/group_data$group.counts[i]))
}

#assess changes in sample counts pre and post run-run scaling
tiff(paste("pre_scaling.tiff",sep=""),width = 1600, height = 800)
print(
ggplot(group_data, aes(x=Sample,y=sample_counts))  + 
  geom_bar(stat="identity", width =0.8) + facet_wrap(~ Run_group, scales="free_x")  +
  theme(axis.text.x = element_text(color = "grey20", size = 14, angle =90, hjust = .5, vjust = .5, face = "plain")) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1)))
)
dev.off()                                                                                                                                              

tiff(paste("post_scaling.tiff",sep=""),width = 1600, height = 800)
print(
  ggplot(group_data, aes(x=Sample,y=scaled_sample_counts))  +
    geom_bar(stat="identity", width =0.8) + facet_wrap(~ Run_group, scales="free_x") +
    theme(axis.text.x = element_text(color = "grey20", size = 14, angle =90, hjust = .5, vjust = .5, face = "plain")) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1)))
)
dev.off()  

#converts NaNs to 0 
group_data$scaled_sample_counts[is.na(group_data$scaled_sample_counts)] <- 0

#adjust species counts using scaled sample counts - multiplied against ratio of raw species counts to total sample counts
#scaled_sample_counts * (Counts/sample_counts)
group_data$scaled_counts<- round(group_data$scaled_sample_counts*(group_data$Counts/group_data$sample_counts))

#converts NaNs to 0 
group_data$scaled_counts[is.na(group_data$scaled_counts)] <- 0

### create new csv files of run-run scaled data 
#select only counts, species & samples
new_counts.tmp <- group_data %>% select(scaled_counts, Species, Sample)

#remove samples with Empty data 
new_counts <- new_counts.tmp %>% filter(Species != "Empty")
new_counts <- subset(new_counts.tmp,new_counts.tmp$Species != "Empty")
# Split dataframe by sample
split_df <- split(new_counts, list(new_counts$Sample))

# Write out separate csv for each sample
for (i in names(split_df)) {
  write.table(split_df[[i]], paste0("../raw_counts_scaled/",i, ".scaled_counts.csv"),row.names = FALSE, col.names = FALSE,sep=",")
}

#########################################################################
setwd("../raw_counts_scaled/")
####Read in files created above and create filtered cutoff counts #####
file_list2=(list.files(pattern=".scaled_counts.csv"))
   
## filter species based on two filters:
# species read cutoff (SRC) - eg. at least 10 reads for species to remain
# species percent cutoff (SPC) - eg. at least 0.2% of total species to remain
for (i in seq_along(file_list2)) {
  filename =file_list2[[i]]
  samplename=str_replace(filename,".csv","")
  sample=str_replace(filename,".scaled_counts.csv","")
  sample_only=sample %>% str_replace("ITS[0-9]*_","") %>% str_replace("no_Hs_","")
  species_read_cutoff = 10
  val=0.002
  percent_cutoff=val*100
  percent_cutoff_nopunc=sub("[.]","_",percent_cutoff)
  if (file.size(file_list2[[i]]) > 0) {
    #read counts into R
    counts=as.data.frame(read.csv(filename, header=FALSE, sep=",", col.names=c("Counts","Species","Sample")))
    #remove species if counts are less than cutoff
    species_read_counts=counts %>% filter(Counts >= species_read_cutoff)
    #check if any data remains
    if (sum(species_read_counts$Counts) > 0) {
      #add a Proportional_counts column
      species_read_counts[4]=species_read_counts["Counts"]/sum(species_read_counts["Counts"])
      colnames(species_read_counts)[4] = "Proportional_counts"
      #filter counts by cutoff
      cutoff_filter_counts=species_read_counts %>% filter(Proportional_counts > val)
      #check if any data remains (again)
      if (sum(cutoff_filter_counts$Counts) > 0) {
        cutoff_filter_counts[4]=cutoff_filter_counts["Counts"]/sum(cutoff_filter_counts["Counts"])
        write.csv(cutoff_filter_counts, paste(samplename,"_SPC_",percent_cutoff_nopunc,"_SRC_", species_read_cutoff,".csv", sep=""), row.names=F)#create a file of cutoff filtered counts
      } else {
        print(paste(sample,"did not pass second filter stage",sep=" "))
        Counts <- 0
        Species <- "ZEmpty"
        Sample <- sample_only
        Proportional_counts <- 1
        test <-data.frame(Counts,Species,Sample,Proportional_counts)
        write.csv(test,paste(samplename,"_SPC_",percent_cutoff_nopunc,"_SRC_", species_read_cutoff,".csv", sep=""),row.names=F)
        
      }
    } else {
      print(paste(sample,"did not pass first filter stage",sep=" "))
      Counts <- 0
      Species <- "ZEmpty"
      Sample <- sample_only
      Proportional_counts <- 1
      test <-data.frame(Counts,Species,Sample,Proportional_counts)
      write.csv(test,paste(samplename,"_SPC_",percent_cutoff_nopunc,"_SRC_", species_read_cutoff,".csv", sep=""),row.names=F)
    }
  }  else {
    print(paste(file_list2[[i]],"has no data"))
    #create empty final counts file (for plotting where empty)
    Counts <- 0
    Species <- "ZEmpty"
    Sample <- sample_only
    Proportional_counts <- 1
    test <-data.frame(Counts,Species,Sample,Proportional_counts)
    write.csv(test,paste(samplename,"_SPC_",percent_cutoff_nopunc,"_SRC_", species_read_cutoff,".csv", sep=""),row.names=F)

  }
  
  
}
# Output: species filtered count files named as eg. NTC03.scaled_counts_scaled_SPC_0_2_SRC_10.csv
#have headers, and an additional 'proportional counts' column
