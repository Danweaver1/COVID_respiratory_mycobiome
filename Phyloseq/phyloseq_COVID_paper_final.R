#### COVID ITS mycobiome study - phyloseq data analysis 

## NB for diversity analyses to work correctly, counts should correspond to read pairs. 
#     Some diversity metrics take singletons into account, if counts represent each read of a pair,
#       there will falsely appear to be no singletons in data.
# ie. If counts correspond to each read, halve the otu data prior to creating phyloseq object

## File requirements:
# count matrix csv file 'out' 
# taxa tsv file 'UNITE_taxa_sp_v2.tsv'
# species level colour table 'colour_table_v5.csv'
# genus level colour table 'genus_colour_table_v2.csv'

### Initial setup ####
#set working dir
setwd("")

#load packages
library("devtools")
library("phyloseq"); packageVersion("phyloseq")
library("plyr")
library("ggplot2"); packageVersion("ggplot2")
library("stringr")
library("microbiomeSeq")
library("readxl")

#set theme
theme_set(theme_bw())
theme_set(theme_light())

##### manually set colours - species #####
colour_table <- read.csv("colour_table_v5.csv",header=TRUE)
colour_table[,1:3] <- lapply(colour_table[,1:3] , as.character)

##### manually set colours - genera #####
genus_colour_table <- read.csv("genus_colour_table_v2.csv",header=TRUE)
genus_colour_table[,1:2]  <- lapply(genus_colour_table[,1:2] , as.character)

### Read files in ####
otu_mat <- as.matrix(read.table("out", sep =",", header=TRUE))
tax_mat <- as.matrix(read.table("UNITE_taxa_sp_v2.tsv", sep ="\t", header=TRUE))

#read in clinical data - only the discrete variables able to compare 
all_samples <- read_excel("clinical_data/COVID_clinical_data_comparison_variables_only_V2.xlsx", sheet = "samples")
all_samples[,1:29]  <- lapply(all_samples[,1:29] , as.character)

#Define row names from otu column "X"
row.names(otu_mat) <- otu_mat[,1]

#Remove the column from the matrix
otu_mat <- otu_mat[,-1]

#Save the data as integer
class(otu_mat) <- "integer"

#Assign taxa rownames for taxa matrix
row.names(tax_mat) <- tax_mat[,1]
#give informative taxa column a name 
colnames(tax_mat)[1] <- "Taxa"

#create a variable containing first section of taxa name 
#   (ie. genus in most cases, but may be a higher tax. rank for some)
#     this is to create a more informative variable indicating genus or the closest classified rank.
GO <- tax_mat[,1]  %>% str_replace("_sp$", " sp") %>% str_remove(("_[a-z]*"))
tax_mat2 <- cbind(GO,tax_mat)
colnames(tax_mat2)[1] <- "Genus_other"

#subset clinical data for samples available 
samples_df <- subset(all_samples, Sample %in% colnames(otu_mat) )
write.table(colnames(otu_mat),file="passing_samples.csv")
#Assign rownames for tax
samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample") 

##Transform into phyloseq object
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat2)
samples = sample_data(samples_df)
physeq = phyloseq(OTU, TAX, samples) 

### raw data plotting ####
dir.create("raw_plots")

#phyloseq function tax_glom merges the OTUs with the same taxonomy, summing the abundances
physeq_family <- tax_glom(physeq,taxrank="Family")

#plot family level taxa for each sample as heatmap
tiff("raw_plots/Family_abundance_heatmap_raw.tiff",width = 2000, height = 6000)
print(phyloseq::plot_heatmap(physeq_family, taxa.label = "Genus_other", method = "NMDS",
                   distance = "bray", taxa.order = "Family",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical")
)
dev.off()

##raw genus level data
physeq_genera <- tax_glom(physeq,taxrank="Genus_other")
#fix names to genus_other
phyloseq::taxa_names(physeq_genera) <- phyloseq::tax_table(physeq_genera)[, "Genus_other"]
#check the change
phyloseq::otu_table(physeq_genera)[1:5, 1:5]

#plot genera level taxa for each sample as heatmap
tiff("raw_plots/Genus_abundance_heatmap_raw.tiff",width = 1200, height = 2000)
print(phyloseq::plot_heatmap(physeq_genera, taxa.label = "Genus_other", method = "NMDS",
                             distance = "bray", taxa.order = "Family",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical")
)
dev.off()

##filter low level species from raw data - keep species occurring above 0.2% in ANY sample
total = median(sample_sums(physeq))
physeq_filt <- filter_taxa(physeq, function(x) sum(x > total*0.002) > 0, TRUE)
physeq_filt_genera <- tax_glom(physeq_filt,taxrank="Genus_other")
phyloseq::taxa_names(physeq_filt_genera) <- phyloseq::tax_table(physeq_filt_genera)[, "Genus_other"]

#### normalise counts - standardise abundances to the median sequencing depth ####
total = median(sample_sums(physeq))
standf = function(x, t=total) round(t * (x / sum(x)))
physeqn = transform_sample_counts(physeq, standf)
### normalised data plotting - NO TAXA FILTER ####
physeqn_family <- tax_glom(physeqn,taxrank="Family")

#plot genera level taxa for each sample as heatmap
tiff("Family_abundance_heatmap_norm.tiff",width = 2000, height = 1600)
print(phyloseq::plot_heatmap(physeqn_family, taxa.label = "Genus_other", method = "NMDS",
                   distance = "bray", taxa.order = "Family",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical")
)
dev.off()

#### filter low level species from normalised data ####
# again, keeping species occurring above 0.2% in ANY sample #
physeq_abund <- filter_taxa(physeqn, function(x) sum(x > total*0.002) > 0, TRUE)

#find taxa now with empty values
taxa_to_keep <- names(which(taxa_sums(physeq_abund)>0))
#remove the empty taxa 
physeq_abund <- prune_taxa(taxa_to_keep, physeq_abund)

#create genus level data from the filtered data
physeq_abund_genera <- tax_glom(physeq_abund,taxrank="Genus_other")


##explore raw and processed data
#taxa in raw data
tax_table(physeq)
taxa_names(physeq)
taxa_names(physeq_genera)
#samples & taxa in filtered data
nsamples(physeq_abund)
ntaxa(physeq_abund)
taxa_names(physeq_abund)
taxa_names(physeq_abund_genera)

#total counts per sample (un-normalised, filtered data)
colSums(otu_table(physeq_filt))
median(colSums(otu_table(physeq_filt)))
mean(colSums(otu_table(physeq_filt)))
#calculate range of sample counts (un-normalised, filtered data)
cnts <- colSums(otu_table(physeq_filt))
range(cnts)

#ntaxa (species level) per sample 
taxa_per_sample <- colSums(otu_table(physeq_abund) > 0 )
range(taxa_per_sample)
median(taxa_per_sample)
mean(taxa_per_sample)
#ntaxa (genus level) per sample 
taxa_genera_per_sample <- colSums(otu_table(physeq_abund_genera) > 0 )
ntaxa(physeq_abund_genera)
range(taxa_genera_per_sample)
median(taxa_genera_per_sample)
mean(taxa_genera_per_sample)

#### Species level data ####

#plot species level taxa for each sample as heatmap
tiff("Species_abundance_heatmap.tiff",width = 1800, height = 1800)
print(phyloseq::plot_heatmap(physeq_abund, taxa.label = "Taxa", method = "NMDS",
                   distance = "bray", taxa.order = "Genus",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 28, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 28, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical"))
dev.off()

#### Genus level data ####

#plot genera level taxa for each sample as barplot
tiff("sample_genus_barplot.tiff",width = 1400, height = 800)
print(
  plot_bar(physeq_abund_genera, fill = "Genus_other") + 
  geom_bar(aes(fill=Genus_other), stat="identity", position="stack") +
  #theme(legend.position="none") +
  scale_fill_manual(values= genus_colour_table$Colour,breaks = genus_colour_table$species) +
  ylab("Abundance (Genus level)") + #set colours using colour table, set labels using formatted species names
  theme(panel.background = element_rect(fill = "white"), #change background to white
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 28, angle =90, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 28, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=29),
          legend.text=element_text(size=27),
          legend.direction="vertical",
          legend.position ="bottom") +
  guides(fill = guide_legend(ncol = 4))
)
dev.off()

### Genus level data ####


#Transform to percentages of total available.
physeq_abund_genera.perc = transform_sample_counts(physeq_abund_genera, function(x) 100 * x/sum(x))

#remove taxa at or below 1% in all merged samples (to keep plot legend simple)
physeq_abund_genera.perc <- filter_taxa(physeq_abund_genera.perc, function(x) sum(x > 1) > 0, TRUE)

#Retransform to percentages of total available.
physeq_abund_genera.perc = transform_sample_counts(physeq_abund_genera.perc, function(x) 100 * x/sum(x))


#transform counts to proportions
physeq_abund_genera.prop = transform_sample_counts(physeq_abund_genera, function(x)  x/sum(x))
#check the change
head(otu_table(physeq_abund_genera.prop))

#plot genera level taxa for each sample as heatmap
genus_heatmap <- phyloseq::plot_heatmap(physeq_abund_genera, taxa.label = "Genus_other", method = "NMDS",
                   distance = "bray", taxa.order = "Family",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 28, angle =90, hjust=0.95, vjust=0.2, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 35, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.key.size=unit(1.5,'cm'),
              legend.title=element_text(size=35),
              legend.text=element_text(size=30),
              legend.direction="vertical")

#save as svg (NB: svg won't render properly with phyloseq's plot_heatmap unless dpi is increased)
ggsave("Genus_abundance_heatmap.svg", genus_heatmap, device = "svg", height = 32, width = 32, dpi=2400)

#### Diversity ####
# Alpha diversity on input data (not normalised or taxa filtered)
tiff("richness_raw.tiff",width = 1200, height = 400)
print(plot_richness(physeq,measures=c("Chao1", "Shannon","Observed")))
dev.off()
# Alpha diversity on normalised data (not taxa filtered)
tiff("richness_norm.tiff",width = 1200, height = 400)
print(plot_richness(physeqn,measures=c("Chao1", "Shannon","Observed")))
dev.off()
# Alpha diversity on normalised & taxa filtered data 
tiff("richness_abund.tiff",width = 1200, height = 400)
print(plot_richness(physeq_abund,measures=c("Chao1", "Shannon","Observed"))) 
dev.off()

#### Ordination ####
##OTU ordination
physeq.ord <- ordinate(physeq_abund, "NMDS", "bray")

tiff("OTU_ordination.tiff",width = 500, height = 400)
print(plot_ordination(physeq_abund, physeq.ord, type="taxa", color="Class", title="OTUs"))
dev.off()

##samples ordination
tiff("sample_ordination_CAPA_Corticos.tiff",width = 500, height = 400)
plot_ordination(physeq_abund, physeq.ord, type="samples", color="Corticos", 
                shape="CAPANoCAPA", title="Samples") + geom_point(size=3)
dev.off()
##Network
tiff("sample_network_CAPA_Corticos.tiff",width = 800, height = 400)
plot_net(physeq_abund, distance = "bray", type = "samples", 
         maxdist = 0.7, color="Corticos", point_label="CAPANoCAPA")
dev.off()

tiff("sample_network_norm_CAPA_Corticos.tiff",width = 800, height = 400)
plot_net(physeq_abund, distance = "bray", type = "samples", 
         maxdist = 0.7, color="Corticos", point_label="CAPANoCAPA")
dev.off()

tiff("raw_plots/sample_network_CAPA_Corticos_raw.tiff",width = 800, height = 400)
plot_net(physeq_abund, distance = "bray", type = "samples", 
         maxdist = 0.7, color="Corticos", point_label="CAPANoCAPA")
dev.off()

##### Clinical data comparisons ########
#### Specific comparisons ###########
dir.create("specific_comparisons")
### CAPA comparisons ####
dir.create("specific_comparisons/CAPA")
##full sample species barplots faceted by CAPA
tiff("specific_comparisons/CAPA/CAPANoCAPA_percentabundance_barplot_facet_1percent_cutoff.tiff",width = 2800, height = 1000)
print(
  plot_bar(physeq_abund_genera.perc, fill = "Genus_other") + 
    geom_bar(aes(fill=Genus_other), stat="identity", position="stack") +
    facet_grid(~CAPANoCAPA,scales="free_x", space="free")  +
    scale_fill_manual(values= genus_colour_table$Colour,breaks = genus_colour_table$species) +
    ylab("Abundance (Genus level)") + #set colours using colour table, set labels using formatted species names
    theme(panel.background = element_rect(fill = "white"), #change background to white
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 28, angle =90, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 28, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=29),
          legend.text=element_text(size=27),
          legend.direction="vertical",
          legend.position ="bottom",
          strip.text.x = element_text(size = 35, colour = "black", angle = 0)) +
    guides(fill = guide_legend(ncol = 6))
)
dev.off()

## merge data by CAPA status & plot ###
physeq_abund_CAPAnoCAPA <- merge_samples(physeq_abund, "CAPANoCAPA")

#Transform to percentages of total available.
physeq_abund_CAPAnoCAPA = transform_sample_counts(physeq_abund_CAPAnoCAPA, function(x) 100 * x/sum(x))

#create proportional bar plot of merged data
tiff("specific_comparisons/CAPA/Genus_barplot_merged_CAPAnoCAPA.tiff",width = 1200, height = 500)
print(
  plot_bar(physeq_abund_CAPAnoCAPA, fill = "Genus_other") + coord_flip() + 
    ylab("Percentage of sequences") +
    xlab("CAPA") +
    geom_bar(aes(fill=Genus_other), stat="identity", position="stack") +
    theme(panel.background = element_rect(fill = "white"), #change background to white
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=27),
          legend.text=element_text(size=25, face="italic"),
          legend.direction="vertical",
          legend.position="bottom",
          strip.text.x = element_text(size = 27)) +
    guides(fill=guide_legend(ncol=4)) + #sets legend items to X column
    scale_fill_manual(values= genus_colour_table$Colour,breaks = genus_colour_table$species)
)
dev.off()

## subset data - those with CAPA ###
CAPA_physeq <- subset_samples(physeq_abund_genera, CAPANoCAPA ==1)
#plot genera level taxa for each sample as heatmap
tiff("specific_comparisons/CAPA/CAPA_Genus_abundance_heatmap.tiff",width = 580, height = 800)
print(phyloseq::plot_heatmap(CAPA_physeq, taxa.label = "Genus_other", method = "NMDS",
                   distance = "bray", taxa.order = "Family",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical")
)
dev.off()

##subset species level data 
CAPA_physeq_sp <- subset_samples(physeq_abund, CAPANoCAPA==1)

#plot genera level taxa for each sample as heatmap
tiff("specific_comparisons/CAPA/CAPA_species_abundance_heatmap.tiff",width = 800, height = 1300)
print(phyloseq::plot_heatmap(CAPA_physeq_sp, taxa.label = "Taxa", method = "NMDS",
                   distance = "bray", taxa.order = "Family",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical")
)
dev.off()
#plot noCAPA rank abundance (average counts)
#set number of top taxa to plot
N <- 10
#plot
tiff("specific_comparisons/CAPA/CAPA_species_rank_mean_abundance.tiff",width =500, height = 1100)
par(mar=c(36,10,4,1)+.1)
barplot(sort(taxa_sums(CAPA_physeq_sp), TRUE)[1:N]/nsamples(CAPA_physeq_sp), las=2,cex.lab=2, cex.axis=3,cex.names=3,ylim = c(0,50000)) 
dev.off()

## subset data - those without CAPA ###
noCAPA_physeq <- subset_samples(physeq_abund_genera, CAPANoCAPA ==0)
#plot genera level taxa for each sample as heatmap
tiff("specific_comparisons/CAPA/noCAPA_Genus_abundance_heatmap.tiff",width = 1200, height = 800)
print(phyloseq::plot_heatmap(noCAPA_physeq, taxa.label = "Genus_other", method = "NMDS",
                   distance = "bray", taxa.order = "Family",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical")
)
dev.off()

##subset species level data 
noCAPA_physeq_sp <- subset_samples(physeq_abund, CAPANoCAPA==0)

#plot genera level taxa for each sample as heatmap
tiff("specific_comparisons/CAPA/noCAPA_species_abundance_heatmap.tiff",width = 1300, height = 1400)
print(phyloseq::plot_heatmap(noCAPA_physeq_sp, taxa.label = "Taxa", method = "NMDS",
                   distance = "bray", taxa.order = "Family",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical")
)
dev.off()
#plot noCAPA rank abundance (average counts)
#set number of top taxa to plot
N <- 10
#plot
tiff("specific_comparisons/CAPA/noCAPA_species_rank_mean_abundance.tiff",width = 1600, height = 1100)
par(mar=c(36,10,4,1)+.1)
barplot(sort(taxa_sums(noCAPA_physeq_sp), TRUE)[1:N]/nsamples(noCAPA_physeq_sp), las=2,cex.lab=2, cex.axis=3,cex.names=3,ylim = c(0,50000)) 
dev.off() 


#### plot top 10 species abundances together by CAPANoCAPA ####
top10 <- names(sort(taxa_sums(physeq_abund),TRUE)[1:10])
t10_data <- phyloseq::psmelt(physeq_abund) %>% dplyr::filter(OTU %in% top10 )
{
  tiff("specific_comparisons/CAPA/CAPANoCAPA_top_10_taxa_abundances.tiff",width = 600, height = 600)
  print(
    phyloseq::psmelt(physeq_abund) %>% #remove control symptomatic sample
      dplyr::filter(OTU %in% top10 ) %>%  dplyr::filter(!CAPANoCAPA == "NA" ) %>% #filter data to top 10 taxa
      ggplot(data = ., aes(x = Species, y = Abundance, fill=CAPANoCAPA)) +
      geom_boxplot(width=0.6) +
      ggtitle(paste("Top 10 taxa abundance",sep=" ")) +
      labs(y = "Abundance\n") +
      scale_fill_manual(values=c("#10639B","#808080")) +
      theme(panel.background = element_blank(),#change background to white
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.line = element_line(size = 1, colour = "black"),
            axis.ticks = element_line(size = 2),
            axis.text.x = element_text(color = "grey20", size = 20, angle =45, hjust = 1, vjust = 1, face = "plain"),
            axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = .5, face = "plain"),
            axis.title.x = element_text(face = "bold", color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0),
            axis.title.y = element_text(face = "bold", color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5),
            plot.title = element_text(size = 20, face = "bold",hjust=0.5, vjust=-.5),
            legend.title=element_text(size=27),
            legend.text=element_text(size=26),
            legend.direction="vertical",
            legend.position ="bottom",
            plot.margin = ggplot2::margin(1,1,1,1, "cm")) 
  )
  dev.off()
}

### SystemicAFbeforeatICUAdmin comparisons ####
dir.create("specific_comparisons/sAF")
##full sample species barplots faceted by SystemicAFbeforeatICUAdmin
tiff("specific_comparisons/sAF/SystemicAFbeforeatICUAdmin_percentabundance_barplot_facet_1percent_cutoff.tiff",width = 2800, height = 1000)
print(
  plot_bar(physeq_abund_genera.perc, fill = "Genus_other") + 
    geom_bar(aes(fill=Genus_other), stat="identity", position="stack") +
    facet_grid(~SystemicAFbeforeatICUAdmin,scales="free_x", space="free")  +
    scale_fill_manual(values= genus_colour_table$Colour,breaks = genus_colour_table$species) +
    ylab("Abundance (Genus level)") + #set colours using colour table, set labels using formatted species names
    theme(panel.background = element_rect(fill = "white"), #change background to white
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 28, angle =90, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 28, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=29),
          legend.text=element_text(size=27),
          legend.direction="vertical",
          legend.position ="bottom",
          strip.text.x = element_text(size = 35, colour = "black", angle = 0)) +
    guides(fill = guide_legend(ncol = 6))
)
dev.off()

## merge data by SystemicAFbeforeatICUAdmin status & plot ###
physeq_abund_SystemicAFbeforeatICUAdmin <- merge_samples(physeq_abund, "SystemicAFbeforeatICUAdmin")

#remove the empty taxa 
taxa_to_keep <- names(which(taxa_sums(physeq_abund_SystemicAFbeforeatICUAdmin)>0))

physeq_abund_SystemicAFbeforeatICUAdmin <- prune_taxa(taxa_to_keep, physeq_abund_SystemicAFbeforeatICUAdmin)

#collapse genera
physeq_abund_SystemicAFbeforeatICUAdmin.genera <- tax_glom(physeq_abund_SystemicAFbeforeatICUAdmin,taxrank="Genus_other")
#Transform to percentages of total available.
physeq_abund_SystemicAFbeforeatICUAdmin.genera = transform_sample_counts(physeq_abund_SystemicAFbeforeatICUAdmin.genera, function(x) 100 * x/sum(x))


#create proportional bar plot of merged data
svg("specific_comparisons/sAF/Genus_barplot_merged_SystemicAFbeforeatICUAdmin.svg",width = 16, height = 5.8)
print(
  plot_bar(physeq_abund_SystemicAFbeforeatICUAdmin.genera, fill = "Genus_other") + coord_flip() + 
    ylab("Percentage of sequences") +
    xlab("SystemicAFbeforeatICUAdmin") +
    geom_bar(aes(fill=Genus_other), stat="identity", position="stack") +
    theme(panel.background = element_blank(),#change background to white
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=27),
          legend.text=element_text(size=25, face="italic"),
          legend.direction="vertical",
          legend.position="bottom",
          strip.text.x = element_text(size = 27)) +
    guides(fill=guide_legend(ncol=4)) + #sets legend items to X column
    scale_fill_manual(values= genus_colour_table$Colour,breaks = genus_colour_table$species)
)
dev.off()

## subset data - those with SystemicAFbeforeatICUAdmin ###
SystemicAFbeforeatICUAdmin_physeq <- subset_samples(physeq_abund_genera, SystemicAFbeforeatICUAdmin ==1)
#plot genera level taxa for each sample as heatmap
tiff("specific_comparisons/sAF/SystemicAFbeforeatICUAdmin_Genus_abundance_heatmap.tiff",width = 580, height = 800)
print(phyloseq::plot_heatmap(SystemicAFbeforeatICUAdmin_physeq, taxa.label = "Genus_other", method = "NMDS",
                             distance = "bray", taxa.order = "Family",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical")
)
dev.off()

##subset species level data 
SystemicAFbeforeatICUAdmin_physeq_sp <- subset_samples(physeq_abund, SystemicAFbeforeatICUAdmin==1)
nsamples(SystemicAFbeforeatICUAdmin_physeq)
#plot genera level taxa for each sample as heatmap
tiff("specific_comparisons/sAF/SystemicAFbeforeatICUAdmin_species_abundance_heatmap.tiff",width = 800, height = 1300)
print(phyloseq::plot_heatmap(SystemicAFbeforeatICUAdmin_physeq_sp, taxa.label = "Taxa", method = "NMDS",
                             distance = "bray", taxa.order = "Family",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical")
)
dev.off()
#plot noSystemicAFbeforeatICUAdmin rank abundance (average counts)
#set number of top taxa to plot
N <- 10

#prep names for axis labels
names <- str_replace(names(sort(taxa_sums(SystemicAFbeforeatICUAdmin_physeq_sp), TRUE)[1:N]),"_", " ")

x <- as.data.frame(sort(taxa_sums(SystemicAFbeforeatICUAdmin_physeq_sp), TRUE)[1:N]/nsamples(SystemicAFbeforeatICUAdmin_physeq_sp))
colnames(x)[1] <- "Mean_abundance"
x$Species <- rownames(x)
x$Genus <- str_remove(x$Species, "_.*")
x$Species <- str_replace(x$Species,"_"," ")

##ggplot coloured barplot 
tiff("specific_comparisons/sAF/SystemicAFbeforeatICUAdmin_species_rank_mean_abundance_coloured.tiff",width = 380, height = 600)
print(
ggplot(x,aes(x=reorder(Species, -Mean_abundance),y=Mean_abundance, fill=Genus))  +
  geom_bar(stat="identity", width =0.5) + theme_bw() +
  xlab("Species") +
  scale_fill_manual(values= genus_colour_table$Colour,breaks = genus_colour_table$species) +
  scale_y_continuous(limits = c(0, 35000), labels = function(x) format(x, scientific = TRUE),
    expand = expansion(mult = c(0, 0.1))) +
  theme(panel.background = element_blank(),#change background to white
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 1, colour = "black"),
        axis.ticks = element_line(size = 2),
        axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
        axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
        axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
        plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
        legend.title=element_text(size=27),
        legend.direction="vertical",
        legend.position ="none") 
)
dev.off()
#plot
tiff("specific_comparisons/sAF/SystemicAFbeforeatICUAdmin_species_rank_mean_abundance.tiff",width = 600, height = 1100)
par(mar=c(36,10,4,1)+.1)
barplot(sort(taxa_sums(SystemicAFbeforeatICUAdmin_physeq_sp), TRUE)[1:N]/nsamples(SystemicAFbeforeatICUAdmin_physeq_sp), las=2,cex.lab=2, cex.axis=3,cex.names=3,ylim = c(0,40000),names.arg = names) 
dev.off()

## subset data - those without SystemicAFbeforeatICUAdmin ###
noSystemicAFbeforeatICUAdmin_physeq <- subset_samples(physeq_abund_genera, SystemicAFbeforeatICUAdmin ==0)
nsamples(noSystemicAFbeforeatICUAdmin_physeq)
#plot genera level taxa for each sample as heatmap
tiff("specific_comparisons/sAF/noSystemicAFbeforeatICUAdmin_Genus_abundance_heatmap.tiff",width = 700, height = 800)
print(phyloseq::plot_heatmap(noSystemicAFbeforeatICUAdmin_physeq, taxa.label = "Genus_other", method = "NMDS",
                             distance = "bray", taxa.order = "Family",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical")
)
dev.off()

##subset species level data 
noSystemicAFbeforeatICUAdmin_physeq_sp <- subset_samples(physeq_abund, SystemicAFbeforeatICUAdmin==0)

#plot genera level taxa for each sample as heatmap
tiff("specific_comparisons/sAF/noSystemicAFbeforeatICUAdmin_species_abundance_heatmap.tiff",width = 800, height = 1400)
print(phyloseq::plot_heatmap(noSystemicAFbeforeatICUAdmin_physeq_sp, taxa.label = "Taxa", method = "NMDS",
                             distance = "bray", taxa.order = "Family",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical")
)
dev.off()
#plot noSystemicAFbeforeatICUAdmin rank abundance (average counts)
#set number of top taxa to plot
N <- 10

#prep names for axis labels
names <- str_replace(names(sort(taxa_sums(noSystemicAFbeforeatICUAdmin_physeq_sp), TRUE)[1:N]),"_", " ")

no_x <- as.data.frame(sort(taxa_sums(noSystemicAFbeforeatICUAdmin_physeq_sp), TRUE)[1:N]/nsamples(noSystemicAFbeforeatICUAdmin_physeq_sp))
colnames(no_x)[1] <- "Mean_abundance"
no_x$Species <- rownames(no_x)
no_x$Genus <- str_remove(no_x$Species, "_.*")
no_x$Species <- str_replace(no_x$Species,"_"," ")

##ggplot coloured barplot 
tiff("specific_comparisons/sAF/noSystemicAFbeforeatICUAdmin_species_rank_mean_abundance_coloured.tiff",width = 380, height = 600)
print(
  ggplot(no_x,aes(x=reorder(Species, -Mean_abundance),y=Mean_abundance, fill=Genus))  +
    geom_bar(stat="identity", width =0.5) + theme_bw() +
    xlab("Species") +
    scale_fill_manual(values= genus_colour_table$Colour,breaks = genus_colour_table$species) +
    scale_y_continuous(limits = c(0, 35000), labels = function(x) format(x, scientific = TRUE),
                       expand = expansion(mult = c(0, 0.1))) +
    theme(panel.background = element_blank(),#change background to white
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
          axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=27),
          legend.direction="vertical",
          legend.position ="none") 
)
dev.off()
#plot
tiff("specific_comparisons/sAF/noSystemicAFbeforeatICUAdmin_species_rank_mean_abundance.tiff",width = 1600, height = 1100)
par(mar=c(36,10,4,1)+.1)
barplot(sort(taxa_sums(noSystemicAFbeforeatICUAdmin_physeq_sp), TRUE)[1:N]/nsamples(noSystemicAFbeforeatICUAdmin_physeq_sp), las=2,cex.lab=2, cex.axis=3,cex.names=3,ylim = c(0,50000)) 
dev.off()

##plot +/- sAF rank species data together
no_sAF_rank_sp <- as.data.frame(sort(taxa_sums(noSystemicAFbeforeatICUAdmin_physeq_sp), TRUE)[1:N]/nsamples(noSystemicAFbeforeatICUAdmin_physeq_sp))
#rename column
colnames(no_sAF_rank_sp)[1] <- "no_systemic_antifungals"
#create a species column
no_sAF_rank_sp[2] <- rownames(no_sAF_rank_sp)
colnames(no_sAF_rank_sp)[2] <- "Species"

sAF_rank_sp <- as.data.frame(sort(taxa_sums(SystemicAFbeforeatICUAdmin_physeq_sp), TRUE)[1:N]/nsamples(SystemicAFbeforeatICUAdmin_physeq_sp))
#rename column
colnames(sAF_rank_sp)[1] <- "systemic_antifungals"
#create a species column
sAF_rank_sp[2] <- rownames(sAF_rank_sp)
colnames(sAF_rank_sp)[2] <- "Species"
#merge by species column 
ranked_sp_sAF <- merge(sAF_rank_sp,no_sAF_rank_sp,by="Species", all=TRUE)
#convert to long format
library(tidyr)
ranked_sp_sAF.long <- gather(ranked_sp_sAF, sAF, Counts, 2:3, factor_key=TRUE)
#ranked_sp_sAF.long$Species <-ranked_sp_sAF.long$Species %>% str_replace("_"," ")

ranked_sp_sAF.long$Genus <- ranked_sp_sAF.long$Species %>% str_remove("_.*")
ranked_sp_sAF.long$Species <- ranked_sp_sAF.long$Species %>% str_replace("_"," ")

#sAF rank spp. facet plot
#need tidytext to order within a facet
library(tidytext)
tiff("specific_comparisons/sAF/systemicAF_rank_facet.tiff",width = 1200, height = 800)
print(
  ggplot(ranked_sp_sAF.long, aes(x=reorder_within(Species, -Counts, sAF), y=Counts, fill=Genus)) + 
    geom_bar(stat="identity", width =0.5,alpha=0.8) + theme_bw() +
    scale_x_reordered() +
    facet_wrap(~ sAF, scales="free_x") +
    xlab("Species") +
    scale_fill_manual(values= genus_colour_table$Colour,breaks = genus_colour_table$species) +
    scale_y_continuous(limits = c(0, 35000), labels = function(x) format(x, scientific = TRUE),
                       expand = expansion(mult = c(0, 0.5))) +
    theme(panel.background = element_blank(),#change background to white
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
          axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=27),
          legend.text=element_text(size=23),
          legend.direction="vertical",
          legend.position ="right") 
)
dev.off() 
##single barplot, split by sEOF
#arrange & mutate to order species by survival data only
tiff("specific_comparisons/sAF/systemicAF_rank_grouped_barplot.tiff",width = 1000, height = 600)
print(
  ranked_sp_sAF.long %>% mutate(sAF = factor(sAF, levels =c("no_systemic_antifungals","systemic_antifungals"))) %>% 
    arrange(sAF, desc(Counts)) %>%
    mutate(Species = factor(Species, levels = unique(Species))) %>%
    ggplot() + aes(x=Species,y=Counts,fill=sAF) +
    geom_col(position = "dodge", width=0.6,alpha=0.8) + 
    scale_fill_manual(values=c("#10639B","#808080")) +
    scale_y_continuous(limits = c(0, 35000), labels = function(x) format(x, scientific = TRUE),
                       expand = expansion(mult = c(0, 0.1))) +
    theme(panel.background = element_blank(),#change background to white
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
          axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=27),
          legend.text=element_text(size=23),
          legend.direction="vertical",
          legend.position ="right") 
)
dev.off()


#### plot top 10 species abundances together by SystemicAFbeforeatICUAdmin ####
top10 <- names(sort(taxa_sums(physeq_abund),TRUE)[1:10])
t10_data <- phyloseq::psmelt(physeq_abund) %>% dplyr::filter(OTU %in% top10 )
{
  tiff("specific_comparisons/sAF/SystemicAFbeforeatICUAdmin_top_10_taxa_abundances.tiff",width = 600, height = 600)
  print(
    phyloseq::psmelt(physeq_abund) %>% #remove control symptomatic sample
      dplyr::filter(OTU %in% top10 ) %>%  dplyr::filter(!SystemicAFbeforeatICUAdmin == "NA" ) %>% #filter data to top 10 taxa
      ggplot(data = ., aes(x = Species, y = Abundance, fill=SystemicAFbeforeatICUAdmin)) +
      geom_boxplot(width=0.6) +
      ggtitle(paste("Top 10 taxa abundance",sep=" ")) +
      labs(y = "Abundance\n") +
      scale_fill_manual(values=c("#10639B","#808080")) +
      theme(panel.background = element_blank(),#change background to white
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.line = element_line(size = 1, colour = "black"),
            axis.ticks = element_line(size = 2),
            axis.text.x = element_text(color = "grey20", size = 20, angle =45, hjust = 1, vjust = 1, face = "plain"),
            axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = .5, face = "plain"),
            axis.title.x = element_text(face = "bold", color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0),
            axis.title.y = element_text(face = "bold", color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5),
            plot.title = element_text(size = 20, face = "bold",hjust=0.5, vjust=-.5),
            legend.title=element_text(size=27),
            legend.text=element_text(size=26),
            legend.direction="vertical",
            legend.position ="bottom",
            plot.margin = ggplot2::margin(1,1,1,1, "cm")) 
  )
  dev.off()
}
### SurvivalEOF comparisons ####
dir.create("specific_comparisons/sEOF")

##full sample species barplots faceted by SurvivalEOF
tiff("specific_comparisons/sEOF/SurvivalEOF_percentabundance_barplot_facet_1percent_cutoff.tiff",width = 2800, height = 1000)
print(
  plot_bar(physeq_abund_genera.perc, fill = "Genus_other") + 
    geom_bar(aes(fill=Genus_other), stat="identity", position="stack") +
    facet_grid(~SurvivalEOF,scales="free_x", space="free")  +
    scale_fill_manual(values= genus_colour_table$Colour,breaks = genus_colour_table$species) +
    ylab("Abundance (Genus level)") + #set colours using colour table, set labels using formatted species names
    theme(panel.background = element_rect(fill = "white"), #change background to white
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 28, angle =90, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 28, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=29),
          legend.text=element_text(size=27),
          legend.direction="vertical",
          legend.position ="bottom",
          strip.text.x = element_text(size = 35, colour = "black", angle = 0)) +
    guides(fill = guide_legend(ncol = 6))
)
dev.off()

## merge data by SurvivalEOF status & plot ###
physeq_abund_SurvivalEOF <- merge_samples(physeq_abund, "SurvivalEOF")

#remove the empty taxa 
taxa_to_keep <- names(which(taxa_sums(physeq_abund_SurvivalEOF)>0))

physeq_abund_SurvivalEOF <- prune_taxa(taxa_to_keep, physeq_abund_SurvivalEOF)
#collapse genera
physeq_abund_SurvivalEOF.genera <- tax_glom(physeq_abund_SurvivalEOF,taxrank="Genus_other")
#Transform to percentages of total available.
physeq_abund_SurvivalEOF.genera = transform_sample_counts(physeq_abund_SurvivalEOF.genera, function(x) 100 * x/sum(x))

#create proportional bar plot of merged data
svg("specific_comparisons/sEOF/Genus_barplot_merged_SurvivalEOF.svg",width = 16, height = 7)
print(
  plot_bar(physeq_abund_SurvivalEOF.genera, fill = "Genus_other") + coord_flip() + 
    ylab("Percentage of sequences") +
    xlab("SurvivalEOF") +
    geom_bar(aes(fill=Genus_other), stat="identity", position="stack") +
    theme(panel.background = element_blank(),#change background to white
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=27),
          legend.text=element_text(size=25, face="italic"),
          legend.direction="vertical",
          legend.position="bottom",
          strip.text.x = element_text(size = 27)) +
    guides(fill=guide_legend(ncol=4)) + #sets legend items to X column
    scale_fill_manual(values= genus_colour_table$Colour,breaks = genus_colour_table$species)
)
dev.off()

## subset data - those with SurvivalEOF ###
SurvivalEOF_physeq <- subset_samples(physeq_abund_genera, SurvivalEOF ==1)
#plot genera level taxa for each sample as heatmap
tiff("specific_comparisons/sEOF/SurvivalEOF_Genus_abundance_heatmap.tiff",width = 800, height = 900)
print(phyloseq::plot_heatmap(SurvivalEOF_physeq, taxa.label = "Genus_other", method = "NMDS",
                             distance = "bray", taxa.order = "Family",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical")
)
dev.off()

##subset species level data 
SurvivalEOF_physeq_sp <- subset_samples(physeq_abund, SurvivalEOF==1)
nsamples(SurvivalEOF_physeq_sp)
#plot genera level taxa for each sample as heatmap
tiff("specific_comparisons/sEOF/SurvivalEOF_species_abundance_heatmap.tiff",width = 1000, height = 1400)
print(phyloseq::plot_heatmap(SurvivalEOF_physeq_sp, taxa.label = "Taxa", method = "NMDS",
                             distance = "bray", taxa.order = "Family",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical")
)
dev.off()
#plot noCAPA rank abundance (average counts)
#set number of top taxa to plot
N <- 10

#prep names for axis labels
names <- str_replace(names(sort(taxa_sums(SurvivalEOF_physeq_sp), TRUE)[1:N]),"_", " ")

x <- as.data.frame(sort(taxa_sums(SurvivalEOF_physeq_sp), TRUE)[1:N]/nsamples(SurvivalEOF_physeq_sp))
colnames(x)[1] <- "Mean_abundance"
x$Species <- rownames(x)
x$Genus <- str_remove(x$Species, "_.*")
x$Species <- str_replace(x$Species,"_"," ")

##ggplot coloured barplot 
tiff("specific_comparisons/sEOF/SurvivalEOF_physeq_species_rank_mean_abundance_coloured.tiff",width = 380, height = 600)
print(
  ggplot(x,aes(x=reorder(Species, -Mean_abundance),y=Mean_abundance, fill=Genus))  +
    geom_bar(stat="identity", width =0.5) + theme_bw() +
    xlab("Species") +
    scale_fill_manual(values= genus_colour_table$Colour,breaks = genus_colour_table$species) +
    scale_y_continuous(limits = c(0, 25000), labels = function(x) format(x, scientific = TRUE),
                       expand = expansion(mult = c(0, 0.1))) +
    theme(panel.background = element_blank(),#change background to white
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
          axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=27),
          legend.direction="vertical",
          legend.position ="none") 
)
dev.off()


#plot
tiff("specific_comparisons/sEOF/SurvivalEOF_species_rank_mean_abundance.tiff",width = 1600, height = 1100)
par(mar=c(36,10,4,1)+.1)
barplot(sort(taxa_sums(SurvivalEOF_physeq_sp), TRUE)[1:N]/nsamples(SurvivalEOF_physeq_sp), las=2,cex.lab=2, cex.axis=3,cex.names=3,ylim = c(0,50000)) 
dev.off()
## subset data - those without CAPA ###
noSurvivalEOF_physeq <- subset_samples(physeq_abund_genera, SurvivalEOF ==0)
#plot genera level taxa for each sample as heatmap
tiff("specific_comparisons/noSurvivalEOF_Genus_abundance_heatmap.tiff",width = 1000, height = 900)
print(phyloseq::plot_heatmap(noSurvivalEOF_physeq, taxa.label = "Genus_other", method = "NMDS",
                             distance = "bray", taxa.order = "Family",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical")
)
dev.off()

#### plot top 10 species abundances together by SurvivalEOF ####
top10 <- names(sort(taxa_sums(physeq_abund),TRUE)[1:10])
t10_data <- phyloseq::psmelt(physeq_abund) %>% dplyr::filter(OTU %in% top10 )
{
tiff("specific_comparisons/sEOF/SurvivalEOF_top_10_taxa_abundances.tiff",width = 600, height = 600)
print(
  phyloseq::psmelt(physeq_abund) %>% #remove control symptomatic sample
    dplyr::filter(OTU %in% top10 ) %>%  dplyr::filter(!SurvivalEOF == "NA" ) %>% #filter data to top 10 taxa
    ggplot(data = ., aes(x = Species, y = Abundance, fill=SurvivalEOF)) +
    geom_boxplot(width=0.6) +
    ggtitle(paste("Top 10 taxa abundance",sep=" ")) +
    labs(y = "Abundance\n") +
    scale_fill_manual(values=c("#10639B","#808080")) +
    theme(panel.background = element_blank(),#change background to white
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 20, angle =45, hjust = 1, vjust = 1, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 20, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=27),
          legend.text=element_text(size=26),
          legend.direction="vertical",
          legend.position ="bottom",
          plot.margin = ggplot2::margin(1,1,1,1, "cm")) 
)
dev.off()
}


##subset species level data 
noSurvivalEOF_physeq_sp <- subset_samples(physeq_abund, SurvivalEOF==0)
nsamples(noSurvivalEOF_physeq_sp)
#plot genera level taxa for each sample as heatmap
tiff("specific_comparisons/sEOF/noSurvivalEOF_species_abundance_heatmap.tiff",width = 1100, height = 1400)
print(phyloseq::plot_heatmap(noSurvivalEOF_physeq_sp, taxa.label = "Taxa", method = "NMDS",
                             distance = "bray", taxa.order = "Family",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical")
)
dev.off()
#plot noCAPA rank abundance (average counts)
#set number of top taxa to plot
N <- 10

#prep names for axis labels
names <- str_replace(names(sort(taxa_sums(noSurvivalEOF_physeq_sp), TRUE)[1:N]),"_", " ")

no_x <- as.data.frame(sort(taxa_sums(noSurvivalEOF_physeq_sp), TRUE)[1:N]/nsamples(noSurvivalEOF_physeq_sp))
colnames(no_x)[1] <- "Mean_abundance"
no_x$Species <- rownames(no_x)
no_x$Genus <- str_remove(no_x$Species, "_.*")
no_x$Species <- str_replace(no_x$Species,"_"," ")

##ggplot coloured barplot 
tiff("specific_comparisons/sEOF/noSurvivalEOF_physeq_species_rank_mean_abundance_coloured.tiff",width = 380, height = 600)
print(
  ggplot(no_x,aes(x=reorder(Species, -Mean_abundance),y=Mean_abundance, fill=Genus))  +
    geom_bar(stat="identity", width =0.5) + theme_bw() +
    xlab("Species") +
    scale_fill_manual(values= genus_colour_table$Colour,breaks = genus_colour_table$species) +
    scale_y_continuous(limits = c(0, 25000), labels = function(x) format(x, scientific = TRUE),
                       expand = expansion(mult = c(0, 0.1))) +
    theme(panel.background = element_blank(),#change background to white
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
          axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=27),
          legend.direction="vertical",
          legend.position ="none") 
)
dev.off()

#plot
tiff("specific_comparisons/sEOF/noSurvivalEOF_species_rank_mean_abundance.tiff",width = 1600, height = 1100)
par(mar=c(36,10,4,1)+.1)
barplot(sort(taxa_sums(noSurvivalEOF_physeq_sp), TRUE)[1:N]/nsamples(noSurvivalEOF_physeq_sp), las=2,cex.lab=2, cex.axis=3,cex.names=3,ylim = c(0,50000)) 
dev.off()

##plot +/- survival rank species data together
no_sEOF_rank_sp <- as.data.frame(sort(taxa_sums(noSurvivalEOF_physeq_sp), TRUE)[1:N]/nsamples(noSurvivalEOF_physeq_sp))
#rename column
colnames(no_sEOF_rank_sp)[1] <- "no_survival"
#create a species column
no_sEOF_rank_sp[2] <- rownames(no_sEOF_rank_sp)
colnames(no_sEOF_rank_sp)[2] <- "Species"

sEOF_rank_sp <- as.data.frame(sort(taxa_sums(SurvivalEOF_physeq_sp), TRUE)[1:N]/nsamples(SurvivalEOF_physeq_sp))
#rename column
colnames(sEOF_rank_sp)[1] <- "survival"
#create a species column
sEOF_rank_sp[2] <- rownames(sEOF_rank_sp)
colnames(sEOF_rank_sp)[2] <- "Species"
#merge by species column 
ranked_sp_sEOF <- merge(sEOF_rank_sp,no_sEOF_rank_sp,by="Species", all=TRUE)
#convert to long format
library(tidyr)
ranked_sp_sEOF.long <- gather(ranked_sp_sEOF, sEOF, Counts, 2:3, factor_key=TRUE)
#ranked_sp_sEOF.long$Species <-ranked_sp_sEOF.long$Species %>% str_replace("_"," ")

ranked_sp_sEOF.long$Genus <- ranked_sp_sEOF.long$Species %>% str_remove("_.*")
ranked_sp_sEOF.long$Species <- ranked_sp_sEOF.long$Species %>% str_replace("_"," ")

#survival EOF rank spp. facet plot
#need tidytext to order within a facet
library(tidytext)
tiff("specific_comparisons/sEOF/survivalEOF_rank_facet.tiff",width = 1200, height = 800)
print(
ggplot(ranked_sp_sEOF.long, aes(x=reorder_within(Species, -Counts, sEOF), y=Counts, fill=Genus)) + 
  geom_bar(stat="identity", width =0.5,alpha=0.8) + theme_bw() +
  scale_x_reordered() +
  facet_wrap(~ sEOF, scales="free_x") +
  xlab("Species") +
  scale_fill_manual(values= genus_colour_table$Colour,breaks = genus_colour_table$species) +
  scale_y_continuous(limits = c(0, 25000), labels = function(x) format(x, scientific = TRUE),
                     expand = expansion(mult = c(0, 0.5))) +
  theme(panel.background = element_blank(),#change background to white
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 1, colour = "black"),
        axis.ticks = element_line(size = 2),
        axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
        axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
        axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
        plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
        legend.title=element_text(size=27),
        legend.text=element_text(size=23),
        legend.direction="vertical",
        legend.position ="right") 
)
dev.off() 
##single barplot, split by sEOF
#arrange & mutate to order species by survival data only
tiff("specific_comparisons/sEOF/survivalEOF_rank_grouped_barplot.tiff",width = 700, height = 600)
print(
ranked_sp_sEOF.long %>%  
  arrange(sEOF, desc(Counts)) %>%
  mutate(Species = factor(Species, levels = unique(Species))) %>%
  ggplot() + aes(x=Species,y=Counts,fill=sEOF) +
  geom_col(position = "dodge", width=0.6,alpha=0.8) + 
  scale_fill_manual(values=c("#10639B","#808080")) +
  scale_y_continuous(limits = c(0, 25000), labels = function(x) format(x, scientific = TRUE),
                     expand = expansion(mult = c(0, 0.1))) +
  theme(panel.background = element_blank(),#change background to white
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 1, colour = "black"),
        axis.ticks = element_line(size = 2),
        axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
        axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
        axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
        plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
        legend.title=element_text(size=27),
        legend.text=element_text(size=23),
        legend.direction="vertical",
        legend.position ="right") 
)
dev.off()
#arrange & mutate to order species by no survival data only 
#must reorder sEOF factor first
ranked_sp_sEOF.long %>% mutate(sEOF = factor(sEOF, levels =c("no_survival","survival")))   %>%   
  arrange(sEOF, desc(Counts)) %>%
  mutate(Species = factor(Species, levels = unique(Species))) %>%
  ggplot() + aes(x=Species,y=Counts,fill=sEOF) +
  geom_col(position = "dodge", width=0.5,alpha=0.8) + 
  scale_y_continuous(limits = c(0, 25000), labels = function(x) format(x, scientific = TRUE),
                     expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values=c("#808080","#10639B")) +
  theme(panel.background = element_blank(),#change background to white
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 1, colour = "black"),
        axis.ticks = element_line(size = 2),
        axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
        axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
        axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
        plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
        legend.title=element_text(size=27),
        legend.direction="vertical",
        legend.position ="right") 

### Systemic AF and survival outcome #####

## merge data by survivalEOF status & plot ###
#create merged variable 
variable1 = as.character(get_variable(physeq_abund, "SurvivalEOF"))
variable2 = as.character(get_variable(physeq_abund, "SystemicAFbeforeatICUAdmin"))
sample_data(physeq_abund)$AFoutcome <- mapply(paste0, "S", variable1, "AF", variable2, 
                                              collapse = "_")
write.table(sample_data(physeq_abund),file="sample_data_with_AFoutcome.csv",sep=",",row.names = T)
#merged data by new variable
AFoutcome_data <- merge_samples(physeq_abund, "AFoutcome")
#remove the empty taxa 
remaining_taxa <- names(which(taxa_sums(AFoutcome_data)>0))

AFoutcome_data <- prune_taxa(remaining_taxa, AFoutcome_data)
#collapse genera
AFoutcome_data <- tax_glom(AFoutcome_data,taxrank="Genus_other")
#add sample names as a column in sample data
sample_data(AFoutcome_data)$name <- sample_names(AFoutcome_data)
#filter by new sample name column to remove those with unknown data (NAs)
AFoutcome_data <- subset_samples(AFoutcome_data, name=="S1AF1" |
                                   name=="S1AF0" |
                                   name=="S0AF1" |
                                   name=="S0AF0" )


#Transform to percentages of total available.
AFoutcome_data_perc = transform_sample_counts(AFoutcome_data, function(x) 100 * x/sum(x))
sample_names(AFoutcome_data_perc) <- c("No survival + no systemic AF","No survival + systemic AF","Survival + no systemic AF","Survival + systemic AF")
sample_names(AFoutcome_data) <- c("No survival + no systemic AF","No survival + systemic AF","Survival + no systemic AF","Survival + systemic AF")

#remove taxa at or below 0.1% in all merged samples (to keep plot legend simple)
AFoutcome_data_perc <- filter_taxa(AFoutcome_data_perc, function(x) sum(x > 0.1) > 0, TRUE)

#create proportional bar plot of merged data
svg("specific_comparisons/Genus_barplot_merged_AFoutcome.svg",width = 16, height = 7)
print(
  plot_bar(AFoutcome_data_perc, fill = "Genus_other") + coord_flip() + 
    ylab("Percentage of sequences") +
    xlab("AFoutcome") +
    geom_bar(aes(fill=Genus_other), stat="identity", position="stack") +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=27),
          legend.text=element_text(size=25, face="italic"),
          legend.direction="vertical",
          legend.position="bottom",
          strip.text.x = element_text(size = 27)) +
    guides(fill=guide_legend(ncol=3)) + #sets legend items to X column
    scale_fill_manual(values= genus_colour_table$Colour,breaks = genus_colour_table$species)
)
dev.off()


#plot heatmap 
tiff("specific_comparisons/AFoutcome_genus_heatmap.tiff",width = 700, height = 1200)
print(phyloseq::plot_heatmap(AFoutcome_data, taxa.label = "Genus_other", method = "NMDS",
                             distance = "bray", taxa.order = "Family",  low="beige",mid="blue", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical")
      
      
)
dev.off()

### Corticos comparisons ####
dir.create("specific_comparisons/Corticos")
## merge data by CAPA status & plot ###
physeq_abund_Corticos <- merge_samples(physeq_abund, "Corticos")

#Transform to percentages of total available.
physeq_abund_Corticos = transform_sample_counts(physeq_abund_Corticos, function(x) 100 * x/sum(x))

#create proportional bar plot of merged data
tiff("specific_comparisons/Corticos/Genus_barplot_merged_Corticos.tiff",width = 1200, height = 500)
print(
  plot_bar(physeq_abund_Corticos, fill = "Genus_other") + coord_flip() + 
    ylab("Percentage of sequences") +
    xlab("Corticosteroid use") +
    geom_bar(aes(fill=Genus_other), stat="identity", position="stack") +
    theme(panel.background = element_rect(fill = "white"), #change background to white
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=27),
          legend.text=element_text(size=25, face="italic"),
          legend.direction="vertical",
          legend.position="bottom",
          strip.text.x = element_text(size = 27)) +
    guides(fill=guide_legend(ncol=4)) + #sets legend items to X column
    scale_fill_manual(values= genus_colour_table$Colour,breaks = genus_colour_table$species)
)
dev.off()

## subset those with Corticos ###
Corticos_physeq <- subset_samples(physeq_abund_genera, Corticos==1)
#plot genera level taxa for each sample as heatmap
tiff("specific_comparisons/Corticos/Corticos_Genus_abundance_heatmap.tiff",width = 1200, height = 1000)
print(phyloseq::plot_heatmap(Corticos_physeq, taxa.label = "Genus_other", method = "NMDS",
                   distance = "bray", taxa.order = "Family",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical")
      
      
)
dev.off()

##subset Corticos species level data 
Corticos_physeq_sp <- subset_samples(physeq_abund, Corticos==1)

#plot Corticos rank abundance (average counts)
#set number of top taxa to plot
N <- 10
#plot
tiff("specific_comparisons/Corticos/Corticos_species_rank_mean_abundance.tiff",width = 1600, height = 1100)
par(mar=c(36,10,4,1)+.1)
barplot(sort(taxa_sums(Corticos_physeq_sp), TRUE)[1:N]/nsamples(Corticos_physeq_sp), las=2,cex.lab=2, cex.axis=3,cex.names=3,ylim = c(0,50000)) 
dev.off()

## subset data - those without Corticos ###
noCorticos_physeq <- subset_samples(physeq_abund_genera, Corticos==0)
#plot genera level taxa for each sample as heatmap
tiff("specific_comparisons/Corticos/noCorticos_Genus_abundance_heatmap.tiff",width = 700, height = 1000)
print(phyloseq::plot_heatmap(noCorticos_physeq, taxa.label = "Genus_other", method = "NMDS",
                   distance = "bray", taxa.order = "Family",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical")
      
      
)
dev.off()

##subset noCorticos species level data 
noCorticos_physeq_sp <- subset_samples(physeq_abund, Corticos==0)
#plot noCorticos rank abundance (average counts)
#set number of top taxa to plot
N <- 10
#plot
tiff("specific_comparisons/Corticos/noCorticos_species_rank_mean_abundance.tiff",width = 1600, height = 1100)
par(mar=c(36,10,4,1)+.1)
barplot(sort(taxa_sums(noCorticos_physeq_sp), TRUE)[1:N]/nsamples(noCorticos_physeq_sp), las=2,cex.lab=2, cex.axis=3,cex.names=3,ylim = c(0,50000)) 
dev.off()

##facet by cortico
#plot genera level taxa for each sample as barplot
tiff("specific_comparisons/Corticos/sample_genus_barplot_facet_Corticos.tiff",width = 1800, height = 800)
print(
  plot_bar(physeq_abund_genera, fill = "Genus_other") + 
    geom_bar(aes(fill=Genus_other), stat="identity", position="stack") +
    facet_grid(~Corticos,scales="free_x", space="free") +
    #theme(legend.position="none") +
    scale_fill_manual(values= genus_colour_table$Colour,breaks = genus_colour_table$species) +
    ylab("Abundance (Genus level)") + #set colours using colour table, set labels using formatted species names
    theme(panel.background = element_rect(fill = "white"), #change background to white
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=27),
          legend.text=element_text(size=25, face="italic"),
          legend.direction="vertical",
          legend.position="bottom",
          strip.text.x = element_text(size = 27)) +
    guides(fill=guide_legend(ncol=4)) #sets legend items to X column
  
)
dev.off()

### BAL GM 1 ODI comparisons ####
dir.create("specific_comparisons/BAL_GM1")
## merge data by BAL GM status & plot ###
physeq_abund_BALGM1ODI <- merge_samples(physeq_abund, "BALGM1ODI")

#Transform to percentages of total available.
physeq_abund_BALGM1ODI = transform_sample_counts(physeq_abund_BALGM1ODI, function(x) 100 * x/sum(x))

#create proportional bar plot of merged data
tiff("specific_comparisons/BAL_GM1/Genus_barplot_merged_BALGM1ODI.tiff",width = 1200, height = 500)
print(
  plot_bar(physeq_abund_BALGM1ODI, fill = "Genus_other") + coord_flip() + 
    ylab("Percentage of sequences") +
    xlab("BAL Galactomannan") +
    geom_bar(aes(fill=Genus_other), stat="identity", position="stack") +
    theme(panel.background = element_rect(fill = "white"), #change background to white
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=27),
          legend.text=element_text(size=25, face="italic"),
          legend.direction="vertical",
          legend.position="bottom",
          strip.text.x = element_text(size = 27)) +
    guides(fill=guide_legend(ncol=4)) + #sets legend items to X column
    scale_fill_manual(values= genus_colour_table$Colour,breaks = genus_colour_table$species)
)
dev.off()

#BALGM1ODI positive only 
BALGM1ODI_physeq_sp <- subset_samples(physeq_abund, BALGM1ODI==1) 
#plot rank abundance (average counts)
#set number of top taxa to plot
N <- 10
#plot
tiff("specific_comparisons/BAL_GM1/BALGM1ODI_species_rank_mean_abundance.tiff",width = 1600, height = 1100)
par(mar=c(36,10,4,1)+.1)
barplot(sort(taxa_sums(BALGM1ODI_physeq_sp), TRUE)[1:N]/nsamples(BALGM1ODI_physeq_sp), las=2,cex.lab=2, cex.axis=3,cex.names=3,ylim = c(0,50000)) 
dev.off()


#plot species level taxa for each sample as heatmap
tiff("specific_comparisons/BAL_GM1/BALGM1ODI_Species_abundance_heatmap.tiff",width = 700, height = 1300)
print(phyloseq::plot_heatmap(BALGM1ODI_physeq_sp, taxa.label = "Taxa", method = "NMDS",
                   distance = "bray", taxa.order = "Genus_other",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical")
)
dev.off()

##BALGM1ODI negative
BALGM1ODI_neg_physeq_sp <- subset_samples(physeq_abund, BALGM1ODI==0) 
#plot rank abundance (average counts)
#set number of top taxa to plot
N <- 10
#plot
tiff("specific_comparisons/BAL_GM1/BALGM1ODI_neg_species_rank_mean_abundance.tiff",width = 1600, height = 1100)
par(mar=c(36,10,4,1)+.1)
barplot(sort(taxa_sums(BALGM1ODI_neg_physeq_sp), TRUE)[1:N]/nsamples(BALGM1ODI_neg_physeq_sp), las=2,cex.lab=2, cex.axis=3,cex.names=3,ylim = c(0,50000)) 
dev.off()

#plot species level taxa for each sample as heatmap
tiff("specific_comparisons/BAL_GM1/BALGM1ODI_neg_Species_abundance_heatmap.tiff",width = 1200, height = 1400)
print(phyloseq::plot_heatmap(BALGM1ODI_neg_physeq_sp, taxa.label = "Taxa", method = "NMDS",
                   distance = "bray", taxa.order = "Genus_other",  low="beige", high="red", na.value="white") +
        theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical")
)
dev.off()

###Merge BALGMODI with survival EOF ####
#create merged variable 
variable1 = as.character(get_variable(physeq_abund, "SurvivalEOF"))
variable2 = as.character(get_variable(physeq_abund, "BALGM1ODI"))
sample_data(physeq_abund)$GM_sEOF <- mapply(paste0, "S", variable1, "GM", variable2, 
                                              collapse = "_")
#show sample sizes per group
count(sample_data(physeq_abund)$GM_sEOF)
write.table(sample_data(physeq_abund),file="sample_data_with_GM_sEOF.csv",sep=",",row.names = T)
#merged data by new variable
GM_sEOF_data <- merge_samples(physeq_abund, "GM_sEOF")
#remove the empty taxa 
remaining_taxa <- names(which(taxa_sums(GM_sEOF_data)>0))

GM_sEOF_data <- prune_taxa(remaining_taxa, GM_sEOF_data)
#add sample names as a column in sample data
sample_data(GM_sEOF_data)$name <- sample_names(GM_sEOF_data)
#filter by new sample name column to remove those with unknown data (NAs)
GM_sEOF_data <- subset_samples(GM_sEOF_data, name=="S1GM1" |
                                   name=="S1GM0" |
                                   name=="S0GM1" |
                                   name=="S0GM0" )



#Transform to percentages of total available.
GM_sEOF_data_perc = transform_sample_counts(GM_sEOF_data, function(x) 100 * x/sum(x))
sample_names(GM_sEOF_data_perc) <- c("No survival + negative BAL GM","No survival + positive BAL GM","Survival + negative BAL GM","Survival + positive BAL GM")
sample_names(GM_sEOF_data) <- c("No survival + negative BAL GM","No survival + positive BAL GM","Survival + negative BAL GM","Survival + positive BAL GM")

#remove taxa at or below 0.1% in all merged samples (to keep plot legend simple)
GM_sEOF_data_perc <- filter_taxa(GM_sEOF_data_perc, function(x) sum(x > 0.1) > 0, TRUE)

#create proportional bar plot of merged data
tiff("specific_comparisons/BAL_GM1/Genus_barplot_merged_GM_sEOF.tiff",width = 1400, height = 500)
print(
  plot_bar(GM_sEOF_data_perc, fill = "Genus_other") + coord_flip() + 
    ylab("Percentage of sequences") +
    xlab("BAL GM and survival outcome") +
    geom_bar(aes(fill=Genus_other), stat="identity", position="stack") +
    theme(panel.background = element_blank(),#change background to white
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=22),
          legend.text=element_text(size=18, face="italic"),
          legend.direction="vertical",
          legend.position="bottom",
          strip.text.x = element_text(size = 27)) +
    guides(fill=guide_legend(ncol=3)) + #sets legend items to X column
    scale_fill_manual(values= genus_colour_table$Colour,breaks = genus_colour_table$species)
)
dev.off()

###Merge BALGMODI with systemic antifungals ####
#create merged variable 
variable1 = as.character(get_variable(physeq_abund, "SystemicAFbeforeatICUAdmin"))
variable2 = as.character(get_variable(physeq_abund, "BALGM1ODI"))
sample_data(physeq_abund)$GM_sAF <- mapply(paste0, "AF", variable1, "GM", variable2, 
                                            collapse = "_")
#show sample sizes per group
count(sample_data(physeq_abund)$GM_sAF)
write.table(sample_data(physeq_abund),file="sample_data_with_GM_sAF.csv",sep=",",row.names = T)
#merged data by new variable
GM_sAF_data <- merge_samples(physeq_abund, "GM_sAF")
#remove the empty taxa 
remaining_taxa <- names(which(taxa_sums(GM_sAF_data)>0))

GM_sAF_data <- prune_taxa(remaining_taxa, GM_sAF_data)
#add sample names as a column in sample data
sample_data(GM_sAF_data)$name <- sample_names(GM_sAF_data)
#filter by new sample name column to remove those with unknown data (NAs)
GM_sAF_data <- subset_samples(GM_sAF_data, name=="AF1GM1" |
                                 name=="AF1GM0" |
                                 name=="AF0GM1" |
                                 name=="AF0GM0" )



#Transform to percentages of total available.
GM_sAF_data_perc = transform_sample_counts(GM_sAF_data, function(x) 100 * x/sum(x))
sample_names(GM_sAF_data_perc) <- c("No antifungals + negative BAL GM","No antifungals + positive BAL GM","Antifungals + negative BAL GM")
sample_names(GM_sAF_data) <- c("No antifungals + negative BAL GM","No antifungals + positive BAL GM","Antifungals + negative BAL GM")

#remove taxa at or below 0.1% in all merged samples (to keep plot legend simple)
GM_sAF_data_perc <- filter_taxa(GM_sAF_data_perc, function(x) sum(x > 0.1) > 0, TRUE)

#create proportional bar plot of merged data
tiff("specific_comparisons/BAL_GM1/Genus_barplot_merged_GM_sAF.tiff",width = 1400, height = 420)
print(
  plot_bar(GM_sAF_data_perc, fill = "Genus_other") + coord_flip() + 
    ylab("Percentage of sequences") +
    xlab("Antifungals and BAL GM") +
    geom_bar(aes(fill=Genus_other), stat="identity", position="stack") +
    theme(panel.background = element_blank(),#change background to white
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=22),
          legend.text=element_text(size=18, face="italic"),
          legend.direction="vertical",
          legend.position="bottom",
          strip.text.x = element_text(size = 27)) +
    guides(fill=guide_legend(ncol=3)) + #sets legend items to X column
    scale_fill_manual(values= genus_colour_table$Colour,breaks = genus_colour_table$species)
)
dev.off()

#summarise group sizes for above 
sink("specific_comparisons/group_sizes_summary.txt")
paste("noSurvivalEOF_physeq:",nsamples(noSurvivalEOF_physeq))
paste("SurvivalEOF_physeq:",nsamples(SurvivalEOF_physeq))
paste("SystemicAFbeforeatICUAdmin_physeq:",nsamples(SystemicAFbeforeatICUAdmin_physeq))
paste("noSystemicAFbeforeatICUAdmin_physeq:",nsamples(noSystemicAFbeforeatICUAdmin_physeq))
paste("noCAPA_physeq:",nsamples(noCAPA_physeq))
paste("CAPA_physeq:",nsamples(CAPA_physeq))
paste("noCorticos_physeq:",nsamples(noCorticos_physeq))
paste("Corticos_physeq:",nsamples(Corticos_physeq))
sink()

### Sample source ####
### Influence of sample source 
#samples ordination
tiff("specific_comparisons/sample_ordination_Institution.tiff",width = 700, height = 500)
plot_ordination(physeq_abund, physeq.ord, type="samples", color="Institution", 
                title="Samples") + geom_point(size=3)
dev.off()


### test variety of variables #### ################
dir.create("general_variable_check")

### Pairwise comparisons ####
#read in clinical data - only the pairwise comparison variables 
all_samples_pw <- read_excel("clinical_data/COVID_clinical_data_pairwise_comparison_variables_only.xlsx", sheet = "samples")

#convert to character
all_samples_pw[,1:25]  <- lapply(all_samples_pw[,1:25] , as.character)
#subset clinical data for samples available 
#alternatively, can subset later using subset_samples(physeq, <sample data column> =="Yes")
samples_df_pw <- subset(all_samples_pw, Sample %in% colnames(otu_mat) )

###check group size for each comparison 
grp.sizes <- lapply(samples_df_pw[,c(1,3:25)] , plyr::count)
grop.sizes <- do.call(rbind, grp.sizes)
grop.sizes$var <- rownames(grop.sizes)
grop.sizes$var <- grop.sizes$var %>% str_remove("\\.[1-3]")
library(tidyr)
group.sizes <- spread(grop.sizes,x,freq)
#add column reporting if group sizes are both more than 5 samples
group.sizes$size_pass <- group.sizes$`0` >= 5 & group.sizes$`1` >= 5 
write.table(group.sizes,file="group_sizes.csv",sep=",",row.names = F)
#Assign rownames for tax
samples_df_pw <- samples_df_pw %>% 
  tibble::column_to_rownames("Sample") 

#### Diversity tests & stats  ####
#create richness estimate data
Obs <- estimate_richness(physeq,measures="Observed")
Chao <- estimate_richness(physeq,measures="Chao1")
Shann <- estimate_richness(physeq,measures="Shannon")

Obs_abund <- estimate_richness(physeq_abund,measures="Observed")
Chao_abund <- estimate_richness(physeq_abund,measures="Chao1")
Shann_abund <- estimate_richness(physeq_abund,measures="Shannon")


#prep data for performing wilcox
#pre data - raw 
Obs$Sample <- rownames(Obs)
Chao$Sample <- rownames(Chao)
Shann$Sample <- rownames(Shann)
samples_df_pw$Sample <- rownames(samples_df_pw)
ri_samples_df_pw  <- merge(Obs,samples_df_pw, by="Sample")
ric_samples_df_pw  <- merge(Chao,ri_samples_df_pw, by="Sample")
rich_samples_df_pw  <- merge(Shann,ric_samples_df_pw, by="Sample")

#prep data - filtered
Obs_abund$Sample <- rownames(Obs_abund)
Chao_abund$Sample <- rownames(Chao_abund)
Shann_abund$Sample <- rownames(Shann_abund)
samples_df_pw$Sample <- rownames(samples_df_pw)
ri_samples_df_pw_abund  <- merge(Obs_abund,samples_df_pw, by="Sample")
ric_samples_df_pw_abund  <- merge(Chao_abund,ri_samples_df_pw_abund, by="Sample")
rich_samples_df_pw_abund  <- merge(Shann_abund,ric_samples_df_pw_abund, by="Sample")

#get mean and median values grouped by variables 
#CAPA status
rich_samples_df_pw_abund %>%  dplyr::filter(CAPANoCAPA != "<NA>") %>% 
  group_by(CAPANoCAPA) %>% summarise(Med.Shannon=median(Shannon,na.rm=TRUE),
                                     Med.Chao1=median(Chao1,na.rm=TRUE), 
                                     Med.Observed=median(Observed,na.rm=TRUE))

rich_samples_df_pw_abund %>%  dplyr::filter(CAPANoCAPA != "<NA>") %>% 
  group_by(CAPANoCAPA) %>% summarise(Mean.Shannon=mean(Shannon,na.rm=TRUE),
                                     Mean.Chao1=mean(Chao1,na.rm=TRUE), 
                                     Mean.Observed=mean(Observed,na.rm=TRUE))
#Corticos
rich_samples_df_pw_abund %>%  dplyr::filter(Corticos != "<NA>") %>% 
  group_by(Corticos) %>% summarise(Med.Shannon=median(Shannon,na.rm=TRUE),
                                     Med.Chao1=median(Chao1,na.rm=TRUE), 
                                     Med.Observed=median(Observed,na.rm=TRUE))

rich_samples_df_pw_abund %>%  dplyr::filter(Corticos != "<NA>") %>% 
  group_by(Corticos) %>% summarise(Mean.Shannon=mean(Shannon,na.rm=TRUE),
                                     Mean.Chao1=mean(Chao1,na.rm=TRUE), 
                                     Mean.Observed=mean(Observed,na.rm=TRUE))
#sAF
rich_samples_df_pw_abund %>%  dplyr::filter(SystemicAFbeforeatICUAdmin != "<NA>") %>% 
  group_by(SystemicAFbeforeatICUAdmin) %>% summarise(Med.Shannon=median(Shannon,na.rm=TRUE),
                                   Med.Chao1=median(Chao1,na.rm=TRUE), 
                                   Med.Observed=median(Observed,na.rm=TRUE))

rich_samples_df_pw_abund %>%  dplyr::filter(SystemicAFbeforeatICUAdmin != "<NA>") %>% 
  group_by(SystemicAFbeforeatICUAdmin) %>% summarise(Mean.Shannon=mean(Shannon,na.rm=TRUE),
                                   Mean.Chao1=mean(Chao1,na.rm=TRUE), 
                                   Mean.Observed=mean(Observed,na.rm=TRUE))
#sEOF
rich_samples_df_pw_abund %>%  dplyr::filter(SurvivalEOF != "<NA>") %>% 
  group_by(SurvivalEOF) %>% summarise(Med.Shannon=median(Shannon,na.rm=TRUE),
                                                     Med.Chao1=median(Chao1,na.rm=TRUE), 
                                                     Med.Observed=median(Observed,na.rm=TRUE))

rich_samples_df_pw_abund %>%  dplyr::filter(SurvivalEOF != "<NA>") %>% 
  group_by(SurvivalEOF) %>% summarise(Mean.Shannon=mean(Shannon,na.rm=TRUE),
                                                     Mean.Chao1=mean(Chao1,na.rm=TRUE), 
                                                     Mean.Observed=mean(Observed,na.rm=TRUE))

## pairwise Wilcox ###
#prep empty df to fill with pvals - raw
mat <- matrix(ncol=4,nrow=24)
df <- data.frame(mat)
colnames(df) <- c("Var","Obs","Shannon","Chao1")
df$Var <- colnames(rich_samples_df_pw)[6:29]

#prep empty df to fill with pvals - filtered
mat2 <- matrix(ncol=4,nrow=24)
df2 <- data.frame(mat2)
colnames(df2) <- c("Var","Obs","Shannon","Chao1")
df2$Var <- colnames(rich_samples_df_pw_abund)[6:29]

dir.create("general_variable_check/wilcox")
#report any sig. pairwise comparison stats results in text file
sink("general_variable_check/wilcox/wilcox_summary.txt")
#loop through and perform pairwise comparisons
#report if p val < 0.05
for (i in seq(6,length(colnames(rich_samples_df_pw)))){
  col <- colnames(rich_samples_df_pw)[i]
  ###raw data
  o.wil <- pairwise.wilcox.test(rich_samples_df_pw$Observed, rich_samples_df_pw[[col]], exact=F)
  #add pval to df
  df$Obs[which(df$Var==col)] <- o.wil$p.value
  #note if significant
  if (o.wil$p.value < 0.05 & o.wil$p.value != 0 ){
    print(paste("Observed OTUs pval is sig for ",col, ".(raw data) pval = ",o.wil$p.value, sep=""))
  } 
  s.wil <- pairwise.wilcox.test(rich_samples_df_pw$Shannon, rich_samples_df_pw[[col]])
  #add pval to df
  df$Shannon[which(df$Var==col)] <- s.wil$p.value
  #note if significant
  if (s.wil$p.value < 0.05 & s.wil$p.value != 0 ){
    print(paste("Shannon pval is sig for ",col, ".(raw data) pval = ",s.wil$p.value, sep=""))
  } 
  c.wil <- pairwise.wilcox.test(rich_samples_df_pw$Chao1, rich_samples_df_pw[[col]], exact=F)
  #add pval to df
  df$Chao1[which(df$Var==col)] <- c.wil$p.value
  #note if significant
  if (c.wil$p.value < 0.05 & c.wil$p.value != 0 ){
    print(paste("Chao1 pval is sig for ",col, ".(raw data) pval = ",c.wil$p.value, sep=""))
  } 
  ##remove NA samples for plotting
  var_values <- sample_data(physeq)[[col]]
  clean <- prune_samples(!is.na(var_values), physeq)
  #Grouped alpha diversity on input data (not normalised or taxa filtered)
  tiff(paste("general_variable_check/richness_",col,"_grouped_raw.tiff",sep=""),width = 500, height = 400)
  print(plot_richness(clean,measures=c("Chao1", "Shannon","Observed"), x=col)+ geom_boxplot() +
          theme(panel.background = element_rect(fill = "white"), #change background to white
                axis.line = element_line(size = 1, colour = "black"),
                axis.ticks = element_line(size = 2),
                axis.text.x = element_text(color = "grey20", size = 30, angle =90, hjust = .5, vjust = .5, face = "plain"),
                axis.text.y = element_text(color = "grey20", size = 30, angle = 0, hjust = 1, vjust = .5, face = "plain"),
                axis.title.x = element_text(face = "bold", color = "grey20", size = 32, angle = 0, hjust = .5, vjust = 0),
                axis.title.y = element_text(face = "bold", color = "grey20", size = 32, angle = 90, hjust = .5, vjust = .5),
                plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
                legend.title=element_text(size=27),
                legend.text=element_text(size=25, face="italic"),
                legend.direction="vertical",
                legend.position="bottom",
                strip.text.x = element_text(size = 32, colour = "black"))
        
  )
  dev.off()
  
  ###filtered data ###
  o.wil2 <- pairwise.wilcox.test(rich_samples_df_pw_abund$Observed, rich_samples_df_pw_abund[[col]], exact=F)
  #add pval to df
  df2$Obs[which(df2$Var==col)] <- o.wil2$p.value
  #note if significant
  if (o.wil2$p.value < 0.05 & o.wil2$p.value != 0 ){
    print(paste("Observed OTUs pval is sig for ",col, ".(filtered data) pval = ",o.wil2$p.value, sep=""))
  } 
  s.wil2 <- pairwise.wilcox.test(rich_samples_df_pw_abund$Shannon, rich_samples_df_pw_abund[[col]])
  #add pval to df
  df2$Shannon[which(df2$Var==col)] <- s.wil2$p.value
  #note if significant
  if (s.wil2$p.value < 0.05 & s.wil2$p.value != 0 ){
    print(paste("Shannon pval is sig for ",col, ".(filtered data) pval = ",s.wil2$p.value, sep=""))
  } 
  c.wil2 <- pairwise.wilcox.test(rich_samples_df_pw_abund$Chao1, rich_samples_df_pw_abund[[col]], exact=F)
  #add pval to df
  df2$Chao1[which(df2$Var==col)] <- c.wil2$p.value
  #note if significant
  if (c.wil2$p.value < 0.05 & c.wil2$p.value != 0 ){
    print(paste("Chao1 pval is sig for ",col, ".(filtered data) pval = ",c.wil2$p.value, sep=""))
  } 
  ##remove NA samples for plotting
  var_values <- sample_data(physeq_abund)[[col]]
  clean_abund <- prune_samples(!is.na(var_values), physeq_abund)
  # Grouped alpha diversity on normalised & taxa filtered data 
  #coloured plots as svg
  svg(paste("general_variable_check/richness_",col,"_grouped_abund.svg",sep=""),width = 8, height = 9)
  print(plot_richness(clean_abund,measures=c("Chao1", "Shannon","Observed"), x=col,color=col)+ geom_boxplot()+
          #scale_color_manual(values=c("#10639B","#808080"))   +
          #(values=c("#10639B","#808080")) blue/grey
          scale_color_manual(values=c("#7570B3","#E6AB02")) + #purple/yellow
          theme(panel.background = element_blank(), #change background to white
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                axis.line = element_line(size = 1, colour = "black"),
                axis.ticks = element_line(size = 2),
                axis.text.x = element_blank(),
                axis.text.y = element_text(color = "grey20", size = 40, angle = 0, hjust = 1, vjust = .5, face = "plain"),
                axis.title.x = element_text(face = "bold", color = "grey20", size = 32, angle = 0, hjust = .5, vjust = 0),
                axis.title.y = element_blank(),
                plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
                legend.title=element_text(size=27),
                legend.text=element_text(size=25, face="italic"),
                legend.direction="vertical",
                legend.position="bottom",
                strip.text.x = element_text(size = 27, colour = "black"))
  )
  dev.off()
  #uncoloured plots as svg
  svg(paste("general_variable_check/richness_",col,"_grouped_abund_uncoloured.svg",sep=""),width = 12, height = 9)
  print(plot_richness(clean_abund,measures=c("Chao1", "Shannon","Observed"), x=col)+ geom_boxplot()+
          theme(panel.background = element_blank(), #change background to white
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                axis.line = element_line(size = 1, colour = "black"),
                axis.ticks = element_line(size = 2),
                axis.text.x = element_blank(),
                axis.text.y = element_text(color = "grey20", size = 40, angle = 0, hjust = 1, vjust = .5, face = "plain"),
                axis.title.x = element_text(face = "bold", color = "grey20", size = 32, angle = 0, hjust = .5, vjust = 0),
                axis.title.y = element_blank(),
                plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
                legend.title=element_text(size=27),
                legend.text=element_text(size=25, face="italic"),
                legend.direction="vertical",
                legend.position="bottom",
                strip.text.x = element_text(size = 27, colour = "black"))
  )
  dev.off()
}
sink()
#write pvals to file 
write.csv(df,"general_variable_check/wilcox/wilcox_pvals.csv", row.names=F)
write.csv(df2,"general_variable_check/wilcox/wilcox_pvals_abund.csv", row.names=F)

#### Ordination permanovas ####
dir.create("general_variable_check/ordination_permanova")

wunifrac_dist = phyloseq::distance(physeq_abund, method="bray", weighted=F)
ordination = ordinate(physeq_abund, method="PCoA", distance=wunifrac_dist)
sink("general_variable_check/ordination_permanova/bray_perma_summary.txt")
##Loop through all variables for ordination permanovas
for (i in seq(6,length(colnames(rich_samples_df_pw)))){
  col <- colnames(rich_samples_df_pw)[i]
  
  ###permanova test on bray ordination ###
  tiff(paste("general_variable_check/ordination_permanova/bray_perma_",col,"_grouped_abund.tiff",sep=""),
       width = 700, height = 600)
  print(plot_ordination(physeq_abund, ordination, color=col) + theme(aspect.ratio=1) + geom_point(size=3) +
    theme(panel.background = element_blank(), #change background to white
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 30, angle =0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 30, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 32, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 32, angle = 90, hjust = .5, vjust = 0),
          legend.title=element_text(size=27),
          legend.text=element_text(size=25, face="italic"),
          legend.direction="vertical",
          legend.position = "bottom",
          strip.text.x = element_text(size = 32, colour = "black"))
  )
  dev.off()
  
  library(vegan)
  #### tests ##
  ##remove NA samples for plotting
  var_values <- sample_data(physeq_abund)[[col]]
  clean2 <- prune_samples(!is.na(var_values), physeq_abund)
  
  data <- phyloseq::distance(clean2, method = "bray")
  f <- as.formula(paste("data~ ",col, sep=""))
  metadata <- as(sample_data(clean2), "data.frame")
  res.permanova <- adonis2(f,
                         data = metadata)
  write.table(res.permanova,file=paste("general_variable_check/ordination_permanova/",col,"_permanova_results.csv",sep=""),sep=",", row.names = F)
  if (res.permanova$`Pr(>F)`[1] < 0.05 ) {
    print(paste("Permanova pval is sig for ",col, ".(filtered data) pval = ",res.permanova$`Pr(>F)`, sep=""))
  } else {
    print(paste(" no significance for ", col, sep=""))
  }
}
sink()


#### Deseq - species level ####
library(DESeq2)
library(randomForest)
#set p val cutoff 
alpha = 0.05
#set basemean cutoff
bM = 500
####  Deseq - unnormalised but species filtered data - species ##
#loop through and perform pairwise comparisons
#report if p val < 0.05
dir.create("general_variable_check/deseq")
sink(file="general_variable_check/deseq/no_hit_variables.txt")
for (i in seq(6,length(colnames(rich_samples_df_pw)))){
  #set variable to test
  col <- colnames(rich_samples_df_pw)[i]
  #To test the differences at OTU level between a variable using DESeq2, we need to convert the variable column into factor
  sample_data(physeq_filt)[[col]] <- as.factor(sample_data(physeq_filt)[[col]])
  ##first need to remove samples with 'NA'
  #switched to prune_samples so that can feed in a variable
  var_values <- sample_data(physeq_filt)[[col]]
  clean <- prune_samples(!is.na(var_values), physeq_filt)
  #convert to deseq
  f <- as.formula(paste("~ ",col, sep=""))
  ds <- phyloseq_to_deseq2(clean, f)
  # calculate geometric means prior to estimate size factors (because there are 0s in OTU table)
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(ds), 1, gm_mean)
  ds = estimateSizeFactors(ds, geoMeans = geoMeans)
  ds = DESeq(ds, fitType="local")
  
  res = results(ds, cooksCutoff = FALSE)
  #add taxa names to results table
  res = cbind(as(rownames(res), "matrix"), as(res, "data.frame"), as(tax_table(physeq_filt)[rownames(res), ], "matrix"))
  colnames(res)[1] <- "Taxa"
  write.table(res,file=paste("general_variable_check/deseq/res_",col,".csv",sep=""),sep=",",row.names = F,col.names = T)
  
  #get sig hits
  res[which(res$padj < alpha ), ]
  sig = res[which(res$padj < alpha & res$baseMean > bM ), ]
  if (nrow(sig) > 0) {
    #add variable and taxa names to significant results table
    sig = cbind(as(rep(col,nrow(sig)), "matrix"), as(rownames(sig), "matrix"), as(sig, "data.frame"), as(tax_table(physeq_filt)[rownames(sig), ], "matrix"))
    colnames(sig)[1] <- "Variable"
    colnames(sig)[2] <- "Taxa"
    write.table(sig,file=paste("general_variable_check/deseq/",col,"deseq.csv",sep=""),sep=",",row.names = F,col.names = T)
  } else {
    print(paste(" no significant hits for ", col, sep=""))
  }
  
}
sink()

##read in all deseq sig hits and plot - standard pval
setwd("general_variable_check/deseq/")
sig_hits_files = (list.files(pattern="deseq.csv"))

#combine into table
tbl <- lapply(sig_hits_files, function(x){
  read.table(x, sep=",", as.is=TRUE, header=FALSE)
})
new_tbl <- do.call(rbind, tbl)

#remove duplicates rows
new_tbl2 <- new_tbl[!duplicated(new_tbl),]
colnames(new_tbl2) <- new_tbl2[1,]
new_tbl2 <- new_tbl2[-1,]
new_tbl2 <- new_tbl2[,-c(2:3)]
new_tbl2[,c(2:7)] <- lapply(new_tbl2[,c(2:7)], as.numeric)

unique(new_tbl2$Variable)
subset
#plot each set seperately (axis free)
for (i in seq(1,length(unique(new_tbl2$Variable)))){
  print(i)
  var <- unique(new_tbl2$Variable)[i]
  subset <- new_tbl2[which(new_tbl2$Variable == var),]
  subset$Taxa <- subset$Taxa %>% str_replace("_"," ")
  H <- (nrow(subset) * .5) + 2.5
  svg(paste(var, "sig_hits.svg", sep=""),height = 8, width = H)
  print(
    ggplot(subset,aes(x=Taxa,y=log2FoldChange)) +
      geom_bar(stat="identity", width =0.2) + geom_hline(yintercept = 0) + theme_minimal()  +
      # scale_y_continuous(limits=c(-30, 30), breaks=c(-30,-15,0,15,30)) +
      #coord_flip() + 
      theme(panel.background = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.border = element_blank(),
            axis.line = element_line(size = 1, colour = "black"),
            axis.ticks = element_line(size = 2),
            axis.text.x = element_text(color = "grey20", size = 25, angle =60, hjust = 1, vjust = 1, face = "italic"),
            axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
            axis.title.x = element_text(face = "bold", color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0),
            axis.title.y = element_text(face = "bold", color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5),
            legend.position ="left",
            plot.margin = ggplot2::margin(1,1,1,3, "cm"))
  
  )
  dev.off()
  
}
#plot each set seperately (axis free) - coord flip
for (i in seq(1,length(unique(new_tbl2$Variable)))) {
  print(i)
  var <- unique(new_tbl2$Variable)[i]
  subset <- new_tbl2[which(new_tbl2$Variable == var),]
  subset$Taxa <- subset$Taxa %>% str_replace("_"," ")
  H <- (nrow(subset) * .5) + 1
  svg(paste(var, "sig_hits_cflip.svg", sep=""),height = H, width = 7)
  print(
    ggplot(subset,aes(x=Taxa,y=log2FoldChange)) +
      geom_bar(stat="identity", width =0.2) + geom_hline(yintercept = 0) + theme_minimal()  +
      # scale_y_continuous(limits=c(-30, 30), breaks=c(-30,-15,0,15,30)) +
      coord_flip() + 
      theme(panel.background = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.border = element_blank(),
            axis.line = element_line(size = 1, colour = "black"),
            axis.ticks = element_line(size = 2),
            axis.text.x = element_text(color = "grey20", size = 25, angle =0, hjust = 1, vjust = 1, face = "plain"),
            axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "italic"),
            axis.title.x = element_text(face = "bold", color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0),
            axis.title.y = element_text(face = "bold", color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5),
            legend.position ="right")
    
  )
  dev.off()
}
  #plot each set seperately (axis free) - coord flip & set axis
  for (i in seq(1,length(unique(new_tbl2$Variable)))){
    print(i)
    var <- unique(new_tbl2$Variable)[i]
    subset <- new_tbl2[which(new_tbl2$Variable == var),]
    subset$Taxa <- subset$Taxa %>% str_replace("_"," ")
    H <- (nrow(subset) * 25) + 100
    tiff(paste(var, "sig_hits_cflip_set.tiff", sep=""),height = H, width = 550)
    print(
      ggplot(subset,aes(x=Taxa,y=log2FoldChange)) +
        geom_bar(stat="identity", width =0.2) + geom_hline(yintercept = 0) + theme_minimal()  +
        scale_y_continuous(limits=c(-30, 30), breaks=c(-30,-15,0,15,30)) +
        coord_flip() + 
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_blank(),
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =0, hjust = 1, vjust = 1, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "italic"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5),
              legend.position ="right")
      
    )
    dev.off()
}
  
#plot each set seperately (set axis)
for (i in seq(1,length(unique(new_tbl2$Variable)))){
  print(i)
  var <- unique(new_tbl2$Variable)[i]
  subset <- new_tbl2[which(new_tbl2$Variable == var),]
  H <- (nrow(subset) * 25) + 75
  tiff(paste(var, "sig_hits_set_axis.tiff", sep=""),height = 500, width = H)
  print(
    ggplot(subset,aes(x=Taxa,y=log2FoldChange)) +
      geom_bar(stat="identity", width =0.2) + geom_hline(yintercept = 0) + theme_minimal()  +
      scale_y_continuous(limits=c(-30, 30), breaks=c(-30,-15,0,15,30)) +
      #coord_flip() + 
      theme(panel.background = element_rect(fill = "white"), #change background to white
            axis.line = element_line(size = 1, colour = "black"),
            axis.ticks = element_line(size = 2),
            axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = 1, vjust = .5, face = "italic"),
            axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
            axis.title.x = element_text(face = "bold", color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0),
            axis.title.y = element_text(face = "bold", color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5),
            strip.text.x = element_text(size = 18)) 
  )
  dev.off()
  
}

#remove variables with insufficient group sizes 
new_tbl3 <- new_tbl2 %>% dplyr::filter(Variable != "BAL_LFD") %>% 
  dplyr::filter(Variable != "ECMOnew") %>% dplyr::filter(Variable != "PosBALculture") %>%
  dplyr::filter(Variable != "PosBALPCR") %>% dplyr::filter(Variable != "SerumGM05") %>%
  dplyr::filter(Variable != "Tocilizumab") %>% dplyr::filter(Variable != "UnderlDM") %>%
  dplyr::filter(Variable != "UnderlHemOnc") %>% dplyr::filter(Variable != "UnderlSOT")

#plot all results for variables with sufficient group sizes
tiff("sig_hits.tiff",width = 800, height = 1200)
print(
ggplot(new_tbl3,aes(x=Taxa,y=log2FoldChange)) +
  geom_bar(stat="identity", width =0.2) + geom_hline(yintercept = 0) + theme_minimal()  +
  #scale_y_continuous(limits=c(-30, 30), breaks=c(-30,-15,0,15,30)) +
  facet_wrap(~ Variable, nrow=4, scales="free_y") + coord_flip() + 
  theme(panel.background = element_rect(fill = "white"), #change background to white
      axis.line = element_line(size = 1, colour = "black"),
      axis.ticks = element_line(size = 2),
      axis.text.x = element_text(color = "grey20", size = 18, angle =90, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "grey20", size = 18, angle = 0, hjust = 1, vjust = .5, face = "plain"),
      axis.title.x = element_text(face = "bold", color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0),
      axis.title.y = element_text(face = "bold", color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5),
      strip.text.x = element_text(size = 18)) 
)
dev.off()

write.table(new_tbl2,file="full_deseq_hits.csv",sep=",",row.names = F)

setwd("../..")

# ##assess basemeans to choose cutoff
# bMs <- read.table("deseq/basemeans.csv")
# class(bMs$V1) <- "numeric"
# 
# median(bMs$V1)
# mean(bMs$V1)
# ggplot(bMs,aes(V1)) + geom_histogram(binwidth=10,alpha=5) 

#### Deseq - genus level ####
# unnormalised but species filtered data  
#loop through and perform pairwise comparisons
#report if p val < 0.05
dir.create("general_variable_check/deseq/genus")
sink(file="general_variable_check/deseq/genus/no_hit_variables_Gen.txt")
for (i in seq(6,length(colnames(rich_samples_df_pw)))){
  #set variable to test
  col <- colnames(rich_samples_df_pw)[i]
  #To test the differences at OTU level between a variable using DESeq2, we need to convert the variable column into factor
  sample_data(physeq_filt_genera)[[col]] <- as.factor(sample_data(physeq_filt_genera)[[col]])
  ##first need to remove samples with 'NA'
  #switched to prune_samples so that can feed in a variable
  var_values <- sample_data(physeq_filt_genera)[[col]]
  clean <- prune_samples(!is.na(var_values), physeq_filt_genera)
  #convert to deseq
  f <- as.formula(paste("~ ",col, sep=""))
  ds <- phyloseq_to_deseq2(clean, f)
  # calculate geometric means prior to estimate size factors (because there are 0s in OTU table)
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(ds), 1, gm_mean)
  ds = estimateSizeFactors(ds, geoMeans = geoMeans)
  ds = DESeq(ds, fitType="local")
  
  res = results(ds, cooksCutoff = FALSE)
  #add taxa names to results table
  res = cbind(as(rownames(res), "matrix"), as(res, "data.frame"), as(tax_table(physeq_filt_genera)[rownames(res), ], "matrix"))
  colnames(res)[1] <- "Taxa"
  write.table(res,file=paste("general_variable_check/deseq/genus/res_",col,".csv",sep=""),sep=",",row.names = F,col.names = T)
  
  #get sig hits
  sig = res[which(res$padj < alpha & res$baseMean > bM ), ]
  if (nrow(sig) > 0) {
    #add variable and taxa names to significant results table
    sig = cbind(as(rep(col,nrow(sig)), "matrix"), as(rownames(sig), "matrix"), as(sig, "data.frame"), as(tax_table(physeq_filt_genera)[rownames(sig), ], "matrix"))
    colnames(sig)[1] <- "Variable"
    colnames(sig)[2] <- "Taxa"
    head(sig)
    write.table(sig,file=paste("general_variable_check/deseq/genus/",col,"deseqGen.csv",sep=""),sep=",",row.names = F,col.names = T)
  } else {
    print(paste(" no significant genera hits for ", col, sep=""))
  }
  
}
sink()

#### ANCOMBC #####
#code adapted from: https://www.nicholas-ollberding.com/post/identifying-differentially-abundant-features-in-microbiome-data/
#BiocManager::install("ANCOMBC")
library("ANCOMBC")
dir.create("general_variable_check/ancom")
dir.create("general_variable_check/ancom/genus")

##loop to run ancom for all pairwise comparisons 
sink(file="general_variable_check/ancom/ancom_no_hits.txt")
for (i in seq(6,length(colnames(rich_samples_df_pw)))){
  #set variable to test
  col <- colnames(rich_samples_df_pw)[i]
  #To test the differences at OTU level between a variable, we need to convert the variable column into factor
  sample_data(physeq_filt)[[col]] <- as.factor(sample_data(physeq_filt)[[col]])

  #run ancombc
  ancom_da = ancombc(phyloseq = physeq_filt, formula = col,
                     p_adj_method = "holm", zero_cut = 0.90, lib_cut = 500,
                     group = col, struc_zero = FALSE, neg_lb = FALSE,
                     tol = 1e-5, max_iter = 100, conserve = TRUE,
                     alpha = 0.05, global = TRUE)
  ancom_res_df <- data.frame(
    Species = row.names(ancom_da$res$beta),
    beta = unlist(ancom_da$res$beta),
    se = unlist(ancom_da$res$se),
    W = unlist(ancom_da$res$W),
    p_val = unlist(ancom_da$res$p_val),
    q_val = unlist(ancom_da$res$q_val),
    diff_abn = unlist(ancom_da$res$diff_abn))
  #get sig hits
  fdr_ancom <- ancom_res_df %>%
    dplyr::filter(q_val < alpha)
  if (nrow(fdr_ancom) > 0) {
    write.table(fdr_ancom,file=paste("general_variable_check/ancom/fdr_",col,".csv",sep=""),sep=",",row.names = F,col.names = T)
    write.table(ancom_res_df,file=paste("general_variable_check/ancom/anc_res_",col,".csv",sep=""),sep=",",row.names = F,col.names = T)
  } else {
    print(paste(" no significant hits for ", col, sep=""))
    write.table(ancom_res_df,file=paste("general_variable_check/ancom/anc_res_",col,".csv",sep=""),sep=",",row.names = F,col.names = T)
    
  }
  
  
}
sink()


##loop to run ancom for all pairwise comparisons  - genus
sink(file="general_variable_check/ancom/genus/ancom_no_hits_Gen.txt")
for (i in seq(6,length(colnames(rich_samples_df_pw)))){
  #set variable to test
  col <- colnames(rich_samples_df_pw)[i]
  #To test the differences at OTU level between a variable, we need to convert the variable column into factor
  sample_data(physeq_filt_genera)[[col]] <- as.factor(sample_data(physeq_filt_genera)[[col]])
  #run ancom
  ancom_da = ancombc(phyloseq = physeq_filt_genera, formula = col,
                     p_adj_method = "holm", zero_cut = 0.90, lib_cut = 500,
                     group = col, struc_zero = FALSE, neg_lb = FALSE,
                     tol = 1e-5, max_iter = 100, conserve = TRUE,
                     alpha = 0.05, global = TRUE)
  
  ancom_res_df <- data.frame(
    Genus = row.names(ancom_da$res$beta),
    beta = unlist(ancom_da$res$beta),
    se = unlist(ancom_da$res$se),
    W = unlist(ancom_da$res$W),
    p_val = unlist(ancom_da$res$p_val),
    q_val = unlist(ancom_da$res$q_val),
    diff_abn = unlist(ancom_da$res$diff_abn))

  #get sig hits
  fdr_ancom <- ancom_res_df %>%
    dplyr::filter(q_val < alpha)
  if (nrow(fdr_ancom) > 0) {
    write.table(fdr_ancom,file=paste("general_variable_check/ancom/genus/fdr_",col,".csv",sep=""),sep=",",row.names = F,col.names = T)
    write.table(ancom_res_df,file=paste("general_variable_check/ancom/genus/anc_res_",col,".csv",sep=""),sep=",",row.names = F,col.names = T)
  } else {
    print(paste(" no significant hits for ", col, sep=""))
    write.table(ancom_res_df,file=paste("general_variable_check/ancom/genus/anc_res_",col,".csv",sep=""),sep=",",row.names = F,col.names = T)
    
  }
  
  
}
sink()

### Continuous data plotting #### ####
dir.create("continuous_data_plots")
#for full data, including continuous variables, read in below data:
all_samples <- read_excel("clinical_data/micobiomeECMM CAPA study 2021 04 27 v17short_edit_cut_filled_V2.xlsx", sheet = "samples")

#find column numbers of numeric variables 
numeric_vars <- c( which(colnames(all_samples)=="BMInew"),
which(colnames(all_samples)=="BALGMODI"),
which(colnames(all_samples)=="HGE") )
#set all excluding numeric variables as character
all_samples[,-numeric_vars] <- lapply(all_samples[,-numeric_vars] , as.character)
#set numeric variables to be numeric
all_samples[,numeric_vars] <- lapply(all_samples[,numeric_vars] , as.numeric)

#### recreate phyloseq object ####
#subset clinical data for samples available 
#alternatively, can subset later using subset_samples(physeq, <sample data column> =="Yes")
samples_df <- subset(all_samples, Sample %in% colnames(otu_mat) )
#Assign rownames for tax
samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample") 
#recreate taxa with full info (instead of qiime trimmed one needed for ancom)
TAX = tax_table(tax_mat2)
#recreate object with full clinical data
samples = sample_data(samples_df)
physeq = phyloseq(OTU, TAX, samples) 
#### normalise counts ####
total = median(sample_sums(physeq))
standf = function(x, t=total) round(t * (x / sum(x)))
physeqn = transform_sample_counts(physeq, standf)

#### filter low level species - this is filter for species occuring above 0.2% in ANY sample ####
physeq_abund <- filter_taxa(physeqn, function(x) sum(x > total*0.002) > 0, TRUE)

### Counts per sample ####
#total counts per sample  ##
sample_counts <- colSums(otu_table(physeq))

##Extract A fumigatus reads 
#subset taxa 
Af <- subset_taxa(physeq_abund, Species=="Aspergillus_fumigatus")
Af.counts <- colSums(otu_table(Af))

#filtered genus level data
physeq_abund_genera <- tax_glom(physeq_abund,taxrank="Genus_other")
#fix names to genus_other
phyloseq::taxa_names(physeq_abund_genera) <- phyloseq::tax_table(physeq_abund_genera)[, "Genus_other"]
####Extract Aspergillus reads from genus level data 
#subset taxa 
Asp <- subset_taxa(physeq_abund_genera, Genus_other=="Aspergillus")
Asp.counts <- as.data.frame(colSums(otu_table(Asp)))

####Extract Candida reads from genus level data 
#subset taxa 
Cand <- subset_taxa(physeq_abund_genera, Genus_other=="Candida")
Cand.counts <- as.data.frame(colSums(otu_table(Cand)))


#plot Candida against Aspergillus counts 
AC_counts <- cbind(Asp.counts,Cand.counts)
colnames(AC_counts) <- c("Aspergillus","Candida")
#remove sample with <1000 Aspergillus & Candida reads (combined)
AC_counts_clean <- AC_counts[which(rowSums(AC_counts)>1000),]

library(ggpubr)


### plot regression with full stats values title ###
#function for regression plotting 
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col="darkslategrey") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5))) +
    theme(plot.title = element_text(size=22))
  
}

#run function on Candida & Aspergillus count data
fit1 <- lm(Candida ~ Aspergillus, data = AC_counts_clean)

#set italic labels 
my_x_title <- expression(paste(italic("Aspergillus"), " counts"))
my_y_title <-  expression(paste(italic("Candida"), " counts"))
#plot regression
tiff("continuous_data_plots/Aspergillus_Candida_counts_scatter.tiff",width = 600, height = 600)
print(
  ggplotRegression(fit1)
  + geom_point(size=2) 
  + xlab(my_x_title)
  + ylab(my_y_title)
  + theme(panel.background = element_rect(fill = "white"), #change background to white
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 25, angle =0, hjust = 1, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
          legend.title=element_text(size=27),
          legend.text=element_text(size=25, face="italic"),
          legend.direction="vertical",)
)
dev.off()

### BAL GM values grouped by variable ####

#get balgm odi data
num <- which(colnames(samples_df)=="BALGMODI")
balgm <- sample_data(physeq)[,num]
balgm[2] <- rownames(balgm)
colnames(balgm)[2] <- "Sample"
#get sEOF data
num <- which(colnames(samples_df)=="SurvivalEOF")
sEOF <- sample_data(physeq)[,num]
#get capa data
num <- which(colnames(samples_df)=="CAPANoCAPA")
capastatus <- sample_data(physeq)[,num]
#merge bal gm odi, sEOF and capa data
sEOF_balgm <- cbind(sEOF,balgm,capastatus)

#remove unknown data points
sEOF_balgm_clean <- sEOF_balgm[!is.na(sEOF_balgm$BALGMODI),]

#calculate median and mean GM ODI values by sEOF 
sEOF_balgm_clean %>% group_by(SurvivalEOF) %>% summarise(Median= median(BALGMODI, na.rm=TRUE))
sEOF_balgm_clean %>% group_by(SurvivalEOF) %>% summarise(Mean= mean(BALGMODI, na.rm=TRUE))


#test if bal gm not normally distributed
shapiro.test(sEOF_balgm_clean$BALGMODI)
#data is not normally distributed


#kruskal test
kruskal.test(BALGMODI ~ SurvivalEOF, data = sEOF_balgm_clean )
#p val 0.7517
wilcox.test(BALGMODI ~ SurvivalEOF, data = sEOF_balgm_clean,  exact=FALSE )

svg("continuous_data_plots/sEOF_balgm_point_clean.svg",width = 5, height =7)
print(ggplot(sEOF_balgm_clean, aes(y=BALGMODI,x=SurvivalEOF))
      + geom_point()
      + geom_hline(yintercept=1, linetype="dashed", alpha=0.5)
      #+geom_boxplot()
      # + geom_violin()
      #+ scale_y_log10()
      #+ scale_color_manual(name="AFoutcome", values=c("1"="red","0"="black"),labels=c("CAPA","No CAPA"))
      + theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 20, angle =0, hjust = 1, vjust = 1, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text( color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="left")
)
dev.off()


svg("continuous_data_plots/sEOF_balgm_point_clean_capastatus.svg",width = 4.4, height = 5.2)
print(ggplot(sEOF_balgm_clean, aes(y=BALGMODI,x=SurvivalEOF,color=CAPANoCAPA))
      + geom_point(size=2.5)
      + geom_hline(yintercept=1, linetype="dashed", alpha=0.5)
      #+geom_boxplot()
      # + geom_violin()
      #+ scale_y_log10()
      #+ scale_color_manual(name="AFoutcome", values=c("1"="red","0"="black"),labels=c("CAPA","No CAPA"))
      + scale_color_manual(name="CAPANoCAPA", values=c("1"="red","0"="black"),labels=c("CAPA","No CAPA")) 
      + theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 20, angle =0, hjust = 1, vjust = 1, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="right") 
      + theme(panel.border = element_blank())
)
dev.off()


## grouped by systemic antifungals ##
#get sAF data
num <- which(colnames(samples_df)=="SystemicAFbeforeatICUAdmin")
sAF <- sample_data(physeq)[,num]

#merge bal gm odi, sAF and capa data
sAF_balgm <- cbind(sAF,balgm,capastatus)

#remove unknown data points
sAF_balgm_clean.tmp <- sAF_balgm[!is.na(sAF_balgm$BALGMODI),]
sAF_balgm_clean <- sAF_balgm_clean.tmp[!is.na(sAF_balgm_clean.tmp$SystemicAFbeforeatICUAdmin),]


#calculate median and mean GM ODI values by sEOF 
sAF_balgm_clean %>% group_by(SystemicAFbeforeatICUAdmin) %>% summarise(Median= median(BALGMODI, na.rm=TRUE))
sAF_balgm_clean %>% group_by(SystemicAFbeforeatICUAdmin) %>% summarise(Mean= mean(BALGMODI, na.rm=TRUE))


#kruskal test
kruskal.test(BALGMODI ~ SystemicAFbeforeatICUAdmin, data = sAF_balgm_clean )
#pval 0.055
wilcox.test(BALGMODI ~ SystemicAFbeforeatICUAdmin, data = sAF_balgm_clean, exact=FALSE )

tiff("continuous_data_plots/sAF_balgm_point_clean.tiff",width = 120, height = 300)
print(ggplot(sAF_balgm_clean, aes(y=BALGMODI,x=SystemicAFbeforeatICUAdmin))
 + geom_point()
#+geom_boxplot()
# + geom_violin()
#+ scale_y_log10()
+ geom_hline(yintercept=1, linetype="dashed", alpha=0.5)
#+ scale_color_manual(name="AFoutcome", values=c("1"="red","0"="black"),labels=c("CAPA","No CAPA"))
+ theme(panel.background = element_blank(),#change background to white
panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
axis.line = element_line(size = 1, colour = "black"),
axis.ticks = element_line(size = 2),
axis.text.x = element_text(color = "grey20", size = 20, angle =0, hjust = 1, vjust = 1, face = "plain"),
axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = .5, face = "plain"),
axis.title.x = element_text( color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0),
axis.title.y = element_text(color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5),
plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
legend.title=element_text(size=27),
legend.direction="vertical",
legend.position ="left")
)
dev.off()

svg("continuous_data_plots/sAF_balgm_point_clean_capastatus.svg",width = 4.4, height = 5.2)
print(ggplot(sAF_balgm_clean, aes(y=BALGMODI,x=SystemicAFbeforeatICUAdmin,color=CAPANoCAPA))
      + geom_point(size=2.5)
      + geom_hline(yintercept=1, linetype="dashed", alpha=0.5)
      #+geom_boxplot()
      # + geom_violin()
      #+ scale_y_log10()
      #+ scale_color_manual(name="AFoutcome", values=c("1"="red","0"="black"),labels=c("CAPA","No CAPA"))
      + scale_color_manual(name="CAPANoCAPA", values=c("1"="red","0"="black"),labels=c("CAPA","No CAPA")) 
      + theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 20, angle =0, hjust = 1, vjust = 1, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="right") 
      + theme(panel.border = element_blank())
)
dev.off()

### Counts & BAL GM values ####
# merge total counts and balgm data
balgm_sample_counts <- cbind(sample_counts,balgm)

#GMODI scatter - total counts
tiff("continuous_data_plots/BALGM_sample_counts_scatter.tiff",width = 1200, height = 1200)
print(ggplot(balgm_sample_counts, aes(x=BALGMODI,y=sample_counts))
      + geom_point(size=3) 
      + xlab("BAL GM ODI")
      # + scale_x_log10()
      #+ scale_y_log10()
      #+ scale_y_continuous(limits = c(0, NA))
      #expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
      + ylab("Total fungal counts")
      + geom_vline(xintercept = 1, linetype="dotted")
      + theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              legend.title=element_text(size=27),
              legend.text=element_text(size=25, face="italic"),
              legend.direction="vertical",)
)
dev.off()

#GMODI scatter - total counts
tiff("continuous_data_plots/BALGM_sample_counts_scatter_label.tiff",width = 1200, height = 1200)
print(ggplot(balgm_sample_counts, aes(x=BALGMODI,y=sample_counts,label=Sample))
      + geom_point(size=3) + geom_text(hjust=0, vjust=0,size=8)
      + xlab("BAL GM ODI")
      # + scale_x_log10()
      #+ scale_y_log10()
      #+ scale_y_continuous(limits = c(0, NA))
      #expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
      + ylab("Total fungal counts")
      + geom_vline(xintercept = 1, linetype="dotted")
      + theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              legend.title=element_text(size=27),
              legend.text=element_text(size=25, face="italic"),
              legend.direction="vertical",)
)
dev.off()

### Counts & qPCR HGE values ####
#get HGE data 
num2 <- which(colnames(samples_df)=="HGE")
hge <- sample_data(physeq)[,num2]
#merge total counts and hge (qpcr) data
hge_sample_counts <- cbind(sample_counts,hge)
#HGE scatter - total counts
tiff("continuous_data_plots/HGE_sample_counts_scatter.tiff",width = 1200, height = 1200)
print(ggplot(hge_sample_counts, aes(x=HGE,y=sample_counts))
      + geom_point(size=3)
      + xlab("HGE")
      + scale_x_log10()
      #+ scale_y_log10()
      #+ scale_y_continuous(limits = c(0, NA))
      #expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
      + ylab("Total fungal counts")
      + theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              legend.title=element_text(size=27),
              legend.text=element_text(size=25, face="italic"),
              legend.direction="vertical",)
)
dev.off()

### Af counts & qPCR HGE values ####
#merge Af counts and HGE data
hge_Af_counts <- cbind(Af.counts,hge)
#HGE scatter - total counts
tiff("continuous_data_plots/HGE_Af_counts_scatter.tiff",width = 1200, height = 1200)
print(ggplot(hge_Af_counts, aes(x=HGE,y=Af.counts))
      + geom_point(size=3)
      + xlab("HGE")
      + scale_x_log10()
      #+ scale_y_log10()
      #+ scale_y_continuous(limits = c(0, NA))
      #expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
      + ylab("A. fumigatus counts")
      + theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              legend.title=element_text(size=27),
              legend.text=element_text(size=25, face="italic"),
              legend.direction="vertical",)
)
dev.off()
#with regression 

#run function on qPCR & Aspergillus count data
fit2 <- lm(Af.counts ~ HGE, data = hge_Af_counts)

#set italic labels 
my_x_title <- expression(paste(italic("Aspergillus"), " counts"))
my_y_title <-  expression(paste("HGE"))
#plot regression
tiff("continuous_data_plots/HGE_Af_counts_regression.tiff",width = 800, height = 800)
print(
  ggplotRegression(fit2)
  + geom_point(size=2) 
  + xlab(my_x_title)
  + ylab(my_y_title)
  + scale_x_log10()
  + theme(panel.background = element_rect(fill = "white"), #change background to white
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
          legend.title=element_text(size=27),
          legend.text=element_text(size=25, face="italic"),
          legend.direction="vertical",)
)
dev.off()


## Af counts against HGE, coloured by CAPA status

#merge Af counts, hge and capa data
hge_Af_counts_capanocapa <- cbind(capastatus,hge_Af_counts)

#scatter
tiff("continuous_data_plots/HGE_Af_counts_capastatus_scatter.tiff",width = 1200, height = 1200)
print(ggplot(hge_Af_counts_capanocapa, aes(x=HGE,y=Af.counts,color=CAPANoCAPA))
      + geom_point(size=3)
      + xlab("HGE")
      + scale_x_log10()
      + scale_color_manual(name="CAPANoCAPA", values=c("1"="red","0"="black"),labels=c("CAPA","No CAPA"))
      #+ scale_y_log10()
      #+ scale_y_continuous(limits = c(0, NA))
      #expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
      + ylab("A. fumigatus counts")
      + theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              legend.title=element_text(size=27),
              legend.text=element_text(size=25, face="italic"),
              legend.direction="vertical",)
)
dev.off()

#Outliers in HGE data?
grubbs.test(hge$HGE, type = 10, opposite = FALSE, two.sided = FALSE)
#pval < 2.2e-16 for max outlier - value 10.86
grubbs.test(hge$HGE, type = 10, opposite = TRUE, two.sided = FALSE)
#no significant min outlier

### HGE by mortality and CAPA status ####
mortnum <- which(colnames(samples_df)=="SurvivalEOF")
mortality <- sample_data(physeq)[,mortnum]
mortality_hge <- cbind(cap,hge,mortality)
mortality_hge <- mortality_hge[!is.na(mortality_hge$SurvivalEOF),]

shapiro.test(mortality_hge$HGE)
pairwise.wilcox.test(mortality_hge$HGE, mortality_hge$SurvivalEOF, exact=FALSE)
#0.95

#calculate median HGE values by sEOF
HGE_median_byMort <- mortality_hge %>% group_by(SurvivalEOF) %>% summarise(Median= median(HGE, na.rm=TRUE))
HGE_mean_byMort <- mortality_hge %>% group_by(SurvivalEOF) %>% summarise(Mean= mean(HGE, na.rm=TRUE))

#mortality alone
tiff("continuous_data_plots/HGE_mort_box.tiff",width = 500, height = 900)
print(ggplot(mortality_hge, aes(y=HGE,x=SurvivalEOF))
      + geom_boxplot()+
      # + geom_violin()
      # + scale_y_log10()
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =0, hjust = 1, vjust = 1, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 35, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="right") 
)
dev.off()

tiff("continuous_data_plots/HGE_mort_capastatus_box.tiff",width = 500, height = 900)
print(ggplot(mortality_hge, aes(y=HGE,x=SurvivalEOF,color=CAPANoCAPA))
      + geom_boxplot()
      # + geom_violin()
      # + scale_y_log10()
      + scale_color_manual(name="CAPANoCAPA", values=c("1"="red","0"="black"),labels=c("CAPA","No CAPA")) +
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =0, hjust = 1, vjust = 1, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 35, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="right") 
)
dev.off()

tiff("continuous_data_plots/HGE_mort_capastatus_box_outlier_removed.tiff",width = 500, height = 900)
print(ggplot(mortality_hge, aes(y=HGE,x=SurvivalEOF,color=CAPANoCAPA))
        + geom_boxplot()
      # + geom_violin()
      + scale_y_continuous(limits=c(0, 2))
      # + scale_y_log10()
      + scale_color_manual(name="CAPANoCAPA", values=c("1"="red","0"="black"),labels=c("CAPA","No CAPA")) +
      theme(panel.background = element_blank(),#change background to white
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.line = element_line(size = 1, colour = "black"),
            axis.ticks = element_line(size = 2),
            axis.text.x = element_text(color = "grey20", size = 25, angle =0, hjust = 1, vjust = 1, face = "plain"),
            axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
            axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
            axis.title.y = element_text(face = "bold", color = "grey20", size = 35, angle = 90, hjust = .5, vjust = .5),
            plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
            legend.title=element_text(size=27),
            legend.direction="vertical",
            legend.position ="right") 
)
dev.off()
svg("continuous_data_plots/HGE_mort_capastatus_points.svg",width = 4.8, height = 5.2)
print(ggplot(mortality_hge, aes(y=HGE,x=SurvivalEOF,color=CAPANoCAPA))
      + geom_point(size=2.5)
      + scale_y_log10()
      + ylab(" HGE (log10 scale)")
      + geom_hline(yintercept = 0.1, linetype="dashed",color="grey")
      + scale_color_manual(name="CAPANoCAPA", values=c("1"="red","0"="black"),labels=c("CAPA","No CAPA")) +
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 20, angle =0, hjust = 1, vjust = 1, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="right") 
      + theme(panel.border = element_blank())
)
dev.off()


tiff("continuous_data_plots/HGE_mort_capastatus_points_outlier_removed.tiff",width = 500, height = 900)
print(ggplot(mortality_hge, aes(y=HGE,x=SurvivalEOF,color=CAPANoCAPA))
      + geom_point()
      + scale_y_continuous(limits=c(0, 2))
      + scale_color_manual(name="CAPANoCAPA", values=c("1"="red","0"="black"),labels=c("CAPA","No CAPA")) +
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =0, hjust = 1, vjust = 1, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 35, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="right") 
)
dev.off()

### HGE by systemic AF ####
safnum <- which(colnames(samples_df)=="SystemicAFbeforeatICUAdmin")
saf <- sample_data(physeq)[,safnum]
saf_hge <- cbind(cap,hge,saf)
saf_hge <- saf_hge[!is.na(saf_hge$SystemicAFbeforeatICUAdmin),]

shapiro.test(saf_hge$HGE)
pairwise.wilcox.test(saf_hge$HGE, saf_hge$SystemicAFbeforeatICUAdmin, exact=FALSE)
#0.089

#calculate median HGE values by sAF use 
HGE_median_bysAF <- saf_hge %>% group_by(SystemicAFbeforeatICUAdmin) %>% summarise(Median= median(HGE, na.rm=TRUE))
HGE_mean_bysAF <- saf_hge %>% group_by(SystemicAFbeforeatICUAdmin) %>% summarise(Mean= mean(HGE, na.rm=TRUE))


svg("continuous_data_plots/HGE_saf_capastatus_points.svg",width = 4.8, height = 5.2)
print(ggplot(saf_hge, aes(y=HGE,x=SystemicAFbeforeatICUAdmin,color=CAPANoCAPA))
      + geom_point(size=2.5)
      + scale_y_log10(limits=c(0.01,NA))
      + ylab(" HGE (log10 scale)")
      + geom_hline(yintercept = 0.1, linetype="dashed",color="grey")
      + scale_color_manual(name="CAPANoCAPA", values=c("1"="red","0"="black"),labels=c("CAPA","No CAPA")) +
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 20, angle =0, hjust = 1, vjust = 1, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="right") 
      + theme(panel.border = element_blank())
)
dev.off()

##HGE by systemic AF and survival EOF merged outcome ####
## merge data by SurvivalEOF status & plot ###
#create merged variable 
variable1 = as.character(get_variable(physeq_abund, "SurvivalEOF"))
variable2 = as.character(get_variable(physeq_abund, "SystemicAFbeforeatICUAdmin"))
sample_data(physeq_abund)$AFoutcome <- mapply(paste0, "S", variable1, "AF", variable2, 
                                              collapse = "_")
samples_df$AFoutcome <- mapply(paste0, "S", samples_df$SurvivalEOF, "AF", samples_df$SystemicAFbeforeatICUAdmin, 
                                              collapse = "_")

AFoutcome_num <- which(colnames(samples_df)=="AFoutcome")
outcome <- sample_data(physeq_abund)[,AFoutcome_num]
outcome_hge <- cbind(cap,hge,outcome)
outcome_hge <- outcome_hge[!is.na(outcome_hge$AFoutcome),]

outcome_hge$AFoutcome <- outcome_hge$AFoutcome %>% str_replace("S0AF0","No survival, no sAF") %>% str_replace("S0AF1","No survival + sAF")  %>% str_replace("S0AFNA","No survival, sAF unknown")  %>% str_replace("S1AF0","Survival, no sAF")  %>% str_replace("S1AF1","Survival + sAF")  %>% str_replace("S1AFNA","Survival, sAF unknown")  %>% str_replace("SNAAFNA","Unknown")
tiff("continuous_data_plots/HGE_sAFoutcome_box_outlier_removed.tiff",width = 800, height = 600)
print(ggplot(outcome_hge, aes(y=HGE,x=AFoutcome))
      + geom_boxplot()
      # + geom_violin()
      + scale_y_continuous(limits=c(0, 2))
      # + scale_y_log10()
      #+ scale_color_manual(name="AFoutcome", values=c("1"="red","0"="black"),labels=c("CAPA","No CAPA")) 
       + theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =45, hjust = 1, vjust = 1, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 35, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="none") 
)
dev.off()

tiff("continuous_data_plots/HGE_sAFoutcome_box.tiff",width = 800, height = 600)
print(ggplot(outcome_hge, aes(y=HGE,x=AFoutcome))
      + geom_boxplot()
      # + geom_violin()
      #+ scale_y_log10()
      #+ scale_color_manual(name="AFoutcome", values=c("1"="red","0"="black"),labels=c("CAPA","No CAPA")) 
      + theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =45, hjust = 1, vjust = 1, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 35, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="none") 
)
dev.off()
#remove unknown data points 
outcome_hge_clean <- outcome_hge[-grep("nknown",outcome_hge$AFoutcome),]

svg("continuous_data_plots/HGE_sAFoutcome_box_clean.svg",width = 4, height = 7.5)
print(ggplot(outcome_hge_clean, aes(y=HGE,x=AFoutcome))
      # + geom_point()
      +geom_boxplot()
      # + geom_violin()
      #+ scale_y_log10()
      #+ scale_color_manual(name="AFoutcome", values=c("1"="red","0"="black"),labels=c("CAPA","No CAPA")) 
      + theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =45, hjust = 1, vjust = 1, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text( color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="left",
              plot.margin = ggplot2::margin(1,1,1,1.75, "cm")) 
)
dev.off()

clean
shapiro.test(outcome_hge_clean$HGE)
#data is not normally distributed
#kruskal test
kruskal.test(HGE ~ AFoutcome, data = outcome_hge_clean )

#pval 0.1132

tiff("continuous_data_plots/HGE_sAFoutcome_capastatus_box_outlier_removed.tiff",width = 1200, height = 900)
print(ggplot(outcome_hge, aes(y=HGE,x=AFoutcome, color=CAPANoCAPA))
      + geom_boxplot()
      # + geom_violin()
      + scale_y_continuous(limits=c(0, 2))
      # + scale_y_log10()
      + scale_color_manual(name="AFoutcome", values=c("1"="red","0"="black"),labels=c("CAPA","No CAPA")) 
      + theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =45, hjust = 1, vjust = 1, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 35, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="right") 
)
dev.off()

tiff("continuous_data_plots/HGE_sAFoutcome_points.tiff",width = 800, height = 600)
print(ggplot(outcome_hge, aes(y=HGE,x=AFoutcome))
      + geom_point()
      # + geom_violin()
      + scale_y_log10()
      + ylab(" HGE (log10 scale)")
      #+ scale_color_manual(name="AFoutcome", values=c("1"="red","0"="black"),labels=c("CAPA","No CAPA")) 
      + theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =45, hjust = 1, vjust = 1, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 35, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="none") 
)
dev.off()

svg("continuous_data_plots/HGE_sAFoutcome_capastatus_points_clean.svg",width = 6, height = 6)
print(ggplot(outcome_hge_clean, aes(y=HGE,x=AFoutcome, color=CAPANoCAPA))
      + geom_point(size=2.5)
      + ylab(" HGE ")
      + scale_color_manual(name="CAPANoCAPA", values=c("1"="red","0"="black"),labels=c("CAPA","No CAPA")) +
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 20,angle =45, hjust = 1, vjust = 1, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="right",
              plot.margin = ggplot2::margin(1,1,1,1.75, "cm")) 
      + theme(panel.border = element_blank())
)
dev.off()

### corticos HGE & sample counts plotting ####
#total counts
num <- which(colnames(samples_df)=="Corticos")
cort <- sample_data(physeq)[,num]
cort[2] <- rownames(cort)
colnames(cort)[2] <- "Sample" 

cort_sample_counts <- cbind(sample_counts,cort)


tiff("continuous_data_plots/Corticos_sample_counts_box.tiff",width = 400, height = 400)
print(ggplot(cort_sample_counts, aes(x=Corticos,y=sample_counts))
      # + geom_point(size=3) 
      +geom_boxplot()
      + xlab("Corticosteroid use")
      # + scale_x_log10()
      #+ scale_y_log10()
      #+ scale_y_continuous(limits = c(0, NA))
      #expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
      + ylab("Total fungal counts")
      + theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              legend.title=element_text(size=27),
              legend.text=element_text(size=25, face="italic"),
              legend.direction="vertical",)
)
dev.off()



tiff("continuous_data_plots/Corticos_sample_counts_point.tiff",width = 400, height = 400)
print(ggplot(cort_sample_counts, aes(x=Corticos,y=sample_counts))
      + geom_point(size=2) 
      + xlab("Corticosteroid use")
      # + scale_x_log10()
      #+ scale_y_log10()
      #+ scale_y_continuous(limits = c(0, NA))
      #expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
      + ylab("Total fungal counts")
      + theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              legend.title=element_text(size=27),
              legend.text=element_text(size=25, face="italic"),
              legend.direction="vertical",)
)
dev.off()

### HGE by Corticos ####
num2 <- which(colnames(samples_df)=="HGE")
hge <- sample_data(physeq)[,num2]
hge_corticos <- cbind(cort,hge)

#calculate median and mean HGE values by Corticosteroid use 
HGE_median_byCorticos <- hge_corticos %>% group_by(Corticos) %>% summarise(Median= median(HGE, na.rm=TRUE))
HGE_mean_byCorticos <- hge_corticos %>% group_by(Corticos) %>% summarise(Mean= mean(HGE, na.rm=TRUE))

#GMODI scatter - total counts
tiff("continuous_data_plots/Corticos_HGE_box.tiff",width = 400, height = 400)
print(ggplot(hge_corticos, aes(x=Corticos,y=HGE))
      # + geom_point(size=3) 
      +geom_boxplot()
      + xlab("Corticosteroid use")
      # + scale_x_log10()
      #+ scale_y_log10()
      #+ scale_y_continuous(limits = c(0, NA))
      #expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
      + ylab("HGE")
      + theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              legend.title=element_text(size=27),
              legend.text=element_text(size=25, face="italic"),
              legend.direction="vertical",)
)
dev.off()



###  HGE by CAPA ####
#total counts
num <- which(colnames(samples_df)=="CAPANoCAPA")
cap <- sample_data(physeq)[,num]
cap[2] <- rownames(cap)
colnames(cap)[2] <- "Sample" 

cap_sample_counts <- cbind(sample_counts,cap)


tiff("continuous_data_plots/CAPANoCAPA_sample_counts_box.tiff",width = 300, height = 400)
print(ggplot(cap_sample_counts, aes(x=CAPANoCAPA,y=sample_counts))
      # + geom_point(size=3) 
      +geom_boxplot()
      + xlab("CAPA")
      # + scale_x_log10()
      #+ scale_y_log10()
      #+ scale_y_continuous(limits = c(0, NA))
      #expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
      + ylab("Total fungal counts")
      + theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              legend.title=element_text(size=27),
              legend.text=element_text(size=25, face="italic"),
              legend.direction="vertical",)
)
dev.off()



tiff("continuous_data_plots/CAPANoCAPA_sample_counts_point.tiff",width = 300, height = 400)
print(ggplot(cap_sample_counts, aes(x=CAPANoCAPA,y=sample_counts))
      + geom_point(size=2) 
      + xlab("CAPA")
      # + scale_x_log10()
      #+ scale_y_log10()
      #+ scale_y_continuous(limits = c(0, NA))
      #expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
      + ylab("Total fungal counts")
      + theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              legend.title=element_text(size=27),
              legend.text=element_text(size=25, face="italic"),
              legend.direction="vertical",)
)
dev.off()
shapiro.test(cap_sample_counts$sample_counts)
pairwise.wilcox.test(cap_sample_counts$sample_counts, cap_sample_counts$CAPANoCAPA, exact=FALSE)

num2 <- which(colnames(samples_df)=="HGE")
hge <- sample_data(physeq)[,num2]
hge_capastatus <- cbind(cap,hge)


pairwise.wilcox.test(hge_capastatus$HGE, hge_capastatus$CAPANoCAPA, exact=FALSE)
#calculate median HGE values by CAPA status 
HGE_median_byCAPA <- hge_capastatus %>% group_by(CAPANoCAPA) %>% summarise(Median= median(HGE, na.rm=TRUE))
HGE_mean_byCAPA <- hge_capastatus %>% group_by(CAPANoCAPA) %>% summarise(Mean= mean(HGE, na.rm=TRUE))


#HGE by capa boxplot
svg("continuous_data_plots/CAPANoCAPA_HGE_box.svg",width = 4, height = 6)
print(ggplot(hge_capastatus, aes(x=CAPANoCAPA,y=HGE, color=CAPANoCAPA))
      # + geom_point(size=3) 
      
      + geom_boxplot(outlier.shape = NA)
      + geom_jitter(height = 0, width = .2,size=2.5)
      + xlab("CAPA")
      + scale_color_manual(values=c("#10639B","#808080")) 
      # + scale_x_log10()
      #+ scale_y_log10()
      #+ scale_y_continuous(limits = c(0, NA))
      #expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
      + ylab("HGE")
      + theme(panel.background = element_blank(),#change background to white
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.line = element_line(size = 1, colour = "black"),
            axis.ticks = element_line(size = 2),
            axis.text.x = element_blank(),
            axis.text.y = element_text(color = "grey20", size = 27, angle = 0, hjust = 1, vjust = .5, face = "plain"),
            axis.title.x = element_text(face = "bold", color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0),
            axis.title.y = element_text(face = "bold", color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5),
            plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
            legend.title=element_text(size=10),
            legend.text=element_text(size=20),
            legend.direction="vertical",
            legend.position ="none",
            plot.margin = ggplot2::margin(1,1,1,1.75, "cm")) +
        theme(panel.border = element_blank())
)
dev.off()

#HGE by capa boxplot, without displaying outlier (trim Y axis)
svg("continuous_data_plots/CAPANoCAPA_HGE_box_outlier_not_shown.svg",width = 4, height = 6)
print(ggplot(hge_capastatus, aes(x=CAPANoCAPA,y=HGE, color=CAPANoCAPA))
      # + geom_point(size=3) 
      
      + geom_boxplot(outlier.shape = NA)
      + geom_jitter(height = 0, width = .2,size=2.5)
      + xlab("CAPA")
      + scale_color_manual(values=c("#10639B","#808080")) 
      # + scale_x_log10()
      #+ scale_y_log10()
      + scale_y_continuous(limits = c(0, 2))
      #expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
      + ylab("HGE")
      + theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_blank(),
              axis.text.y = element_text(color = "grey20", size = 27, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=10),
              legend.text=element_text(size=20),
              legend.direction="vertical",
              legend.position ="none",
              plot.margin = ggplot2::margin(1,1,1,1.75, "cm")) +
        theme(panel.border = element_blank())
)
dev.off()

#HGE by capa boxplot, with y axis break between outlier and other data
library(ggbreak)
svg("continuous_data_plots/CAPANoCAPA_HGE_box_axis_break.svg",width = 4, height = 7)
print(ggplot(hge_capastatus, aes(x=CAPANoCAPA,y=HGE, color=CAPANoCAPA))
      # + geom_point(size=3) 
      + scale_y_continuous(limits = c(0, 12), breaks = c(1,2))
      + scale_y_break(c(2, 10), ticklabels=c(10, 11),space=.7) 
      + geom_boxplot(outlier.shape = NA)
      + geom_jitter(height = 0, width = .2,size=2.5)
      + xlab("CAPA")
      + scale_color_manual(values=c("#10639B","#808080")) 
      + ylab("HGE")
      + theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_blank(),
              axis.text.y = element_text(color = "grey20", size = 27, angle = 0, hjust = .5, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=10),
              legend.text=element_text(size=20),
              legend.direction="vertical",
              legend.position ="none",
              plot.margin = ggplot2::margin(1,1,1,1.75, "cm")) +
        theme(panel.border = element_blank())
)
dev.off()

#get full GM, qPCR, CAPA and Af count table
num4 <- which(colnames(samples_df)=="BALGM1ODI")
GMresult <- sample_data(physeq)[,num4]

GM_hge_af_capa <- cbind(hge_Af_counts_capanocapa,balgm,GMresult)
write.table(GM_hge_af_capa,file="full_GM_Af_qPCR_info.csv",sep=",",row.names = F)

### qPCR HGE & BAL GM values ####
#total counts

num2 <- which(colnames(samples_df)=="HGE")
hge <- sample_data(physeq)[,num2]

balgm_hge <- cbind(hge,balgm)
#GMODI scatter - HGE
tiff("continuous_data_plots/BALGM_HGE_scatter.tiff",width = 1200, height = 1200)
print(ggplot(balgm_hge, aes(x=BALGMODI,y=HGE,colour=cut(HGE,c(-Inf,0,20))))
      + geom_point(size=3)
      + xlab("BAL GM ODI")
      # + scale_x_log10()
      # + scale_y_log10()
      + scale_color_manual(name="HGE", values=c("(-Inf,0]"="red","(0,20]"="black"),labels=c("qPCR negative","qPCR positive"))
      + scale_x_continuous(limits = c(0, NA),
                           expand = expansion(mult = c(0, 0.1)))
      + scale_y_continuous(limits = c(0, NA),
      expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
      + ylab("HGE")
      + geom_vline(xintercept = 1, linetype="dotted")
      + theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              legend.title=element_text(size=27),
              legend.text=element_text(size=25, face="italic"),
              legend.direction="vertical",)
)
dev.off()

#GMODI scatter - HGE
tiff("continuous_data_plots/BALGM_HGE_scatter_label.tiff",width = 1200, height = 1200)
print(ggplot(balgm_hge, aes(x=BALGMODI,y=HGE,label=Sample,colour=cut(HGE,c(-Inf,0,20))))
      + geom_point(size=3) + geom_text(hjust=0, vjust=0,size=8)
      + xlab("BAL GM ODI")
      # + scale_x_log10()
      # + scale_y_log10()
      + scale_color_manual(name="HGE", values=c("(-Inf,0]"="red","(0,20]"="black"),labels=c("qPCR negative","qPCR positive"))
      + scale_x_continuous(limits = c(0, NA),
                           expand = expansion(mult = c(0, 0.1)))
      + scale_y_continuous(limits = c(0, NA),
                           expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
      + ylab("HGE")
      + geom_vline(xintercept = 1, linetype="dotted")
      + theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              legend.title=element_text(size=27),
              legend.text=element_text(size=25, face="italic"),
              legend.direction="vertical",)
)
dev.off()

### Species loop ##
#take top 20 rank sum taxa
rank_20 <- sort(taxa_sums(physeq_abund), TRUE)[1:20]/nsamples(physeq_abund)
for (i in seq(1,length(names(rank_20)))){
  spe <- names(rank_20)[i]
  #subset taxa 
  physeq_subset <- subset_taxa(physeq_abund, Taxa==spe)
  #Extract species reads 
  subset.counts <- colSums(otu_table(physeq_subset))
  #Merge with GM data 
  balgm_species_counts <- cbind(subset.counts,balgm)
  
  #GMODI & species counts scatter  
  tiff(paste("continuous_data_plots/BALGM_", spe,"_counts_scatter.tiff",sep=""),width = 1200, height = 1200)
  print(ggplot(balgm_species_counts, aes(x=BALGMODI,y=subset.counts))
        + geom_point(size=3)
        + xlab("BAL GM ODI")
        # + scale_x_log10()
        # + scale_y_log10()
        #+ scale_y_continuous(limits = c(0, NA))
        #expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
        + ylab(paste(spe," counts",sep=""))
        + theme(panel.background = element_rect(fill = "white"), #change background to white
                axis.line = element_line(size = 1, colour = "black"),
                axis.ticks = element_line(size = 2),
                axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
                axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
                axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
                axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
                legend.title=element_text(size=27),
                legend.text=element_text(size=25, face="italic"),
                legend.direction="vertical",)
  )
  #scatter with labelled samples 
  tiff(paste("continuous_data_plots/BALGM_", spe,"_counts_scatter_label.tiff",sep=""),width = 1200, height = 1200)
  print(ggplot(balgm_species_counts, aes(x=BALGMODI,y=subset.counts,label=Sample))
        + geom_point(size=3)  + geom_text(hjust=0, vjust=0,size=8)
        + xlab("BAL GM ODI")
        # + scale_x_log10()
        # + scale_y_log10()
        #+ scale_y_continuous(limits = c(0, NA))
        #expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
        + ylab(paste(spe," counts",sep=""))
        + theme(panel.background = element_rect(fill = "white"), #change background to white
                axis.line = element_line(size = 1, colour = "black"),
                axis.ticks = element_line(size = 2),
                axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
                axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
                axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
                axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
                legend.title=element_text(size=27),
                legend.text=element_text(size=25, face="italic"),
                legend.direction="vertical",)
  )
  dev.off()
}

#Aspergillus 
balgm_Asp_counts <- cbind(Asp.counts,balgm)
colnames(balgm_Asp_counts)[1] <- "Asp.counts"
#GMODI scatter
tiff("continuous_data_plots/BALGM_Asp_counts_scatter.tiff",width = 1200, height = 1200)
print(ggplot(balgm_Asp_counts, aes(x=BALGMODI,y=Asp.counts))
      + geom_point(size=3)
      + xlab("BAL GM ODI")
      # + scale_x_log10()
      # + scale_y_log10()
      #+ scale_y_continuous(limits = c(0, NA))
      #expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
      + ylab("Aspergillus counts")
      + theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              legend.title=element_text(size=27),
              legend.text=element_text(size=25, face="italic"),
              legend.direction="vertical",)
)
dev.off()

### richness & BAL GM values ####
balgm_shannon <- cbind(Shann,balgm[,-2])
#GMODI scatter - total counts
tiff("continuous_data_plots/BALGM_shannon_scatter.tiff",width = 1200, height = 1200)
print(ggplot(balgm_shannon, aes(x=BALGMODI,y=Shannon))
      + geom_point(size=3)
      + xlab("BAL GM ODI")
      # + scale_x_log10()
      #+ scale_y_log10()
      #+ scale_y_continuous(limits = c(0, NA))
      #expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
      + ylab("Shannon")
      + theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              legend.title=element_text(size=27),
              legend.text=element_text(size=25, face="italic"),
              legend.direction="vertical",)
)
dev.off()

##shannon and duration ICU
which(colnames(samples_df)=="DurationICUDays")
var <- sample_data(physeq)[,50]
var_shannon <- cbind(Shann,var)
#GMODI scatter - total counts
tiff("continuous_data_plots/DurationICUDays_shannon_scatter.tiff",width = 1200, height = 1200)
print(ggplot(var_shannon, aes(x=DurationICUDays,y=Shannon))
      + geom_point(size=3)
      + xlab("DurationICUDays")
      + geom_smooth(method = "lm", se = TRUE)
      + stat_cor(aes(label=..rr.label..), label.x=30)
      # + scale_x_log10()
      #+ scale_y_log10()
      #+ scale_y_continuous(limits = c(0, NA))
      #expand = expansion(mult = c(0, 0.1))) #makes Y axis start from origin and expands so that 1 is not cut off
      + ylab("Shannon")
      + theme(panel.background = element_rect(fill = "white"), #change background to white
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 25, angle =90, hjust = .5, vjust = .5, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 27, angle = 90, hjust = .5, vjust = .5),
              legend.title=element_text(size=27),
              legend.text=element_text(size=25, face="italic"),
              legend.direction="vertical",)
)
dev.off()


##sample batch ######################
###permanova test on bray ordination ###
wunifrac_dist = phyloseq::distance(physeq_abund, method="bray", weighted=F)
ordination = ordinate(physeq_abund, method="PCoA", distance=wunifrac_dist)

tiff("Batch_ordination.tiff",width = 600, height = 600)
plot_ordination(physeq_abund, ordination, color="Batch") + theme(aspect.ratio=1)  +
  geom_point(size=3) + 
  theme(panel.background = element_blank(), #change background to white
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.line = element_line(size = 1, colour = "black"),
            axis.ticks = element_line(size = 2),
            axis.text.x = element_text(color = "grey20", size = 30, angle =0, hjust = .5, vjust = .5, face = "plain"),
            axis.text.y = element_text(color = "grey20", size = 30, angle = 0, hjust = 1, vjust = .5, face = "plain"),
            axis.title.x = element_text(face = "bold", color = "grey20", size = 32, angle = 0, hjust = .5, vjust = 0),
            axis.title.y = element_text(face = "bold", color = "grey20", size = 32, angle = 90, hjust = .5, vjust = 0),
            legend.title=element_text(size=27),
            legend.text=element_text(size=25, face="italic"),
            legend.direction="vertical",
            legend.position = "bottom",
            strip.text.x = element_text(size = 32, colour = "black"))
dev.off()
library(vegan)
## tests###
##remove NA samples for plotting
col4 <- which(colnames(samples_df)=="Batch")
var_values <- sample_data(physeq_abund)[[col4]]
clean2 <- prune_samples(!is.na(var_values), physeq_abund)

metadata <- as(sample_data(clean2), "data.frame")

res.permanova <- adonis2(phyloseq::distance(clean2, method = "bray") ~ Batch,
       data = metadata)
write.table(res.permanova,file="Batch_permanova_results.csv",sep=",", row.names = F)

### Species abundance boxplots ####
## function to melt phyloseq object into long format to make ggplot happy
x11()
phyloseq::psmelt(physeq_abund_genera) %>%
  ggplot(data = ., aes(x = SurvivalEOF, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")


phyloseq::psmelt(physeq_filt_genera) %>%
  ggplot(data = ., aes(x = SurvivalEOF, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")

## boxplot selected taxa abundances
dir.create("specific_comparisons/abundance_boxplots")
library(dplyr)
#plot abundances for top 5 species by variable - raw data -with colours, CAPA
for (i in seq(1,(length(names(rank_20))/4))){
  spe <- names(rank_20)[i]
  #species
  p <- phyloseq::psmelt(physeq) %>% dplyr::filter(Species == spe ) %>% dplyr::filter(CAPANoCAPA != "<NA>") %>%
    ggplot(data = ., aes(x = CAPANoCAPA, y = Abundance, color=CAPANoCAPA)) +
    geom_boxplot(outlier.shape = NA, width=0.7) +
    geom_jitter(height = 0, width = .2,size=2.5) +
    labs(y = "Abundance\n") +
    scale_color_manual(values=c("#10639B","#808080"))  +
    #(values=c("#7570B3","#E6AB02"))  purple/yellow
    #(values=c("#7570B3","#00846b")) purple/green
    #(values=c("#10639B","#808080")) blue/grey
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    theme(panel.background = element_blank(),#change background to white
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = "grey20", size = 27, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=10),
          legend.direction="vertical",
          legend.position ="none",
          plot.margin = ggplot2::margin(1,1,1,1.75, "cm")) +
  theme(panel.border = element_blank())
  
  svg(paste("specific_comparisons/abundance_boxplots/",spe,"_","CAPANoCAPA","_abundance.svg",sep=""),width = 5, height = 6)
  print(p)
  dev.off() 
}
#check if normal
Af_only_data <- phyloseq::psmelt(physeq) %>% dplyr::filter(Species == "Aspergillus_fumigatus" )
shapiro.test(Af_only_data$Abundance)

#find median and mean A.f abundance by capa status
phyloseq::psmelt(physeq) %>% dplyr::filter(Species == "Aspergillus_fumigatus" ) %>%
  dplyr::filter(CAPANoCAPA != "<NA>") %>% group_by(CAPANoCAPA) %>% summarise(Med=median(Abundance,na.rm=TRUE))

phyloseq::psmelt(physeq) %>% dplyr::filter(Species == "Aspergillus_fumigatus" ) %>%
  dplyr::filter(CAPANoCAPA != "<NA>") %>% group_by(CAPANoCAPA) %>% summarise(Mean=mean(Abundance,na.rm=TRUE))

#plot abundances for top 5 species by variable - raw data -with colours, Corticos
for (i in seq(1,(length(names(rank_20))/4))){
  spe <- names(rank_20)[i]
  #species
  p <- phyloseq::psmelt(physeq) %>% dplyr::filter(Species == spe ) %>% dplyr::filter(Corticos != "<NA>") %>%
    ggplot(data = ., aes(x = Corticos, y = Abundance, color=Corticos)) +
    geom_boxplot(outlier.shape = NA, width=0.7) +
    geom_jitter(height = 0, width = .2,size=2.5) +
    labs(y = "Abundance\n") +
    scale_color_manual(values=c("#7570B3","#E6AB02"))   +
    #(values=c("#7570B3","#E6AB02"))  purple/yellow
    #(values=c("#7570B3","#00846b")) purple/green
    #(values=c("#10639B","#808080")) blue/grey
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    theme(panel.background = element_blank(),#change background to white
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = "grey20", size = 27, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=10),
          legend.direction="vertical",
          legend.position ="none",
          plot.margin = ggplot2::margin(1,1,1,1.75, "cm")) +
    theme(panel.border = element_blank())
  
  svg(paste("specific_comparisons/abundance_boxplots/",spe,"_","Corticos","_abundance.svg",sep=""),width = 5, height = 6)
  print(p)
  dev.off() 
}

#find median and mean A.f abundance by Corticos
phyloseq::psmelt(physeq) %>% dplyr::filter(Species == "Aspergillus_fumigatus" ) %>%
  dplyr::filter(Corticos != "<NA>") %>% group_by(Corticos) %>% summarise(Med=median(Abundance,na.rm=TRUE))

phyloseq::psmelt(physeq) %>% dplyr::filter(Species == "Aspergillus_fumigatus" ) %>%
  dplyr::filter(Corticos != "<NA>") %>% group_by(Corticos) %>% summarise(Mean=mean(Abundance,na.rm=TRUE))

#plot abundances for top 5 species by variable - raw data  - no colours - sEOF
for (i in seq(1,(length(names(rank_20))/4))){
  spe <- names(rank_20)[i]
  #species
  p <- phyloseq::psmelt(physeq) %>% dplyr::filter(Species == spe ) %>% dplyr::filter(SurvivalEOF != "<NA>") %>%
    ggplot(data = ., aes(x = SurvivalEOF, y = Abundance)) +
    geom_boxplot(outlier.shape = NA, width=0.7) +
    geom_jitter(height = 0, width = .2,size=2.5) +
    labs(y = "Abundance\n") +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    theme(panel.background = element_blank(),#change background to white
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = "grey20", size = 27, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=10),
          legend.direction="vertical",
          legend.position ="none",
          plot.margin = ggplot2::margin(1,1,1,1.75, "cm")) +
    theme(panel.border = element_blank())
  
  svg(paste("specific_comparisons/abundance_boxplots/",spe,"_","SurvivalEOF","_abundance_uncoloured.svg",sep=""),width = 4, height = 7)
  print(p)
  dev.off() 
}


#find median and mean A.f abundance by SurvivalEOF
phyloseq::psmelt(physeq) %>% dplyr::filter(Species == "Aspergillus_fumigatus" ) %>%
  dplyr::filter(SurvivalEOF != "<NA>") %>% group_by(SurvivalEOF) %>% 
  summarise(Med=median(Abundance,na.rm=TRUE))

phyloseq::psmelt(physeq) %>% dplyr::filter(Species == "Aspergillus_fumigatus" ) %>%
  dplyr::filter(SurvivalEOF != "<NA>") %>% group_by(SurvivalEOF) %>%
  summarise(Mean=mean(Abundance,na.rm=TRUE))

#plot abundances for top 5 species by variable - raw data  - no colours - sAF
for (i in seq(1,(length(names(rank_20))/4))){
  spe <- names(rank_20)[i]
  #species
  p <- phyloseq::psmelt(physeq) %>% dplyr::filter(Species == spe ) %>% dplyr::filter(SystemicAFbeforeatICUAdmin != "<NA>") %>%
    ggplot(data = ., aes(x = SystemicAFbeforeatICUAdmin, y = Abundance)) +
    geom_boxplot(outlier.shape = NA, width=0.7) +
    geom_jitter(height = 0, width = .2,size=2.5) +
    labs(y = "Abundance\n") +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    theme(panel.background = element_blank(),#change background to white
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = "grey20", size = 27, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=10),
          legend.direction="vertical",
          legend.position ="none",
          plot.margin = ggplot2::margin(1,1,1,1.75, "cm")) +
    theme(panel.border = element_blank())
  
  svg(paste("specific_comparisons/abundance_boxplots/",spe,"_","SystemicAFbeforeatICUAdmin","_abundance_uncoloured.svg",sep=""),width = 4, height = 6.5)
  print(p)
  dev.off() 
}

Af_only_data_sAF <- phyloseq::psmelt(physeq) %>% dplyr::filter(Species == "Aspergillus_fumigatus" ) %>%
  dplyr::filter(SystemicAFbeforeatICUAdmin != "<NA>") 
shapiro.test(Af_only_data_sAF$Abundance)

#find median and mean A.f abundance by sAF
phyloseq::psmelt(physeq) %>% dplyr::filter(Species == "Aspergillus_fumigatus" ) %>%
  dplyr::filter(SystemicAFbeforeatICUAdmin != "<NA>") %>% group_by(SystemicAFbeforeatICUAdmin) %>% 
  summarise(Med=median(Abundance,na.rm=TRUE))

phyloseq::psmelt(physeq) %>% dplyr::filter(Species == "Aspergillus_fumigatus" ) %>%
  dplyr::filter(SystemicAFbeforeatICUAdmin != "<NA>") %>% group_by(SystemicAFbeforeatICUAdmin) %>%
  summarise(Mean=mean(Abundance,na.rm=TRUE))

#plot abundances for top 5 species by variable - physeq_abund data 
for (i in seq(1,(length(names(rank_20))/4))){
  spe <- names(rank_20)[i]
  #species
  p <- phyloseq::psmelt(physeq_abund) %>% dplyr::filter(Species == spe ) %>%
    ggplot(data = ., aes(x = SurvivalEOF, y = Abundance)) +
    #geom_violin()+
    geom_boxplot() +
    geom_jitter(height = 0, width = .2) +
    labs(y = "Abundance\n") +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    theme(panel.background = element_blank(),#change background to white
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 25, angle =0, hjust = 1, vjust = 1, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 35, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=27),
          legend.direction="vertical",
          legend.position ="right") 
  
  tiff(paste("specific_comparisons/abundance_boxplots/",spe,"_","SurvivalEOF","norm_abundance.tiff",sep=""),width = 400, height = 520)
  #width = 340, height = 400)
  print(p)
  dev.off() 
}

#genus
p <- phyloseq::psmelt(physeq_genera) %>% dplyr::filter(Genus_other =="Candida" ) %>%
  ggplot(data = ., aes(x = SystemicAFbeforeatICUAdmin, y = Abundance)) +
  geom_boxplot() +
  geom_jitter(height = 0, width = .2) +
  labs(y = "Abundance\n") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(panel.background = element_blank(),#change background to white
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 1, colour = "black"),
        axis.ticks = element_line(size = 2),
        axis.text.x = element_text(color = "grey20", size = 25, angle =0, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = .5, face = "plain"),
        axis.title.x = element_text(face = "bold", color = "grey20", size = 27, angle = 0, hjust = .5, vjust = 0),
        axis.title.y = element_text(face = "bold", color = "grey20", size = 35, angle = 90, hjust = .5, vjust = .5),
        plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
        legend.title=element_text(size=27),
        legend.direction="vertical",
        legend.position ="right") 

tiff("specific_comparisons/abundance_boxplots/SystemicAFbeforeatICUAdmin_Candida_abundance.tiff",width = 400, height = 510)
print(p)
dev.off() 


#wilcox of raw abundances 
Af_data <- phyloseq::psmelt(physeq_abund) %>% dplyr::filter(Species == "Aspergillus_fumigatus" )
shapiro.test(Af_data$Abundance)
pairwise.wilcox.test(Af_data$Abundance, Af_data$CAPANoCAPA, exact=FALSE)
pairwise.wilcox.test(Af_data$Abundance, Af_data$Corticos, exact=FALSE)
pairwise.wilcox.test(Af_data$Abundance, Af_data$SystemicAFbeforeatICUAdmin, exact=FALSE)
pairwise.wilcox.test(Af_data$Abundance, Af_data$SurvivalEOF, exact=FALSE)
kruskal.test(Abundance ~ AFoutcome, data =Af_data)

### Taxa prevalence ####
#### species prevalence - full data ####
##(fractional) ##
#define 5% of median counts as cutoff 
perc1 <- round(median(sample_sums(physeq_abund)) * 0.05)
##positive group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf = apply(X = otu_table(physeq_abund),
                 MARGIN = ifelse(taxa_are_rows(physeq_abund), yes = 1, no = 2),
                 FUN = function(x){sum(x > perc1)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(physeq_abund),
                      tax_table(physeq_abund))

prevdf <- prevdf[,c("Species","Prevalence","TotalAbundance")]

#keep taxa found in more than 10% of samples only ##
keep <- prevdf$Species[which(prevdf$Prevalence > 0.1)]

clean_prevdf <- prevdf[which(prevdf$Species %in% keep),]

#prep taxa names for plotting
clean_prevdf$Species <- str_replace(clean_prevdf$Species,"_"," ")
#set width of plot based on number of taxa
num_taxa <- length(unique(clean_prevdf$Species))
plot.width <- (num_taxa * 50) + 100

##ggplot coloured barplot 
tiff("Species_prevalence.tiff",width = plot.width, height = 500)
print(ggplot(clean_prevdf,aes(x=reorder(Species, -Prevalence),y=Prevalence))  +
        geom_bar(stat="identity", width =0.5,position = "dodge") + theme_bw() +
        xlab("Species") +
        scale_y_continuous(limits = c(0, 1),
                           expand = expansion(mult = c(0, 0.1))) +
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
              axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="bottom",
              plot.margin = margin(1,1,1,1.75, "cm")) 
)
dev.off()


#### genera prevalence - full data ####
##(fractional) ##
#define 1% of median counts as cutoff 
perc1 <- round(median(sample_sums(physeq_abund)) * 0.01)
##positive group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf = apply(X = otu_table(physeq_abund),
               MARGIN = ifelse(taxa_are_rows(physeq_abund), yes = 1, no = 2),
               FUN = function(x){sum(x > perc1)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(physeq_abund),
                    tax_table(physeq_abund))

prevdf <- prevdf[,c("Genus_other","Prevalence","TotalAbundance")]

#keep taxa found in more than 10% of samples only ##
keep <- prevdf$Genus_other[which(prevdf$Prevalence > 0.1)]

clean_prevdf <- prevdf[which(prevdf$Genus_other %in% keep),]

#prep taxa names for plotting
clean_prevdf$Genus_other <- str_replace(clean_prevdf$Genus_other,"_"," ")
#set width of plot based on number of taxa
num_taxa <- length(unique(clean_prevdf$Genus_other))
plot.width <- (num_taxa * 50) + 100

##ggplot coloured barplot 
tiff("Genera_prevalence.tiff",width = plot.width, height = 500)
print(ggplot(clean_prevdf,aes(x=reorder(Genus_other, -Prevalence),y=Prevalence))  +
        geom_bar(stat="identity", width =0.5,position = "dodge") + theme_bw() +
        xlab("Genus_other") +
        scale_y_continuous(limits = c(0, 1),
                           expand = expansion(mult = c(0, 0.1))) +
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
              axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="bottom",
              plot.margin = margin(1,1,1,1.75, "cm")) 
)
dev.off()

#### species prevalence - SurvivalEOF ####
##(fractional) ##
#define 5% of median counts as cutoff 
perc1 <- round(median(sample_sums(noSurvivalEOF_physeq_sp)) * 0.05)
##positive group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_1 = apply(X = otu_table(SurvivalEOF_physeq_sp),
                 MARGIN = ifelse(taxa_are_rows(SurvivalEOF_physeq_sp), yes = 1, no = 2),
                 FUN = function(x){sum(x > perc1)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_1 = data.frame(Prevalence = prevdf_1,
                      TotalAbundance = taxa_sums(SurvivalEOF_physeq_sp),
                      tax_table(SurvivalEOF_physeq_sp))
#add descriptor for data source
prevdf_1$type <- "pos"
##negative group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_0 = apply(X = otu_table(noSurvivalEOF_physeq_sp),
                 MARGIN = ifelse(taxa_are_rows(noSurvivalEOF_physeq_sp), yes = 1, no = 2),
                 FUN = function(x){sum(x > perc1)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_0 = data.frame(Prevalence = prevdf_0,
                      TotalAbundance = taxa_sums(noSurvivalEOF_physeq_sp),
                      tax_table(noSurvivalEOF_physeq_sp))
#add descriptor for data source
prevdf_0$type <- "neg"

prevdf <- rbind(prevdf_0,prevdf_1)

prevdf <- prevdf[,c("Species","Prevalence","TotalAbundance","type")]

#keep taxa found in more than 10% of samples only (in either group) ##
keep <- prevdf$Species[which(prevdf$Prevalence > 0.1)]

clean_prevdf <- prevdf[which(prevdf$Species %in% keep),]
#prep taxa names for plotting
clean_prevdf$Species <- str_replace(clean_prevdf$Species,"_"," ")
#set width of plot based on number of taxa
num_taxa <- length(unique(clean_prevdf$Species))
plot.width <- (num_taxa * 50) + 100

##ggplot coloured barplot 
tiff("specific_comparisons/sEOF/SurvivalEOF_Species_prevalence.tiff",width = plot.width, height = 500)
print(ggplot(clean_prevdf,aes(x=reorder(Species,-Prevalence),y=Prevalence, fill=type))  +
    geom_bar(stat="identity", width =0.5,position = "dodge") + theme_bw() +
    xlab("Species") +
    scale_y_continuous(limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(values=c("#10639B","#808080")) +
    theme(panel.background = element_blank(),#change background to white
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(size = 1, colour = "black"),
          axis.ticks = element_line(size = 2),
          axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
          axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
          axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
          axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
          plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
          legend.title=element_text(size=27),
          legend.direction="vertical",
          legend.position ="bottom",
          plot.margin = margin(1,1,1,1.75, "cm")) 
)
dev.off()
plot.height <- (num_taxa ) + 1.5
##ggplot coloured barplot - vertical
svg("specific_comparisons/sEOF/SurvivalEOF_Species_prevalence_vert.svg",width = 8, height = plot.height)
print(ggplot(clean_prevdf,aes(x=reorder(Species,Prevalence),y=Prevalence, fill=type))  +
        geom_bar(stat="identity", width =0.5,position = "dodge") + theme_bw() +
        xlab("Species") +
        coord_flip() +
        scale_y_continuous(limits = c(0, 1),
                           expand = expansion(mult = c(0, 0.1))) +
        scale_fill_manual(values=c("#10639B","#808080")) +
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.y = element_text(color = "grey20", size = 29, angle =0, hjust = 1, vjust = .5, face = "italic"),
              axis.text.x = element_text(color = "grey20", size = 29, angle = 40, hjust = 1, vjust = 1, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="bottom",
              plot.margin = margin(1,1,1,1.75, "cm")) +
        theme(panel.border = element_blank())
)
dev.off()
## prevalence of dominance (fractional) ##
#define 50% of median counts as cutoff for dominating a sample
perc50 <- round(median(sample_sums(SurvivalEOF_physeq_sp)) * 0.5)
##positive group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_1.dom = apply(X = otu_table(SurvivalEOF_physeq_sp),
               MARGIN = ifelse(taxa_are_rows(SurvivalEOF_physeq_sp), yes = 1, no = 2),
               FUN = function(x){sum(x > perc50)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_1.dom = data.frame(Prevalence = prevdf_1.dom,
                    TotalAbundance = taxa_sums(SurvivalEOF_physeq_sp),
                    tax_table(SurvivalEOF_physeq_sp))
#add descriptor for data source
prevdf_1.dom$type <- "pos"

##negative group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_0.dom = apply(X = otu_table(noSurvivalEOF_physeq_sp),
               MARGIN = ifelse(taxa_are_rows(noSurvivalEOF_physeq_sp), yes = 1, no = 2),
               FUN = function(x){sum(x > perc50)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_0.dom = data.frame(Prevalence = prevdf_0.dom,
                      TotalAbundance = taxa_sums(noSurvivalEOF_physeq_sp),
                      tax_table(noSurvivalEOF_physeq_sp))

#add descriptor for data source
prevdf_0.dom$type <- "neg"

prevdf.dom <- rbind(prevdf_0.dom,prevdf_1.dom)

prevdf.dom <- prevdf.dom[,c("Species","Prevalence","TotalAbundance","type")]

#keep taxa found in more than 10% of samples only (in either group) #
keep <- prevdf.dom$Species[which(prevdf.dom$Prevalence > 0.1)]

clean_prevdf.dom <- prevdf.dom[which(prevdf.dom$Species %in% keep),]
#prep taxa names for plotting
clean_prevdf.dom$Species <- str_replace(clean_prevdf.dom$Species,"_"," ")
#set width of plot based on number of taxa
num_taxa <- length(unique(clean_prevdf.dom$Species))
plot.width <- (num_taxa * 50) + 150
##ggplot coloured barplot 
tiff("specific_comparisons/sEOF/SurvivalEOF_Species_prevalence_dominance.tiff",width = plot.width, height = 500)
print(ggplot(clean_prevdf.dom,aes(x=reorder(Species,-Prevalence),y=Prevalence, fill=type))  +
        geom_bar(stat="identity", width =0.5,position = "dodge") + theme_bw() +
        xlab("Species") +
        scale_y_continuous(limits = c(0, 1),
                           expand = expansion(mult = c(0, 0.1))) +
        scale_fill_manual(values=c("#10639B","#808080")) +
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
              axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="bottom",
              plot.margin = margin(1,1,1,1.75, "cm")) 
)
dev.off()

#### species prevalence - SystemicAFbeforeatICUAdmin ####
##(fractional) ##
#define 5% of median counts as cutoff 
perc1 <- round(median(sample_sums(noSystemicAFbeforeatICUAdmin_physeq_sp)) * 0.05)
##positive group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_1 = apply(X = otu_table(SystemicAFbeforeatICUAdmin_physeq_sp),
                 MARGIN = ifelse(taxa_are_rows(SystemicAFbeforeatICUAdmin_physeq_sp), yes = 1, no = 2),
                 FUN = function(x){sum(x > perc1)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_1 = data.frame(Prevalence = prevdf_1,
                      TotalAbundance = taxa_sums(SystemicAFbeforeatICUAdmin_physeq_sp),
                      tax_table(SystemicAFbeforeatICUAdmin_physeq_sp))
#add descriptor for data source
prevdf_1$type <- "pos"
##negative group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_0 = apply(X = otu_table(noSystemicAFbeforeatICUAdmin_physeq_sp),
                 MARGIN = ifelse(taxa_are_rows(noSystemicAFbeforeatICUAdmin_physeq_sp), yes = 1, no = 2),
                 FUN = function(x){sum(x > perc1)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_0 = data.frame(Prevalence = prevdf_0,
                      TotalAbundance = taxa_sums(noSystemicAFbeforeatICUAdmin_physeq_sp),
                      tax_table(noSystemicAFbeforeatICUAdmin_physeq_sp))
#add descriptor for data source
prevdf_0$type <- "neg"

prevdf <- rbind(prevdf_0,prevdf_1)

prevdf <- prevdf[,c("Species","Prevalence","TotalAbundance","type")]

#keep taxa found in more than 10% of samples only (in either group) ##
keep <- prevdf$Species[which(prevdf$Prevalence > 0.1)]

clean_prevdf <- prevdf[which(prevdf$Species %in% keep),]
#prep taxa names for plotting
clean_prevdf$Species <- str_replace(clean_prevdf$Species,"_"," ")
#set width of plot based on number of taxa
num_taxa <- length(unique(clean_prevdf$Species))
plot.width <- (num_taxa * 50) + 100
##ggplot coloured barplot 
tiff("specific_comparisons/sAF/SystemicAFbeforeatICUAdmin_Species_prevalence.tiff",width = plot.width, height = 600)
print(ggplot(clean_prevdf,aes(x=reorder(Species,-Prevalence),y=Prevalence, fill=type))  +
        geom_bar(stat="identity", width =0.5,position = "dodge") + theme_bw() +
        xlab("Species") +
        scale_y_continuous(limits = c(0, 1),
                           expand = expansion(mult = c(0, 0.1))) +
        scale_fill_manual(values=c("#10639B","#808080")) +
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
              axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="bottom",
              plot.margin = margin(1,1,1,1.75, "cm")) 
)
dev.off()
plot.height <- (num_taxa ) +.2
##ggplot coloured barplot - vertical
svg("specific_comparisons/sAF/SystemicAFbeforeatICUAdmin_Species_prevalence_vert.svg",width = 10, height = plot.height)
print(ggplot(clean_prevdf,aes(x=reorder(Species,Prevalence),y=Prevalence, fill=type))  +
        geom_bar(stat="identity", width =0.5,position = "dodge") + theme_bw() +
        xlab("Species") +
        coord_flip() +
        scale_y_continuous(limits = c(0, 1),
                           expand = expansion(mult = c(0, 0.1))) +
        scale_fill_manual(values=c("#10639B","#808080")) +
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.y = element_text(color = "grey20", size = 32, angle =0, hjust = 1, vjust = .5, face = "italic"),
              axis.text.x = element_text(color = "grey20", size = 32, angle = 40, hjust = 1, vjust = 1, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="bottom",
              plot.margin = margin(1,1,1,1.75, "cm")) +
        theme(panel.border = element_blank())
)
dev.off()
### prevalence of dominance (fractional) ##
#define 50% of median counts as cutoff for dominating a sample
perc50 <- round(median(sample_sums(SystemicAFbeforeatICUAdmin_physeq_sp)) * 0.5)
##positive group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_1.dom = apply(X = otu_table(SystemicAFbeforeatICUAdmin_physeq_sp),
                     MARGIN = ifelse(taxa_are_rows(SystemicAFbeforeatICUAdmin_physeq_sp), yes = 1, no = 2),
                     FUN = function(x){sum(x > perc50)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_1.dom = data.frame(Prevalence = prevdf_1.dom,
                          TotalAbundance = taxa_sums(SystemicAFbeforeatICUAdmin_physeq_sp),
                          tax_table(SystemicAFbeforeatICUAdmin_physeq_sp))
#add descriptor for data source
prevdf_1.dom$type <- "pos"

##negative group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_0.dom = apply(X = otu_table(noSystemicAFbeforeatICUAdmin_physeq_sp),
                     MARGIN = ifelse(taxa_are_rows(noSystemicAFbeforeatICUAdmin_physeq_sp), yes = 1, no = 2),
                     FUN = function(x){sum(x > perc50)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_0.dom = data.frame(Prevalence = prevdf_0.dom,
                          TotalAbundance = taxa_sums(noSystemicAFbeforeatICUAdmin_physeq_sp),
                          tax_table(noSystemicAFbeforeatICUAdmin_physeq_sp))

#add descriptor for data source
prevdf_0.dom$type <- "neg"

prevdf.dom <- rbind(prevdf_0.dom,prevdf_1.dom)

prevdf.dom <- prevdf.dom[,c("Species","Prevalence","TotalAbundance","type")]

#keep taxa found in more than 10% of samples only (in either group) #
keep <- prevdf.dom$Species[which(prevdf.dom$Prevalence > 0.1)]

clean_prevdf.dom <- prevdf.dom[which(prevdf.dom$Species %in% keep),]
#prep taxa names for plotting
clean_prevdf.dom$Species <- str_replace(clean_prevdf.dom$Species,"_"," ")
#set width of plot based on number of taxa
num_taxa <- length(unique(clean_prevdf.dom$Species))
plot.width <- (num_taxa * 50) + 150
##ggplot coloured barplot 
tiff("specific_comparisons/sAF/SystemicAFbeforeatICUAdmin_Species_prevalence_dominance.tiff",width = plot.width, height = 600)
print(ggplot(clean_prevdf.dom,aes(x=reorder(Species,-Prevalence),y=Prevalence, fill=type))  +
        geom_bar(stat="identity", width =0.5,position = "dodge") + theme_bw() +
        xlab("Species") +
        scale_y_continuous(limits = c(0, 1),
                           expand = expansion(mult = c(0, 0.1))) +
        scale_fill_manual(values=c("#10639B","#808080")) +
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
              axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="bottom",
              plot.margin = margin(1,1,1,1.75, "cm")) 
)
dev.off()

#### species prevalence - CAPANoCAPA ####

##(fractional) ##
#define 5% of median counts as cutoff 
perc1 <- round(median(sample_sums(noCAPA_physeq_sp)) * 0.05)
##positive group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_1 = apply(X = otu_table(CAPA_physeq_sp),
                 MARGIN = ifelse(taxa_are_rows(CAPA_physeq_sp), yes = 1, no = 2),
                 FUN = function(x){sum(x > perc1)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_1 = data.frame(Prevalence = prevdf_1,
                      TotalAbundance = taxa_sums(CAPA_physeq_sp),
                      tax_table(CAPA_physeq_sp))
#add descriptor for data source
prevdf_1$type <- "pos"
##negative group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_0 = apply(X = otu_table(noCAPA_physeq_sp),
                 MARGIN = ifelse(taxa_are_rows(noCAPA_physeq_sp), yes = 1, no = 2),
                 FUN = function(x){sum(x > perc1)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_0 = data.frame(Prevalence = prevdf_0,
                      TotalAbundance = taxa_sums(noCAPA_physeq_sp),
                      tax_table(noCAPA_physeq_sp))
#add descriptor for data source
prevdf_0$type <- "neg"

prevdf <- rbind(prevdf_0,prevdf_1)

prevdf <- prevdf[,c("Species","Prevalence","TotalAbundance","type")]

#keep taxa found in more than 10% of samples only (in either group) ##
keep <- prevdf$Species[which(prevdf$Prevalence > 0.1)]

clean_prevdf <- prevdf[which(prevdf$Species %in% keep),]
#prep taxa names for plotting
clean_prevdf$Species <- str_replace(clean_prevdf$Species,"_"," ")
#set width of plot based on number of taxa
num_taxa <- length(unique(clean_prevdf$Species))
plot.width <- (num_taxa * 50) + 100
##ggplot coloured barplot 
tiff("specific_comparisons/CAPA/CAPANoCAPA_Species_prevalence.tiff",width = plot.width, height = 600)
print(ggplot(clean_prevdf,aes(x=reorder(Species,-Prevalence),y=Prevalence, fill=type))  +
        geom_bar(stat="identity", width =0.5,position = "dodge") + theme_bw() +
        xlab("Species") +
        scale_y_continuous(limits = c(0, 1),
                           expand = expansion(mult = c(0, 0.1))) +
        scale_fill_manual(values=c("#10639B","#808080")) +
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
              axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="bottom",
              plot.margin = margin(1,1,1,1.75, "cm")) 
)
dev.off()

### prevalence of dominance (fractional) ##
#define 50% of median counts as cutoff for dominating a sample
perc50 <- round(median(sample_sums(CAPA_physeq_sp)) * 0.5)
##positive group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_1.dom = apply(X = otu_table(CAPA_physeq_sp),
                     MARGIN = ifelse(taxa_are_rows(CAPA_physeq_sp), yes = 1, no = 2),
                     FUN = function(x){sum(x > perc50)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_1.dom = data.frame(Prevalence = prevdf_1.dom,
                          TotalAbundance = taxa_sums(CAPA_physeq_sp),
                          tax_table(CAPA_physeq_sp))
#add descriptor for data source
prevdf_1.dom$type <- "pos"

##negative group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_0.dom = apply(X = otu_table(noCAPA_physeq_sp),
                     MARGIN = ifelse(taxa_are_rows(noCAPA_physeq_sp), yes = 1, no = 2),
                     FUN = function(x){sum(x > perc50)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_0.dom = data.frame(Prevalence = prevdf_0.dom,
                          TotalAbundance = taxa_sums(noCAPA_physeq_sp),
                          tax_table(noCAPA_physeq_sp))

#add descriptor for data source
prevdf_0.dom$type <- "neg"

prevdf.dom <- rbind(prevdf_0.dom,prevdf_1.dom)

prevdf.dom <- prevdf.dom[,c("Species","Prevalence","TotalAbundance","type")]

#keep taxa found in more than 10% of samples only (in either group) #
keep <- prevdf.dom$Species[which(prevdf.dom$Prevalence > 0.1)]

clean_prevdf.dom <- prevdf.dom[which(prevdf.dom$Species %in% keep),]
#prep taxa names for plotting
clean_prevdf.dom$Species <- str_replace(clean_prevdf.dom$Species,"_"," ")
#set width of plot based on number of taxa
num_taxa <- length(unique(clean_prevdf.dom$Species))
plot.width <- (num_taxa * 50) + 150
##ggplot coloured barplot 
tiff("specific_comparisons/CAPA/CAPANoCAPA_Species_prevalence_dominance.tiff",width = plot.width, height = 550)
print(ggplot(clean_prevdf.dom,aes(x=reorder(Species,-Prevalence),y=Prevalence, fill=type))  +
        geom_bar(stat="identity", width =0.5,position = "dodge") + theme_bw() +
        xlab("Species") +
        scale_y_continuous(limits = c(0, 1),
                           expand = expansion(mult = c(0, 0.1))) +
        scale_fill_manual(values=c("#10639B","#808080")) +
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
              axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="bottom",
              plot.margin = margin(1,1,1,1.75, "cm")) 
)
dev.off()

#### genera prevalence - SurvivalEOF ####
##(fractional) ##
#define 1% of median counts as cutoff 
perc1 <- round(median(sample_sums(noSurvivalEOF_physeq)) * 0.01)
##positive group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_1 = apply(X = otu_table(SurvivalEOF_physeq),
                 MARGIN = ifelse(taxa_are_rows(SurvivalEOF_physeq), yes = 1, no = 2),
                 FUN = function(x){sum(x > perc1)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_1 = data.frame(Prevalence = prevdf_1,
                      TotalAbundance = taxa_sums(SurvivalEOF_physeq),
                      tax_table(SurvivalEOF_physeq))
#add descriptor for data source
prevdf_1$type <- "pos"
##negative group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_0 = apply(X = otu_table(noSurvivalEOF_physeq),
                 MARGIN = ifelse(taxa_are_rows(noSurvivalEOF_physeq), yes = 1, no = 2),
                 FUN = function(x){sum(x > perc1)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_0 = data.frame(Prevalence = prevdf_0,
                      TotalAbundance = taxa_sums(noSurvivalEOF_physeq),
                      tax_table(noSurvivalEOF_physeq))
#add descriptor for data source
prevdf_0$type <- "neg"

prevdf <- rbind(prevdf_0,prevdf_1)

prevdf <- prevdf[,c("Genus_other","Prevalence","TotalAbundance","type")]

#keep taxa found in more than 10% of samples only (in either group) ##
keep <- prevdf$Genus_other[which(prevdf$Prevalence > 0.1)]

clean_prevdf <- prevdf[which(prevdf$Genus_other %in% keep),]
#prep taxa names for plotting
clean_prevdf$Species <- str_replace(clean_prevdf$Species,"_"," ")
#set width of plot based on number of taxa
num_taxa <- length(unique(clean_prevdf$Genus_other))
plot.width <- (num_taxa * 50) + 100

##ggplot coloured barplot 
tiff("specific_comparisons/sEOF/SurvivalEOF_Genera_prevalence.tiff",width = plot.width, height = 500)
print(ggplot(clean_prevdf,aes(x=Genus_other,y=Prevalence, fill=type))  +
        geom_bar(stat="identity", width =0.5,position = "dodge") + theme_bw() +
        xlab("Species") +
        scale_y_continuous(limits = c(0, 1),
                           expand = expansion(mult = c(0, 0.1))) +
        scale_fill_manual(values=c("#10639B","#808080")) +
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
              axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="bottom",
              plot.margin = margin(1,1,1,1.75, "cm")) 
)
dev.off()

### prevalence of dominance (fractional) ##
#define 50% of median counts as cutoff for dominating a sample
perc50 <- round(median(sample_sums(SurvivalEOF_physeq)) * 0.5)
##positive group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_1.dom = apply(X = otu_table(SurvivalEOF_physeq),
                     MARGIN = ifelse(taxa_are_rows(SurvivalEOF_physeq), yes = 1, no = 2),
                     FUN = function(x){sum(x > perc50)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_1.dom = data.frame(Prevalence = prevdf_1.dom,
                          TotalAbundance = taxa_sums(SurvivalEOF_physeq),
                          tax_table(SurvivalEOF_physeq))
#add descriptor for data source
prevdf_1.dom$type <- "pos"

##negative group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_0.dom = apply(X = otu_table(noSurvivalEOF_physeq),
                     MARGIN = ifelse(taxa_are_rows(noSurvivalEOF_physeq), yes = 1, no = 2),
                     FUN = function(x){sum(x > perc50)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_0.dom = data.frame(Prevalence = prevdf_0.dom,
                          TotalAbundance = taxa_sums(noSurvivalEOF_physeq),
                          tax_table(noSurvivalEOF_physeq))

#add descriptor for data source
prevdf_0.dom$type <- "neg"

prevdf.dom <- rbind(prevdf_0.dom,prevdf_1.dom)

prevdf.dom <- prevdf.dom[,c("Genus_other","Prevalence","TotalAbundance","type")]

#keep taxa found in more than 10% of samples only (in either group) #
keep <- prevdf.dom$Genus_other[which(prevdf.dom$Prevalence > 0.1)]

clean_prevdf.dom <- prevdf.dom[which(prevdf.dom$Genus_other %in% keep),]

#set width of plot based on number of taxa
num_taxa <- length(unique(clean_prevdf.dom$Genus_other))
plot.width <- (num_taxa * 50) + 150
##ggplot coloured barplot 
tiff("specific_comparisons/sEOF/SurvivalEOF_Genera_prevalence_dominance.tiff",width = plot.width, height = 500)
print(ggplot(clean_prevdf.dom,aes(x=Genus_other,y=Prevalence, fill=type))  +
        geom_bar(stat="identity", width =0.5,position = "dodge") + theme_bw() +
        xlab("Species") +
        scale_y_continuous(limits = c(0, 1),
                           expand = expansion(mult = c(0, 0.1))) +
        scale_fill_manual(values=c("#10639B","#808080")) +
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
              axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="bottom",
              plot.margin = margin(1,1,1,1.75, "cm")) 
)
dev.off()

#### genera prevalence - SystemicAFbeforeatICUAdmin ####
##(fractional) ##
#define 1% of median counts as cutoff 
perc1 <- round(median(sample_sums(noSystemicAFbeforeatICUAdmin_physeq)) * 0.01)
##positive group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_1 = apply(X = otu_table(SystemicAFbeforeatICUAdmin_physeq),
                 MARGIN = ifelse(taxa_are_rows(SystemicAFbeforeatICUAdmin_physeq), yes = 1, no = 2),
                 FUN = function(x){sum(x > perc1)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_1 = data.frame(Prevalence = prevdf_1,
                      TotalAbundance = taxa_sums(SystemicAFbeforeatICUAdmin_physeq),
                      tax_table(SystemicAFbeforeatICUAdmin_physeq))
#add descriptor for data source
prevdf_1$type <- "pos"
##negative group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_0 = apply(X = otu_table(noSystemicAFbeforeatICUAdmin_physeq),
                 MARGIN = ifelse(taxa_are_rows(noSystemicAFbeforeatICUAdmin_physeq), yes = 1, no = 2),
                 FUN = function(x){sum(x > perc1)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_0 = data.frame(Prevalence = prevdf_0,
                      TotalAbundance = taxa_sums(noSystemicAFbeforeatICUAdmin_physeq),
                      tax_table(noSystemicAFbeforeatICUAdmin_physeq))
#add descriptor for data source
prevdf_0$type <- "neg"

prevdf <- rbind(prevdf_0,prevdf_1)

prevdf <- prevdf[,c("Genus_other","Prevalence","TotalAbundance","type")]

#keep taxa found in more than 10% of samples only (in either group) ##
keep <- prevdf$Genus_other[which(prevdf$Prevalence > 0.1)]

clean_prevdf <- prevdf[which(prevdf$Genus_other %in% keep),]
#prep taxa names for plotting
clean_prevdf$Species <- str_replace(clean_prevdf$Species,"_"," ")
#set width of plot based on number of taxa
num_taxa <- length(unique(clean_prevdf$Genus_other))
plot.width <- (num_taxa * 50) + 100
##ggplot coloured barplot 
tiff("specific_comparisons/sAF/SystemicAFbeforeatICUAdmin_Genera_prevalence.tiff",width = plot.width, height = 500)
print(ggplot(clean_prevdf,aes(x=Genus_other,y=Prevalence, fill=type))  +
        geom_bar(stat="identity", width =0.5,position = "dodge") + theme_bw() +
        xlab("Species") +
        scale_y_continuous(limits = c(0, 1),
                           expand = expansion(mult = c(0, 0.1))) +
        scale_fill_manual(values=c("#10639B","#808080")) +
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
              axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="bottom",
              plot.margin = margin(1,1,1,1.75, "cm")) 
)
dev.off()

### prevalence of dominance (fractional) ##
#define 50% of median counts as cutoff for dominating a sample
perc50 <- round(median(sample_sums(SystemicAFbeforeatICUAdmin_physeq)) * 0.5)
##positive group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_1.dom = apply(X = otu_table(SystemicAFbeforeatICUAdmin_physeq),
                     MARGIN = ifelse(taxa_are_rows(SystemicAFbeforeatICUAdmin_physeq), yes = 1, no = 2),
                     FUN = function(x){sum(x > perc50)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_1.dom = data.frame(Prevalence = prevdf_1.dom,
                          TotalAbundance = taxa_sums(SystemicAFbeforeatICUAdmin_physeq),
                          tax_table(SystemicAFbeforeatICUAdmin_physeq))
#add descriptor for data source
prevdf_1.dom$type <- "pos"

##negative group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_0.dom = apply(X = otu_table(noSystemicAFbeforeatICUAdmin_physeq),
                     MARGIN = ifelse(taxa_are_rows(noSystemicAFbeforeatICUAdmin_physeq), yes = 1, no = 2),
                     FUN = function(x){sum(x > perc50)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_0.dom = data.frame(Prevalence = prevdf_0.dom,
                          TotalAbundance = taxa_sums(noSystemicAFbeforeatICUAdmin_physeq),
                          tax_table(noSystemicAFbeforeatICUAdmin_physeq))

#add descriptor for data source
prevdf_0.dom$type <- "neg"

prevdf.dom <- rbind(prevdf_0.dom,prevdf_1.dom)

prevdf.dom <- prevdf.dom[,c("Genus_other","Prevalence","TotalAbundance","type")]

#keep taxa found in more than 10% of samples only (in either group) #
keep <- prevdf.dom$Genus_other[which(prevdf.dom$Prevalence > 0.1)]

clean_prevdf.dom <- prevdf.dom[which(prevdf.dom$Genus_other %in% keep),]
#set width of plot based on number of taxa
num_taxa <- length(unique(clean_prevdf.dom$Genus_other))
plot.width <- (num_taxa * 50) + 150
##ggplot coloured barplot 
tiff("specific_comparisons/sAF/SystemicAFbeforeatICUAdmin_Genera_prevalence_dominance.tiff",width = plot.width, height = 500)
print(ggplot(clean_prevdf.dom,aes(x=Genus_other,y=Prevalence, fill=type))  +
        geom_bar(stat="identity", width =0.5,position = "dodge") + theme_bw() +
        xlab("Species") +
        scale_y_continuous(limits = c(0, 1),
                           expand = expansion(mult = c(0, 0.1))) +
        scale_fill_manual(values=c("#10639B","#808080")) +
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
              axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="bottom",
              plot.margin = margin(1,1,1,1.75, "cm")) 
)
dev.off()

#### genera prevalence - CAPANoCAPA ####

##(fractional) ##
#define 1% of median counts as cutoff 
perc1 <- round(median(sample_sums(noCAPA_physeq)) * 0.01)
##positive group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_1 = apply(X = otu_table(CAPA_physeq),
                 MARGIN = ifelse(taxa_are_rows(CAPA_physeq), yes = 1, no = 2),
                 FUN = function(x){sum(x > perc1)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_1 = data.frame(Prevalence = prevdf_1,
                      TotalAbundance = taxa_sums(CAPA_physeq),
                      tax_table(CAPA_physeq))
#add descriptor for data source
prevdf_1$type <- "pos"
##negative group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_0 = apply(X = otu_table(noCAPA_physeq),
                 MARGIN = ifelse(taxa_are_rows(noCAPA_physeq), yes = 1, no = 2),
                 FUN = function(x){sum(x > perc1)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_0 = data.frame(Prevalence = prevdf_0,
                      TotalAbundance = taxa_sums(noCAPA_physeq),
                      tax_table(noCAPA_physeq))
#add descriptor for data source
prevdf_0$type <- "neg"

prevdf <- rbind(prevdf_0,prevdf_1)

prevdf <- prevdf[,c("Genus_other","Prevalence","TotalAbundance","type")]

#keep taxa found in more than 10% of samples only (in either group) ##
keep <- prevdf$Genus_other[which(prevdf$Prevalence > 0.1)]

clean_prevdf <- prevdf[which(prevdf$Genus_other %in% keep),]
#prep taxa names for plotting
clean_prevdf$Species <- str_replace(clean_prevdf$Species,"_"," ")
#set width of plot based on number of taxa
num_taxa <- length(unique(clean_prevdf$Genus_other))
plot.width <- (num_taxa * 50) + 100
##ggplot coloured barplot 
tiff("specific_comparisons/CAPA/CAPANoCAPA_Genera_prevalence.tiff",width = 500, height = 500)
print(ggplot(clean_prevdf,aes(x=Genus_other,y=Prevalence, fill=type))  +
        geom_bar(stat="identity", width =0.5,position = "dodge") + theme_bw() +
        xlab("Species") +
        scale_y_continuous(limits = c(0, 1),
                           expand = expansion(mult = c(0, 0.1))) +
        scale_fill_manual(values=c("#10639B","#808080")) +
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
              axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="bottom",
              plot.margin = margin(1,1,1,1.75, "cm")) 
)
dev.off()

### prevalence of dominance (fractional) ##
#define 50% of median counts as cutoff for dominating a sample
perc50 <- round(median(sample_sums(CAPA_physeq)) * 0.5)
##positive group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_1.dom = apply(X = otu_table(CAPA_physeq),
                     MARGIN = ifelse(taxa_are_rows(CAPA_physeq), yes = 1, no = 2),
                     FUN = function(x){sum(x > perc50)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_1.dom = data.frame(Prevalence = prevdf_1.dom,
                          TotalAbundance = taxa_sums(CAPA_physeq),
                          tax_table(CAPA_physeq))
#add descriptor for data source
prevdf_1.dom$type <- "pos"

##negative group
# Compute prevalence of each feature (as fraction), store as data.frame
prevdf_0.dom = apply(X = otu_table(noCAPA_physeq),
                     MARGIN = ifelse(taxa_are_rows(noCAPA_physeq), yes = 1, no = 2),
                     FUN = function(x){sum(x > perc50)/length(x)})

# Add taxonomy and total read counts to this data.frame
prevdf_0.dom = data.frame(Prevalence = prevdf_0.dom,
                          TotalAbundance = taxa_sums(noCAPA_physeq),
                          tax_table(noCAPA_physeq))

#add descriptor for data source
prevdf_0.dom$type <- "neg"

prevdf.dom <- rbind(prevdf_0.dom,prevdf_1.dom)

prevdf.dom <- prevdf.dom[,c("Genus_other","Prevalence","TotalAbundance","type")]

#keep taxa found in more than 10% of samples only (in either group) #
keep <- prevdf.dom$Genus_other[which(prevdf.dom$Prevalence > 0.1)]

clean_prevdf.dom <- prevdf.dom[which(prevdf.dom$Genus_other %in% keep),]

#set width of plot based on number of taxa
num_taxa <- length(unique(clean_prevdf.dom$Genus_other))
plot.width <- (num_taxa * 50) + 150
##ggplot coloured barplot 
tiff("specific_comparisons/CAPA/CAPANoCAPA_Genera_prevalence_dominance.tiff",width = plot.width, height = 500)
print(ggplot(clean_prevdf.dom,aes(x=Genus_other,y=Prevalence, fill=type))  +
        geom_bar(stat="identity", width =0.5,position = "dodge") + theme_bw() +
        xlab("Species") +
        scale_y_continuous(limits = c(0, 1),
                           expand = expansion(mult = c(0, 0.1))) +
        scale_fill_manual(values=c("#10639B","#808080")) +
        theme(panel.background = element_blank(),#change background to white
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.line = element_line(size = 1, colour = "black"),
              axis.ticks = element_line(size = 2),
              axis.text.x = element_text(color = "grey20", size = 23, angle =60, hjust = 1, vjust = 1, face = "italic"),
              axis.text.y = element_text(color = "grey20", size = 23, angle = 0, hjust = 1, vjust = .5, face = "plain"),
              axis.title.x = element_text(face = "bold", color = "grey20", size = 30, angle = 0, hjust = .5, vjust = 0),
              axis.title.y = element_text(face = "bold", color = "grey20", size = 30, angle = 90, hjust = .5, vjust = .5),
              plot.title = element_text(size = 29, face = "bold",hjust=0.5, vjust=-.5),
              legend.title=element_text(size=27),
              legend.direction="vertical",
              legend.position ="bottom",
              plot.margin = margin(1,1,1,1.75, "cm")) 
)
dev.off()

### Save data ####
save.image(file='COVID_paper_phyloseq_data.RData')





