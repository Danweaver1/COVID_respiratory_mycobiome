# COVID_respiratory_mycobiome
Respiratory ITS1 mycobiome study of critically ill COVID-19 patients. Sequencing data is available at NCBI SRA accession: PRJNA905224.

## Summary of main scripts
**ITS_PE_pipeline_v2_COVID_final.sh:** Performs read QC and mapping to UNITE database.

**Full_CAPA_ITS_count_analysis_final.sh:** Performs run-run count scaling and filters species based on two filters:
1. species read cutoff (SRC) - eg. at least 10 reads for species to remain
2. species percent cutoff (SPC) - eg. at least 0.2% of total species to remain

**Fetch_data.sh:** Retreives scaled and filtered count data using a key file which indicates samples passing the initial qPCR-based QC of extraction batches.

**create_CAPA_phyloseq_sample_filtered_otu_final.R:** Converts the count data into a matrix compatible with Phyloseq package. Also performs further sample filtering based on:
1. A total sample counts cutoff (500)
2. The sample yield (total counts) being higher than that of the corresponding negative extraction control 

**phyloseq_COVID_paper_final.R:** Full analysis of the mycobiome data using Phyloseq and several other packages. 

**Generate_patient_clinical_data_table_final.R:** A standalone script for generating a clinical data table appropriate for publishing:
