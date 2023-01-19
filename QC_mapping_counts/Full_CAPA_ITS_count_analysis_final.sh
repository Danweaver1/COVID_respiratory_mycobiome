##### Analysis of full ITS CAPA count dataset
### requires following R scripts:
##generate_scaled_cutoff_counts_v1.R
#also requires: ITS_run_group_key_batches_passing.csv


### To define batches to be processed, modify "Extraction*batches_passing.csv" file:
#file is extraction group key: "<sample ID>,<batch number>" 
#reduced to only those batches which passed qPCR NEC filter
#this file is used to fetch the count files

### Run to run scaling 
##before creating usual cutoff count files, samples are first subject to run-run scaling 
#creating scale counts: eg. "NEC13.scaled_counts.csv"
#and scaled cutoff counts: eg. "NEC13.scaled_counts_SPC_0_2_SRC_10.csv"
##this requires a key file indicating which sequencing run each sample belongs to ("ITS_run_group_key_batches_passing.csv")

##species counts file dir ********** SET DIRECTORY OF COUNT FILES HERE ***********
data=""

##count file labelling
count=".Rspecies_counts.csv"

##clean up any previous analyses
if [ -d "raw_counts" ]
  then
    rm -r raw_counts/
    rm -r merged_cutoff_data/
    rm *merge.csv
    rm *.tiff
fi

##fetch all ITS count data (only data specified by the batches passing file)
mkdir raw_counts
cd raw_counts
while read p
do
  file=$( echo $p | cut -d"," -f1 | sed 's/"//g' )
  if [ -e $data/"$file""$count" ]
  then
    cp $data/"$file""$count" .
  else
    echo "$data/"$file""$count" not found"
fi
done <../../Key_files/Extraction*batches_passing.csv

##create %species cutoff counts (raw)
Rscript ../generate_scaled_cutoff_counts_v1.R  > Routput.txt
#move into directory of scaled counts
cd ../raw_counts_scaled

#create dir for merged data
mkdir ../merged_cutoff_data
cat *.scaled_counts_SPC_0_2_SRC_10.csv | grep -v Counts > ../merged_cutoff_data/merge.csv
