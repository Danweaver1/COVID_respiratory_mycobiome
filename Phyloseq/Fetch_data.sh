###Fetches data based on key file 
#To define batches to be processed, modify "*batches_passing.csv" file:
#file is extraction group key: "<sample ID>,<batch number>" 
#reduced to only those batches which passed qPCR NEC filter
#this file is used to fetch the count files

##species count file dir
data="../../R_analysis/raw_counts_scaled"

##count file labelling
count=".scaled_counts.csv" 

##clean up any previous analyses
if [ -d "scaled_counts" ]
  then
    rm -r scaled_counts/
fi

##fetch all ITS count data (only data specified by the batches passing file)
mkdir scaled_counts
cd scaled_counts
while read p
do
  file=$( echo $p | cut -d"," -f1 | sed 's/"//g' )
  if [ -e $data/"$file""$count" ]
  then
    cp $data/"$file""$count" .
  else
    echo "$data/"$file""$count" not found"
fi
done <../*batches_passing.csv
