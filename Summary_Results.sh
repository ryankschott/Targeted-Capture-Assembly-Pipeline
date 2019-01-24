#!/bin/bash

# Targeted capture of complete coding regions across divergent species
# Ryan K Schott, Bhawandeep Panesar, Daren C Card, Matthew Preston, Todd A Castoe, Belinda SW Chang

# Reads stats_*.txt, Summary_CompletenessUpper.csv, Summary_CompletenessACGT.csv
# and Summary_Depth.csv, parses desired components, and outputs them to $1

USAGE='Usage: bash Summary_Results.sh [-o output.csv]'

# Calculates average of array $1
calc_average () {
  awk '
    BEGIN {
	  for (i = 1; i <= ARGC; i++){
	    total+=ARGV[i]
	  }
	  printf "%f", total/(ARGC-1)
	}
  ' $@
}

output="Summary_Results.csv"

while getopts 'ho:' flag; do
  case $flag in
    h) 
      echo "$USAGE"
      exit 1;
      ;;
    o)
      output=$OPTARG
      ;;
    *)
      error "Unexpected option $flag"
      ;;
  esac
done

# Parses given stats_.*.txt into desired components and returns them in arrays
sequence=()
bases_read=()
percent_mapped=()
for stats_file in $( ls -R | grep "stats_.*\.txt" ); do
  entry1=$( echo $stats_file | grep -o 'stats_.*\.txt' )
  entry1=$( echo $entry1 | sed 's/stats_//' | sed 's/.txt//' )
  entry2=$( head $entry1/$stats_file -n 1 | cut -d ' ' -f 1 )
  entry3=$( head $entry1/$stats_file -n 3 | cut -d $'\n' -f 3 | cut -d ' ' -f 5 )
  entry3=$( echo $entry3 | cut -d '%' -f 1 | tr -cd [:digit:]. )
  sequence+=("$entry1")
  bases_read+=("$entry2")
  percent_mapped+=("$entry3")
done

# Parses given file with last line being averages and returns them in an array
uppercase=$( tail -n 1 "Summary_CompletenessUpper.csv" | cut -d ',' -f 2- )
IFS=',' read -r -a uppercase <<< ${uppercase//\"}
acgt=$( tail -n 1 "Summary_CompletenessACGT.csv" | cut -d ',' -f 2- )
IFS=',' read -r -a acgt <<< ${acgt//\"}
depth=$( tail -n 1 "Summary_Depth.csv" | cut -d ',' -f 2- )
IFS=',' read -r -a depth <<< ${depth//\"}

# Calculate averages for each column and store in array
averages=("Average")
averages+=( $( calc_average ${bases_read[@]} ) )
averages+=( $( calc_average ${percent_mapped[@]} ) )
averages+=( $( calc_average ${uppercase[@]} ) )
averages+=( $( calc_average ${acgt[@]} ) )
averages+=( $( calc_average ${depth[@]} ) )

# Output ze stats to file
echo "Sequence Reads,Bases Read,Percentage Mapped," \
     "Average Completeness Uppercase,Average Completeness ACGT," \
     "Average Depth of Coverage" > $output
count=0
while [ $count -lt ${#sequence[@]} ]; do
  echo "${sequence[$count]},${bases_read[$count]},${percent_mapped[$count]}," \
       "${uppercase[$count]},${acgt[$count]},${depth[$count]}" >> $output
  count=$(( count+1 ))
done
echo >> $output
count=0
while [ $count -lt ${#averages[@]} ]; do
  echo -n "${averages[$count]}," >> $output
  count=$(( count+1 ))
done
