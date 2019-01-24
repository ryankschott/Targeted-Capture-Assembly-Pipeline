#!/bin/bash

# Targeted capture of complete coding regions across divergent species
# Ryan K Schott, Bhawandeep Panesar, Daren C Card, Matthew Preston, Todd A Castoe, Belinda SW Chang

# Searches directories for stats_.*.csv, parses, and outputs to .csv

USAGE='Usage: bash parse_stats.sh [-o output.csv]'

# Parses given stats_.*.txt and writes to given arg as output
parser () {
  output=$1
  shift
  while (($#)); do
    entry1=$( echo $1 | grep -o 'stats_.*\.txt' )
    entry1=$( echo $entry1 | sed 's/stats_//' | sed 's/.txt//' )
    cd $entry1
    entry2=$( head $1 -n 1 | cut -d ' ' -f 1 )
    entry3=$( head $1 -n 3 | cut -d $'\n' -f 3 | cut -d ' ' -f 5 )
    entry3=$( echo $entry3 | cut -d '%' -f 1 | tr -cd [:digit:]. )
    cd ..
    echo $entry1,$entry2,$entry3 >> $output
    shift
  done
}

output="parse_stats.csv"

while getopts 'ho:' flag; do
  case "${flag}" in
    h) 
      echo "$USAGE"
      exit 1;
      ;;
    o)
      output="${OPTARG}"
      ;;
    *)
      error "Unexpected option ${flag}"
      ;;
  esac
done

echo "Sequence Reads,Bases Read,Percentage Mapped" > $output
params=$( ls -R | grep "stats_.*\.txt" )
parser $output $params
