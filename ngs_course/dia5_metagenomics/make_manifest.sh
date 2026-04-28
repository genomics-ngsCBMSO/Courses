#!/bin/bash

usage()
{
    echo "Usage: $0 -r <TYPE_OF_READS>"
    echo -e "\n"
    echo "  -r | --reads     	Type of reads: single or paired (Default = paired) (OPTIONAL)"
    echo "  -d | --directory    Path to the directory containing the FASTQ files (Default = $PWD) (OPTIONAL)"
    echo "  -o | --output    	Output file name (Default = manifest.csv) (OPTIONAL)"
    echo "  -h | --help    	Display help"

}

TYPE_OF_READS='paired'
DIRECTORY=$PWD 
OUTPUT_FILE='manifest.csv'

#Get parameters

while [[ "$1" > 0 ]]
do
  case $1 in
   -r| --reads)
    shift
    TYPE_OF_READS=$1
    shift
  ;;
   -d| --directory)
    shift
    DIRECTORY=$(realpath "$1")
    shift
  ;;
   -o| --output)
    shift
    OUTPUT_FILE=$1
    shift
  ;;
    -h|--help)
    usage
    exit
  ;;
  *)
    echo "Wrong parameter!!!"
  ;;
  esac
done

#Writting Manifest file

echo 'sample-id,absolute-filepath,direction' > header.txt

if [ $TYPE_OF_READS = "paired" ]; then

  echo -e "\n"

  echo "Making a manifest file of all the FASTQ files in $DIRECTORY (as paired-end reads)"

  echo -e "\n"

  ls -1 $DIRECTORY | grep .fastq | awk 'BEGIN {FS="_"} {print $1","}' | awk '!x[$0]++' > sample_names.txt

  ls -d "$DIRECTORY"/*fastq | awk '{if($1 ~ "_R1") {print $1",""forward"}}' > forward.paths

  ls -d "$DIRECTORY"/*fastq | awk '{if($1 ~ "_R2") {print $1",""reverse"}}' > reverse.paths

  paste -d "" sample_names.txt forward.paths > forward
  paste -d "" sample_names.txt reverse.paths > reverse

  cat header.txt forward reverse > $OUTPUT_FILE

  rm header.txt sample_names.txt forward.paths reverse.paths forward reverse

else
  
  echo -e "\n"

  echo "Making a manifest file of all the FASTQ files in $DIRECTORY (as single-end reads)"

  echo -e "\n"

  ls -1 $DIRECTORY | grep .fastq | awk 'BEGIN {FS="_"} {print $1","}' | awk '!x[$0]++' > sample_names.txt

  ls -d "$DIRECTORY"/*fastq | awk '{print $1",""forward"}' > forward.paths

  paste -d "" sample_names.txt forward.paths > forward
 
  cat header.txt forward > $OUTPUT_FILE

  rm header.txt sample_names.txt forward.paths forward

fi

echo "JOB DONE!!!!"
echo -e "\n"
echo "
           ___              |\            .---.             _
           ( o )            |'_\           \ V /            | |
           _| |_           _| |_           _| |_           _| |_
         .'_____'.       .'_____'.       .'_____'.       .'_____'.
       |\ /     \ /|   |\ /     \ /|   |\ /     \ /|   |\ /     \ /|
       |||  @ @  |||   |||  9 9  |||   |||  6 6  |||   |||  o o  |||
       \_\   =   /_/   \_\   -   /_/   \_\   o   /_/   \_\  ._.  /_/
        .-'-----'-.     .-'-----'-.     .-'-----'-.     .-'-----'-.
       (_   ___   _)   (_   ___   _)   (_   ___   _)   (_   ___   _)
         | |___| |       | |___| |       | |___| |       | |___| |
         |       |       |       |       |       |       |       |
         (___|___)       (___|___)       (___|___)       (___|___)
"
echo -e "\n"


