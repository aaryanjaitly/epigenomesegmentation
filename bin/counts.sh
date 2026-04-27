#!/bin/bash
script_name="counts.sh"
set -e
set -o pipefail

usage() {
  echo "Usage: $0 -t <input_sheet> -f <output_stem> -o <output_directory>"
  exit 1
}
timestamp() {
	  date +"%T"
}

printHelp() {
   echo -e "${bold}Description:${normal}"
   echo -e "This script generates count matrix for the ChIP-seq data. The inputs of this script are binsize and tab delimited file where first column is the name of histone mark and second column is the path to the corresponding bam file. It outputs a count matrix for provided histone marks in a non-overlapping window of given binsize. All histone marks with same name would be treated as replicates and merged together."
   echo -e "* For all the histone marks please provide corresponding .bam file for each mark."
   echo -e ""
   echo -e "${bold}Usage:${normal}"
   echo -e "bash $script_name  ${bold}$(tput setaf 1)-t Input sheet containing epigenetic marks  -o output directory -f output file stem name [-b binsize ] [-d input directory] [-p paired end ignore, filter, midpoint] [-c number of cores] [-g genome assembly]"
   echo -e ""
   echo -e "${bold}${red}Please give always absolute paths!${normal}"
   echo -e "${bold}Mandatory:$(tput sgr0)"
   echo -e "  -t Input sheet- Epigenetic marks label (first column) and corresponding histone mark bam file (second column )as tab seperated sheet." 
   echo -e "	* See sample input sheet "
   echo -e "	* Either provide complete path to the files or "
   echo -e "	* use -d argument to provide the complete path to the folder containing all marks and provide only file names in the input sheet. "
   echo -e "  -f output file stem name or sample name which will be used as prefix for the files generated during mixute model segmentation."
   echo -e "  -o complete path to the output directory where the output files would be generated. "
   echo -e "${bold}Optional:${normal}"
   echo -e "  -b binsize in base pairs. Default:200"
   echo -e "  -d complete path to the directory containing input files"
   echo -e "  -s number of basepairs to be shifted from 5' direction while counting reads for a particular bin. Only necessary for single end. Default: 75 bp"
   echo -e "  -p ignore, filter, midpoint Default: ignore"
   echo -e "		ignore - treat like single end"
   echo -e "		filter - 5’-end of first read in a properly aligned pair"
   echo -e "		midpoint - consider the midpoint of an aligned fragment)"
   echo -e "  -c number of cores to be used for generating count matrices. Default: 4"
   echo -e "  -g choosing apppropriate genome. Example inputs: "
   echo -e "		mm10 - mouse genome version 10"
   echo -e "		hg19 - Human genome version 19 (Ensembl-37)"
   echo -e "		hg38 - Human genome version 38 (Ensembl-38)"
   echo -e "  -r regions files that should be used for creating count matrix. Contains three column:"
   echo -e "		first column  - name of the chromosome. Eg. chr1"
   echo -e "		second column - start of the chromosome. Eg. 1"
   echo -e "		third column  - length of the chromosome in bp. Eg. 248956422"
}

# defiult values 
binsize=200
shift_bp=75
paired_end=ignore
cores=4

while getopts ":ht:o:c:b:s:p:f:d:g:r:" OPTION
do
        case $OPTION in
                h)
                        printHelp; exit 0 ;;
                t)
                        tab_file="$OPTARG" ;;
                o)
                        output_dir=$OPTARG ;;
                c)
                        cores=$OPTARG ;;
                b)
                        binsize=$OPTARG ;;
                s)
                        shift_bp=$OPTARG ;;
                p)
                        paired_end=$OPTARG ;;
                f) 
                        filename=$OPTARG ;; 
                d) 
                        input_dir=$OPTARG ;;
                g) 
                        genome=$OPTARG ;;
                r)
                        regions_chromosome=$OPTARG ;;
  esac
done

echo -e "${bold}${magenta}The following command was used:${normal} ${bold}bash mixture_model_container.sh" $@ ${normal}

# Output folder check if the user did not create the output folder than it would be created 
if [ ! -d $output_dir ]
	then
		mkdir -p ${output_dir}
	else
		echo -e "Failed to create an output directory! Directory may exist, in this case files will be overwritten." 
fi


if [ -z ${genome} ]
then 
        output_file=$(realpath -s ${output_dir}/${filename})
else    
        output_file=$(realpath -s ${output_dir}/${filename}_${genome})
fi


# if user provided input directory and paths of the bam file then create a new tab file 
if [ -z ${input_dir} ]
then
        echo -e "Input file to be used: ${tab_file}"
else 
	
        while read marks files 
        do 
                path=$(realpath -s ${input_dir}/${files})
                echo -e "${marks}\t${path}"
        done < ${tab_file} >  ${output_file}_episegmix_input.txt
        tab_file=$(realpath -s ${output_file}_episegmix_input.txt)
        echo -e "Input file to be used: ${tab_file}"
fi

echo -e "$(timestamp) Checking for the input files"
# 1. Check all the paths and the bam files works 
bam_files=($(cut -f 2 ${tab_file}))
number_of_bams=${#bam_files[@]}
j=0
while [ "$j" -le "$((number_of_bams-1))" ] 
        do      
                bamfile=${bam_files[$j]}
                [ -f ${bamfile} ] || { echo "File '${bamfile}' not found."; exit 1; }
                [ -f ${bamfile}.bai ] || { samtools index ${bamfile} ; } 
                j=$(($j + 1))
        done	
# 2. Create regions file if user did not provide
echo -e "$(timestamp) Finalizing the regions for the count matrix"
if [ -z ${regions_chromosome} ]
then 
        samtools view -H ${bamfile} 	\
                | grep "@SQ" \
                |cut -f 2,3 \
                |cut -d ":" -f 2,3 \
                |awk -F"\t|:" '{print $1,$3}' \
                |awk -vOFS='\t' '{print $1,1,$2}' \
                | awk '$1~/^chr[1-9XY][0-9]$/ || $1~/^chr[1-9XY]$/ ||$1~/^[1-9XY][0-9]$/ || $1~/^[1-9XY]$/' \
                > ${output_file}_chrom_sizes.bed
        chrom_sizes=${output_file}_chrom_sizes.bed
else    
        echo "user provided chromosome regions will be used for generating count matrices"
        chrom_sizes=${regions_chromosome}
fi


# 4. Run the script and generate the count matrix in the output folder
echo -e "$(timestamp) Reading bam files and generating count matrix"
counts.R ${tab_file} ${chrom_sizes} ${binsize} ${shift_bp} ${output_file} ${cores} ${paired_end} 

# 5. Create a count matrix for the EpiSegMix input
counts_txt=${output_file}"_counts.txt"
refine_chr_regions=${output_file}"_refine_chr_regions.bed"
refine_regions=${output_file}".bed"
episegmix_input=${output_file}"_refined_counts.txt"
bedtools makewindows -b ${refine_chr_regions} -w ${binsize} > ${refine_regions}

paste <(awk -v OFS="\t" 'BEGIN{print "chr\tstart\tend"} {gsub(/chr/,""); print$0}' ${refine_regions}) \
        <(awk -vOFS="\t" '{print}' ${counts_txt})  \
	> ${episegmix_input}

# 6. Remove all the temporary files created in between
rm ${counts_txt} ${chrom_sizes} ${refine_chr_regions}


echo -e "$(timestamp) Count matrix generated and saved as ${episegmix_input}"