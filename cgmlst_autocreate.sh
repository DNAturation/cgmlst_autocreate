#!/bin/bash

if test $# -eq 0; then
    echo "No arguments provided!"
    exit 1
fi

SRC_DIR=$(cd "$(dirname "$0")/.."; pwd)
SCRIPT_DIR=$(cd "$(dirname "$0")"; pwd)
T="$(date +%s)"
function fullpath {

    echo "$( cd "$(dirname "$1")"; pwd)/$(basename "$1")"
}


while test $# -gt 0; do

    case "$1" in
        
        -h|--help)

            echo ""
            echo "Required arguments:"
            echo "--work-dir         Working directory for this script."
            echo "--genomes          Path to directory containing genomes as FASTAs."
            echo ""
            echo "Optional arguments:"
            echo "--cores            Number of CPU cores to use where possible. (default = all)"
            echo "--fragment         Fragment size (bp) for Panseq to generate a pangenome. (default = 500)"
            echo "--percentid        Percent identity cutoff for pangenome. (default = 85%)"
            echo "--carriage         Percent presence for a gene to be considered core. (default = 99.5)"
            echo ""
            exit 0
            ;;
        
        --work-dir)
            shift
            if test $# -gt 0; then
                
                export WORKDIR=$(echo $(fullpath "$1"))/
            else
                echo "You need to give me a work directory."
                exit 1
            fi
            shift
            ;;

        --genomes)
            shift
            if test $# -gt 0; then
                export GENOMES=$(echo $(fullpath "$1"))/
            else
                echo "No genomes specified!"
                exit 1
            fi
            shift
            ;;

        --cores)
            shift
            if test $# -gt 0; then
                export CORES=$1
            else
                export CORES=$( grep -ci processor /proc/cpuinfo )
            fi
            echo "Using $CORES cores."
            shift
            ;;

        --fragment)
            shift
            if test $# -gt 0; then
                export FRAGMENT=$1
            else
                export FRAGMENT=500
                echo "Using default fragment size of 500 bp."
            fi
            shift
            ;;

        --percentid)
            shift
            if test $# -gt 0; then
                export PERCENTID=$1
            else
                export PERCENTID=85
                echo "Using default percent identity of 85."
            fi
            shift
            ;;
            
        --carriage)
            shift
            if test $# -gt 0; then
                export CARRIAGERATE=$1
            else
                export CARRIAGERATE=99.5
                echo "Using default carriage threshold of 99.5%"
            fi
            shift
            ;;

    esac
done

# Make sure mandatory args are given
if [ -z ${GENOMES+x} ] || [ -z ${WORKDIR+x} ]; then
    echo "Required arguments unset!"
    exit 1
fi

# Set up and test Panseq
if [ ! -f ${SCRIPT_DIR}/panseq ]; then
    echo "Setting up panseq..."
    git clone https://github.com/chadlaing/panseq ${SCRIPT_DIR}/panseq/
    perl ${SCRIPT_DIR}/panseq/Build.PL
    ${SCRIPT_DIR}/panseq/Build installdeps
    perl ${SCRIPT_DIR}/panseq/t/output.t || exit 1
fi

function fasta_rename {

    FASTANAME=$( head -n 1 $1 | cut -d' ' -f1 | sed 's/>//')
    FASTANAME=${FASTANAME}.fasta
    sed -i 's/^\(>.*\)/>1/' $1
    mv $1 $FASTANAME

}

function split_args {

    OIFS=$IFS
    IFS=','
    arr=$1
    IFS=$OIFS 
    echo arr
}

### Set up directories ###

mkdir $WORKDIR
cd $WORKDIR 

mkdir alleles blast_out jsons msa temp

### Get non-redundant gene set ###

### Format arguments and run panseq ###

printf "\nRunning Panseq to find pan-genome\n"

$PANSEQ_ARGS = ${SCRIPT_DIR}/settings.txt

python format_panseq_args.py --output     $PANSEQ_ARGS \
                             --fragment   $FRAGMENT \
                             --querydir   $GENOMES \
                             --basedir    panseq_results \
                             --cores      $CORES \
                             --percent_id $PERCENTID \
                             --carriage   $CARRIAGERATE

perl ${SCRIPT_DIR}/panseq/lib/panseq.pl $PANSEQ_ARGS

### Tidy up Panseq output for MIST friendliness ###
sed -i 's/_([0-9]*\.*[0-9]*)//g' panseq_results/panGenomeFragments.fasta

### Create .markers file for MIST ###
printf "\nSplitting to discrete fastas.\n"
csplit --quiet --prefix alleles/ -z panseq_results/panGenomeFragments.fasta '/>/' '{*}'

cd alleles/

for i in $( ls ); do
    fasta_rename $i
done

cd $WORKDIR

printf "\nBuilding reference genome .markers file\n"
python ${SCRIPT_DIR}/marker_maker.py --fastas alleles/ \
                       --out wgmlst.markers \
                       --test wgmlst

### run MIST ###
printf "\nRunning MIST in parallel.\n"

# A hack that will find the shortest path to MIST
    # The notion is that it won't accidentally find 
    # debugging binaries buried in the project directory
parallel mono $( locate MIST.exe | sort -n | tail -n 1 ) \
         -t wgmlst.markers \
         -T temp/ \
         -a alleles/ \
         -b -j jsons/{/.}.json \
         {} ::: ${GENOMES}/*.fasta

### Update allele definitions ### 
printf "\nUpdating allele definitions.\n"
python ${SCRIPT_DIR}/update_definitions.py --alleles alleles/ \
                                           --jsons   jsons/ \
                                           --test    wgmlst

### Align genes with clustalo ###
printf "\nAligning genes.\n"
parallel clustalo -i {} -o msa/{/} ::: alleles/*.fasta

### Divide Reference-based calls into core, genome, accessory schemes ###
printf "\nParsing JSONs.\n"
python ${SCRIPT_DIR}/json2csv.py --jsons jsons/ \
                                 --test  wgmlst \
                                 --out   wgmlst_calls.csv

printf "\nDividing markers into core and accessory schemes.\n"
Rscript ${SCRIPT_DIR}/divide_schemes.R wgmlst_calls.csv wgmlst.markers

printf "\nScript complete at `date`\n"

T="$(($(date +%s)-T))"
printf "Total run time: %02d:%02d:%02d:%02d\n" "$((T/86400))" "$((T/3600%24))" "$((T/60%60))" "$((T%60))"

