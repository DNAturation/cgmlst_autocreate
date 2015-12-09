#!/bin/bash

if test $# -eq 0; then
    echo "No arguments provided!"
    exit 1
fi

SRC_DIR=$(cd "$(dirname "$0")/.."; pwd)
SCRIPT_DIR=$(cd "$(dirname "$0")"; pwd)

while test $# -gt 0; do

    case "$1" in
        
        -h|--help)
            echo "Required arguments:"
            echo "--work-dir         Working directory for this script."
            echo "--reference        Fasta filename (not path) within --genomes path to be used as reference."
            echo "--genomes          Path to directory containing genomes as FASTAS."
            exit 0
            ;;
        
        --work-dir)
            shift
            if test $# -gt 0; then
                
                export WORKDIR=${SRC_DIR}/$1
            else
                echo "You need to give me a work directory."
                exit 1
            fi
            shift
            ;;

        --genomes)
            shift
            if test $# -gt 0; then
                export GENOMES=$1
            else
                echo "No genomes specified!"
                exit 1
            fi
            shift
            ;;

        --reference)
            shift
            if test $# -gt 0; then
                export REFERENCE=$1
            else
                echo "No reference genomes provided!"
                exit 1
            fi
            shift
            ;;

    esac
done

function fasta_rename {

    FASTANAME=$( head -n 1 $1 | cut -d' ' -f1 | sed 's/>//')
    FASTANAME=${FASTANAME}.fasta
    sed -i 's/^\(>.*\)/>1/' $1
    mv $1 $FASTANAME

}

### Set up directories ###

mkdir $WORKDIR
cd $WORKDIR 

mkdir alleles blast_out jsons temp

### Get non-redundant gene set ###

PROKKA_PREFIX=$( echo $REFERENCE | sed 's/\(.*\)\..*/\1/' ) 

prokka --outdir prokka_out/ \
       --prefix $PROKKA_PREFIX \
       --locustag $PROKKA_PREFIX \
       --cpus 0 \
       ${SRC_DIR}/${GENOMES}${REFERENCE}

### all-vs-all BLAST search to filter homologues
makeblastdb -in prokka_out/${PROKKA_PREFIX}.ffn -dbtype nucl \
            -out temp/${PROKKA_PREFIX}_db
blastn -query prokka_out/${PROKKA_PREFIX}.ffn \
       -db temp/${PROKKA_PREFIX}_db \
       -num_threads $( nproc ) \
       -outfmt 10 \
       -out blast_out/all_vs_all.csv

# AVA filters the all-vs-all search for homologues
    # Thresholds are (currently) hardcoded as 90% PID and 50% length 
    # Only the longest variant is kept
python ${SCRIPT_DIR}/ava.py --seq prokka_out/${PROKKA_PREFIX}.ffn \
                            --result blast_out/all_vs_all.csv \
                            --out blast_out/non_redundant.fasta

### TODO Maybe add filter to get genes for a cgMLST, wgMLST, or accessory scheme

### Create .markers file for MIST ###

csplit --quiet --prefix alleles/ -z prokka_out/${PROKKA_PREFIX}.ffn '/>/' '{*}'

cd alleles/

for i in $( ls ); do
    fasta_rename $i
done

cd $WORKDIR

python ${SCRIPT_DIR}/marker_maker.py --fastas alleles/ \
                       --out cgmlst.markers \
                       --test cgmlst

### run MIST ###

# A hack that will find the shortest path to MIST
    # The notion is that it won't accidentally find 
    # debugging binaries buried in the project directory
parallel mono $( locate MIST.exe | sort -n | tail -n 1 ) \
         -t cgmlst.markers \
         -T temp/ \
         -a alleles/ \
         -b -j jsons/{/.}.json \
         {} ::: ${SRC_DIR}/${GENOMES}/*.fasta

### Update allele definitions ### 
python ${SCRIPT_DIR}/update_definitions.py --alleles alleles/ \
                                           --jsons jsons/ \
                                            --test cgmlst

