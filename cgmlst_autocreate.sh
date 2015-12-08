#!/bin/bash

if test $# -eq 0; then
    echo "No arguments provided!"
    exit 1
fi

while test $# -gt 0; do

    case "$1" in
        
        -h|--help)
            echo "Required arguments:"
            echo "--work-dir         Working directory for this script."
            echo "--reference        Genome within --genomes path to be used as reference."
            echo "--genomes          Path to directory containing genomes as FASTAS."
            exit 0
            ;;
        
        --work-dir)
            shift
            if test $# -gt 0; then
                export WORKDIR=$1
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

    FASTANAME=$( head -n 1 $1 | cut -d' ' -f1 )
    FASTANAME=${FASTANAME}.fasta
    sed -i 's/^\(>.*\)/>1/' $1
    mv $1 $FASTANAME

}

### Set up directories ###

mkdir $WORKDIR
cd $WORKDIR 

mkdir alleles blast_out jsons temp

### Get non-redundant gene set ###

PROKKA_PREFIX=$( echo $REFERENCE | sed 's/\(.*\)\..*/1/' ) 

prokka --outdir prokka_out/ \
       --prefix $PROKKA_PREFIX \
       --locustag $PROKKA_PREFIX \
       --cpus 0
       $GENOMES/$REFERENCE


### all-vs-all BLAST search to filter homologues
makeblastdb -in $PROKKA_PREFIX.ffn -dbtype nucl \
            -out temp/${PROKKA_PREFIX}_db

blastn -db temp/${PROKKA_PREFIX}_db \
       -query prokka_out/${PROKKA_PREFIX}.ffn \ 
       -num_threads $( nproc ) \ 
       -outfmt 6 \
       -out blast_out/all_vs_all.tsv

python ava.py prokka_out/${PROKKA_PREFIX}.ffn blast_out/all_vs_all.tsv blast_out/non_redundant.fasta

### Create .markers file for MIST ###

csplit --prefix alleles/ -z prokka_out/${PROKKA_PREFIX}.ffn '/>/' '{*}'
cd alleles/

for i in $( ls ); do
    fasta_rename $i
done

cd $WORKDIR

python marker_maker.py --fastas alleles/ \ 
                       --out cgmlst.markers \
                       --test cgmlst

### run MIST ###

parallel mono $( locate MIST.exe | sort -n | tail -n 1 ) \ 
         -t cgmlst.markers \
         -T temp/ \
         -a alleles/ \
         -b -j jsons/{/.}.json \
         {} ::: ${GENOMES}/*.fasta


