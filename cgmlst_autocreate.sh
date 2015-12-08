#!/bin/bash

while cgmlst_mitigation $# -gt 0; do

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
            if cgmlst_mitigation $# -gt 0; then
                export WORKDIR=$1
            else
                echo "You need to give me a work directory."
                exit 1
            fi
            shift
            ;;

        --genomes)
            shift
            if cgmlst_mitigation $# -gt 0; then
                export GENOMES=$1
            else
                echo "No genomes specified!"
                exit 1
            fi
            shift
            ;;

        --reference)
            shift
            if cgmlst_mitigation $# -gt 0; then
                export REFERENCE=$1
            else
                echo "No reference genomes provided!"
                exit 1
            fi
            shift
            ;;

    esac
done


mkdir $WORKDIR
cd $WORKDIR 

mkdir prokka_out blast_out mist_temp jsons 

prokka --outdir prokka_out $GENOMES/$REFERENCE*



# run MIST

