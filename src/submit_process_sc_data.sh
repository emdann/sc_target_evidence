#!/bin/bash

conda activate sc-target-evidence-env
cd /nfs/team205/ed6/bin/sc_target_evidence/src/
outdir=/nfs/team205/ed6/bin/sc_target_evidence/data/

disease_ids=$(cat /nfs/team205/ed6/bin/sc_target_evidence/data/all_diseaseIDs.txt)

for d in $(cat failed_mondos.txt); do \
# for d in $disease_ids; do \
    echo "python process_sc_data.py ${d} --output_dir ${outdir} --keep_all_genes True"  | \
        bsub -G teichlab -o logfile-pbulk-%J.out -e logfile-pbulk-%J.err -M300000 -R "select[mem>300000] rusage[mem=300000]" 
done 

# Check which diseases have completed
for d in $disease_ids; do \
    d2=$(echo $d | sed 's/:/_/')
    if [ ! -f ${outdir}/cellxgene_targets_${d2}.pbulk_all_genes.h5ad ]; then
        echo "${d2} [ ]"
        echo "${d2}" >> failed_mondos.txt
    else 
        echo "${d2} [X]"
    fi
done

for d in $(cat failed_mondos.txt); do \
    echo "python process_sc_data.py ${d} --output_dir ${outdir}"  | \
        bsub -G teichlab -o logfile-failedpbulk-%J.out -e logfile-failedpbulk-%J.err -M250000 -R "select[mem>250000] rusage[mem=250000]" 
done 

## Make diagnostic plots
for d in $disease_ids; do \
    d2=$(echo $d | sed 's/:/_/')
    if [ ! -f ${outdir}/cellxgene_targets_${d2}.pbulk_all_OT_targets.h5ad ]; then
        echo "${d} missing"
    else 
        echo "python plot_diagnostics.py ${d} --plot_dir /nfs/team205/ed6/bin/sc_target_evidence/data/plots/" | \
            bsub -G teichlab -o logfile-diagnostic-%J.out -e logfile-diagnostic-%J.err -M50000 -R "select[mem>50000] rusage[mem=50000]" 
    fi
done

## Run DE analysis 
# for d in $(cat failed_mondos.txt); do \
for d in $disease_ids; do \
    d2=$(echo $d | sed 's/:/_/')
    if [ ! -f ${outdir}/DE_diseasebulk_${d2}.all_targets.csv ]; then
      echo "python run_de.py ${d}" | \
            bsub -G teichlab -o logfile-de-%J.out -e logfile-de-%J.err -M50000 -R "select[mem>50000] rusage[mem=50000]"
    fi
done