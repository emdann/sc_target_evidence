#!/bin/bash

conda activate sc-target-evidence-env
cd /nfs/team205/ed6/bin/sc_target_evidence/src/
outdir=/nfs/team205/ed6/bin/sc_target_evidence/data/

disease_ids=$(cat /nfs/team205/ed6/bin/sc_target_evidence/data/all_diseaseIDs.txt)

for d in $disease_ids; do \
    echo "python process_sc_data.py ${d} --output_dir ${outdir}"  | \
        bsub -G teichlab -o logfile-pbulk-%J.out -e logfile-pbulk-%J.err -M50000 -R "select[mem>50000] rusage[mem=50000]" 
done 

# Check which diseases have completed
for d in $disease_ids; do \
    d2=$(echo $d | sed 's/:/_/')
    if [ ! -f ${outdir}/cellxgene_targets_${d2}.pbulk_all_OT_targets.h5ad ]; then
        echo "${d2} [ ]"
    else 
        echo "${d2} [X]"
    fi
done
