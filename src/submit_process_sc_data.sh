#!/bin/bash

conda activate sc-target-evidence-env
cd /nfs/team205/ed6/bin/sc_target_evidence/src/

outdir=/nfs/team205/ed6/bin/sc_target_evidence/data/

disease_ids=$(cat /nfs/team205/ed6/bin/sc_target_evidence/data/all_diseaseIDs.txt)

for d in $disease_ids; do \
    echo "python process_sc_data.py ${d} --output_dir ${outdir} --keep_all_genes True"  | \
        bsub -G teichlab -o logfile-pbulk-%J.out -e logfile-pbulk-%J.err -M100000 -R "select[mem>100000] rusage[mem=100000]" 
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

## Make diagnostic plots
for d in $disease_ids; do \
    d2=$(echo $d | sed 's/:/_/')
    if [ ! -f ${outdir}/cellxgene_targets_${d2}.pbulk_all_genes.h5ad ]; then
        echo "${d} missing"
    else 
        echo "python plot_diagnostics.py ${d} --plot_dir /nfs/team205/ed6/bin/sc_target_evidence/data/plots/" | \
            bsub -G teichlab -o logfile-diagnostic-%J.out -e logfile-diagnostic-%J.err -M50000 -R "select[mem>50000] rusage[mem=50000]" 
    fi
done

de_ready_disease_ids=$(ls $outdir/cellxgene_*pbulk_all_genes.h5ad | cut -f 9 -d '/' | sed 's/cellxgene_targets_//' | sed 's/.pbulk_all_genes.h5ad//')
for d in $de_ready_disease_ids; do \
    mv /nfs/team205/ed6/bin/sc_target_evidence/data/plots/cellxgene_${d}.celltype_harmonization.* /nfs/team205/ed6/bin/sc_target_evidence/data/plots/${d}_*
done

## Run DE analysis 
for d in $de_ready_disease_ids; do \
    echo "python run_de_ct.py ${d}" | bsub -G teichlab -q long -o logfile-de-ct-%J.out -e logfile-de-ct-%J.err -M50000 -R "select[mem>50000] rusage[mem=50000]"
        echo "python run_de_disease.py ${d}" | bsub -G teichlab -q long -o logfile-de-disease-%J.out -e logfile-de-disease-%J.err -M50000 -R "select[mem>50000] rusage[mem=50000]"
done














