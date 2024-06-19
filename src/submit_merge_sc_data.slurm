#!/bin/bash
#SBATCH --partition=owners,pritch
#SBATCH --time=0:30:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/groups/pritch/emma/slurm-merge_sc_data-%A_%a.out
#SBATCH --error=/scratch/groups/pritch/emma/slurm-merge_sc_data-%A_%a.err

# RUN WITH: sbatch --array=1-${#input_array[@]} submit_merge_sc_data.slurm

cd /oak/stanford/groups/pritch/users/emma/bin/sc_target_evidence/src/

datadir=/oak/stanford/groups/pritch/users/emma/bin/sc_target_evidence/data/
tissue_table='tissue2dataset_id.csv'
tissue_ids=$(cut -d',' -f 2 "${datadir}${tissue_table}" | tail -n +2 | sort | uniq)
readarray -t input_array <<< "$tissue_ids"

# Get the tissue for this task
t=${input_array[$SLURM_ARRAY_TASK_ID]}

# Run the Python script with the current tissue
output_file="${datadir}/cellxgene_targets_${t}.pbulk_all_genes.h5ad"
if [ ! -f "$output_file" ]; then
    python merge_sc_data_tissue.py ${t} --data_dir ${datadir} --target_table_file filtered_nelson_disease_relevant_tissues_06172024.csv
else
    echo "Output file $output_file already exists. Skipping..."
fi
