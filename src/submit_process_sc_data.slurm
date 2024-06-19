#!/bin/bash
#SBATCH --partition=owners,pritch
#SBATCH --time=8:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=3
#SBATCH --output=/scratch/groups/pritch/emma/slurm-process_sc_data_v2-%A_%a.out
#SBATCH --error=/scratch/groups/pritch/emma/slurm-process_sc_data_v2-%A_%a.err

cd /oak/stanford/groups/pritch/users/emma/bin/sc_target_evidence/src/

datadir=/oak/stanford/groups/pritch/users/emma/bin/sc_target_evidence/data/
filters_table='cellxgene_filters.ct_marker_evidence.csv'
tissue_dataset_ids=$(cut -d',' -f 1 "${datadir}${filters_table}" | tail -n +2)
readarray -t dataset_array <<< "$tissue_dataset_ids"

# Get the dataset for this task
d=${dataset_array[$SLURM_ARRAY_TASK_ID]}

# Run the Python script with the current dataset
output_file="${datadir}/cellxgene_targets_${d}.pbulk_all_genes.h5ad"
if [ ! -f "$output_file" ]; then
    python process_sc_data_v2.py ${d} --data_dir ${datadir}
else
    echo "Output file $output_file already exists. Skipping..."
fi
