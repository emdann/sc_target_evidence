
To extract single-cell evidence for each disease, we run the following:

```bash
disease_id=MONDO:0000000
outdir=/path/to/output/dir
plotdir=/path/to/plot/dir

## ---- Download relevant data from cellxgene ---- ##
# saves {outdir}/cellxgene_targets_{disease_id}.pbulk_all_genes.h5ad
python process_sc_data_v2.py ${disease_id} --output_dir ${outdir} --keep_all_genes True

## ---- Make diagnostic plots ---- ##
# saves plots in {plotdir}/{disease_id}_{disease_name}/
python plot_diagnostics.py ${disease_id} --plot_dir ${plotdir}

## ---- Run differential expression analysis ---- ##
# saves {outdir}/DE_celltype_{disease_id}.hvgs.csv
python run_de_ct.py ${disease_id} --plot_dir ${plotdir}
# saves {outdir}/DE_diseasecelltype_{disease_id}.hvgs.csv
# saves {outdir}/DE_diseasebulk_{disease_id}.hvgs.csv
python run_de_disease.py ${disease_id} --plot_dir ${plotdir}
```