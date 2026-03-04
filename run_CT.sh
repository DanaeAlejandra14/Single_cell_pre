
python ~/Single_cell_pre/celltypist_export/Celltypi_anno.py \
  --base /STORAGE/csbig/sc_ADers/celltypist/unassigned_only_minimal \
  --model "/STORAGE/csbig/sc_ADers/Adult_Human_PrefrontalCortex.pkl" \
  --majority_voting \
  --min_prob 0.6 \
  --out_csv /STORAGE/csbig/sc_ADers/celltypist/unassigned_only_minimal/unassigned_celltypist_predictions.csv
