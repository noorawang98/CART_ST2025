#python batch_destvi_bin50_cli_opt.py --obs-csv /home/ST_Data/scData/RNA/obs_coarse_label.csv --out-dir /home/ST_Data/cellbin/destvi --st-dir  /home/ST_Data/cellbin --obs-index-is-barcode --st-glob *cellbin*h5ad --label-col leiden --reference-h5ad /home/ST_Data/scData/RNA/destvi_reference/scrna_reference_for_destvi.h5ad  > /home/ST_Data/cellbin/master.stdout.log 2>&1 
#python batch_destvi_bin50_cli_opt.py --obs-csv /home/ST_Data/scData/RNA/obs_coarse_label.csv --out-dir /home/ST_Data/cellbin/destvi_bin50_2_leiden --st-dir  /home/ST_Data/bin50 --obs-index-is-barcode --label-col leiden --reference-h5ad /home/ST_Data/scData/RNA/destvi_reference/scrna_reference_for_destvi.h5ad  > /home/ST_Data/bin50/master.stdout.log 2>&1 
#python batch_destvi_bin50_cli_opt.py --obs-csv /home/ST_Data/scData/RNA/obs_coarse_label.csv --out-dir /home/ST_Data/ljx/destvi_cellbin --st-dir  /home/ST_Data/ljx/human/st --obs-index-is-barcode  --label-col leiden --reference-h5ad /home/ST_Data/scData/RNA/destvi_reference/scrna_reference_for_destvi.h5ad  > /home/ST_Data/ljx/human/st/master.stdout.log 2>&1

python batch_destvi_stmode_scrna.py \
  --st-mode cellbin \
  --st-dir /home/ST_Data/cellbin \
  --out-dir /home/ST_Data/cellbin/destvi_leiden \
  --scrna-input h5ad \
  --scrna-h5ad /home/ST_Data/scData/RNA/destvi_reference/scrna_reference_for_destvi.h5ad \
  --obs-csv /home/ST_Data/scData/RNA/obs_coarse_label.csv \ 
  --label-col leiden \
  --cond-label-source original \
  --batch-col sample \
  --st-symbol-col real_gene_name \
  --st-gene-name-source symbol_col \
  --retrain-reference 
#--st-gene-name-source var_names \
#  --retrain-reference

python batch_destvi_stmode_scrna.py \
  --st-mode cellbin \
  --st-dir /home/ST_Data/ljx/human/st \
  --out-dir /home/ST_Data/ljx/human/st/destvi_leiden \
  --scrna-input h5ad \
  --scrna-h5ad /home/ST_Data/scData/RNA/destvi_reference/scrna_reference_for_destvi.h5ad \
  --obs-csv /home/ST_Data/scData/RNA/obs_coarse_label.csv \
  --label-col leiden \
  --cond-label-source original \
  --batch-col sample \
  --st-symbol-col real_gene_name \
  --st-gene-name-source symbol_col \
  --retrain-reference
 # --st-gene-name-source var_names \
  #--retrain-reference

#/home/ST_Data/ljx/human/st
