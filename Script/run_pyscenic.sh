export DASK_DISTRIBUTED__WORKER__MEMORY__TERMINATE=0.8


dir=/data/nas1/mengqingyao_OD/project/Program-179/18_scRNA_SCENIC/cisTarget_databases

tfs=$dir/hs_hgnc_tfs.txt
f_db_500bp=$dir/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
f_db_10kb=$dir/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
tbl=$dir/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
input_loom=$dir/scobj.loom

echo "清理残留Dask进程..."
pkill -f "dask-scheduler" 2>/dev/null
pkill -f "dask-worker" 2>/dev/null

# 列出文件路径，确保它们存在
ls $tfs $f_db_500bp $f_db_10kb  $tbl

# 第一步：GRN推断
pyscenic grn \
--seed 777 \
--num_workers 10 \
--output adj.sample.tsv \
--method grnboost2 \
$input_loom $tfs
--client-timeout 360000
if [ ! -s adj.sample.tsv ]; then
 echo "Error: adj.sample.tsv is empty or not generated!"
exit 1
 fi
 
 
# 第二步：Cistarget分析
pyscenic ctx \
adj.sample.tsv \
$f_db_500bp $f_db_10kb \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 15 \
--mask_dropouts
if [[ -s reg.csv ]]; then
    echo "Cistarget成功生成reg.csv"
else
    echo "致命错误：reg.csv生成失败！请执行以下诊断步骤："
    ls -lh /tmp/dask-worker-*.log
    grep -A 50 'Traceback' /tmp/dask-worker-*.log
    exit 1
fi

# 第三步：AUCell分析
pyscenic aucell \
$input_loom \
reg.csv \
--output out_SCENIC.loom \
--num_workers 15 \
  --seed 777