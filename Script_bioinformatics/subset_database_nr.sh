taxonkit list -j 2 --ids 2 --indent "" --data-dir /data/mengqy/database/taxonomy > Bacteria.list &

cat /data/mengqy/database/taxonomy/prot.accession2taxid |csvtk -t grep -f taxid -P Bacteria.list |csvtk -t cut -f accession.version >Bacteria.taxid.acc.txt &

blastdb_aliastool -seqidlist Bacteria.taxid.acc.txt -db /data/liang/nr.database/nr -out nr_Bacteria -title nr_Bacteria &

blastdbcmd -db /mnt/nas2/database/ncbi/Nr/nr_Bacteria -entry all -dbtype prot -out nr_Bacteria.fa
