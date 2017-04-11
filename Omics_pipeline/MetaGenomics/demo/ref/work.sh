# some basic tools required, you may add the our Binboard dir into your PATH:
# export PATH=$PATH:/ifs1/ST_MD/PMO/script/bin/

# make ref demo (100 gene)
#zcat $MLIB/760MetaHIT_139HMP_368PKU_511Bac_IGC.fa.gz|perl -e '$/=">";<>;for(1..100){$_=<>;chomp;print ">$_"}' > IGC.100.demo.fa
# make index
gzip -d IGC.100.demo.fa
2bwt-builder IGC.100.demo.fa

grep ">" IGC.100.demo.fa|sed 's/>//'|awk '{print $1}' > gene.name.lst

gzip -d human.chr11.fa.gz
2bwt-builder human.chr11.fa
