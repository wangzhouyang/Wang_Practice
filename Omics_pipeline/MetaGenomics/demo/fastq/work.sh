#dir="/ifs1/ST_MD/USER/fangchao/BGISEQ/igA.n20.CTRL.Nov2016/00.DataSource/igA.10.TEST/soapA/"
#q=`sed 's/$/\|/' ../ref/gene.name.lst|tr -d '\n'|sed 's/|$//'`
#zcat $dir/NEW10BF0.gene.build/*gz|grep -E "$q"|cut -f 1 > hit.id.tmp
#hitId=`sed 's/$/\|//' hit.id.tmp|tr -d '\n'|sed 's/|$//'`
#perl -e 'while(<>){chomp;@a=split;$path=`ls /ifs1/ST_MD/USER/fangchao/BGISEQ/igA.n20.CTRL.Nov2016/00.DataSource/igA.10.TEST/soapA/$a[1].gene.build/*gz`;chomp($path);print "$a[0]\t$path\t$a[2]\n"}' sample.tmp.lst > sam.tmp.lst
#zcat /ifs1/ST_MD/USER/fangchao/BGISEQ/igA.n20.CTRL.Nov2016/00.DataSource/igA.10.TEST/soapA/ILM22BX0.gene.build/ILM22BX0.soap.SE.se.gz | sort -k1,1 | gzip > ILM22BX0.soap.SE.sort.gz &
#zcat /ifs1/ST_MD/USER/fangchao/BGISEQ/igA.n20.CTRL.Nov2016/00.DataSource/igA.10.TEST/soapA/ILM01BX0.gene.build/ILM01BX0.soap.SE.se.gz| sort -k1,1 | gzip > ILM01BX0.soap.SE.sort.gz &
perl getFqDemo.pl ../ref/gene.name.lst ID.soap.fq.lst
