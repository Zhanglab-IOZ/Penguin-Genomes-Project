source ~/.bashrc

###prepare a size.len_ltr.txt file!

genome='ADL.v6b.fa'
concensus='ADL.earlGrey.families.ERVK-full.fa'
species='ADL.v6b'

faSize -tab -detailed ${concensus} >${concensus}.size
#less ${concensus}.size | sed "s/#/\t/g" | cut -f 1,3 | sort -k1,1V >${concensus}.size.txt

ln -s ../${concensus}.size
ln -s ../${concensus}.size.len_ltr.txt

awk '{FS=OFS="\t"}FILENAME==ARGV[1]{aaa[$1]=$2"\t"$3}FILENAME==ARGV[2]{if($1 in aaa) print $0,aaa[$1]}' ${concensus}.size.len_ltr.txt ${species}.full-ERVK.rmblastn.out | awk '{if($7<0.5*$14 && $13-$8<0.5*$14) print }' | sed "s/#/\t/g" | sort -k3,3V -k10,10n | awk '{FS=OFS="\t"}{if($10<$11)print $3,$10,$11,"+",$1"."NR,$0; else print $3,$11,$10,"-",$1"."NR,$0}' >${species}.full-ERVK.rmblastn.out.full-len-ltr.bed1.txt

bedtools merge -i ${species}.full-ERVK.rmblastn.out.full-len-ltr.bed1.txt | bedtools intersect -wao -a - -b ${species}.full-ERVK.rmblastn.out.full-len-ltr.bed1.txt | awk '{if($2==$5 && $3==$6) print }' | cut -f 4- >${species}.full-ERVK.rmblastn.out.full-len-ltr.no-overlap.bed1.txt

less ${species}.full-ERVK.rmblastn.out.full-len-ltr.no-overlap.bed1.txt | sort -k1,1V -k2,2n -k9,9nr | sort -k1,1V -k2,2n -u >${species}.full-ERVK.rmblastn.out.full-len-ltr.no-overlap.sort-u.bed1.txt
less ${species}.full-ERVK.rmblastn.out.full-len-ltr.bed1.txt | sort -k1,1V -k2,2n -k9,9nr | sort -k1,1V -k2,2n -u >${species}.full-ERVK.rmblastn.out.full-len-ltr.sort-u.bed1.txt
