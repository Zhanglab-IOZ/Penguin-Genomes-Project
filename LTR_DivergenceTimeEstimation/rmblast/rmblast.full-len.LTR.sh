source ~/.bashrc

genome='ADL.v6b.fa'
concensus='ADL.earlGrey.families.ERVK-full.fa'
species='ADL.v6b'

mkdir -p fas

less ${species}.full-ERVK.rmblastn.out.full-len-ltr.no-overlap.sort-u.bed1.txt | while read ccc sss eee strnd id aaa; do samtools faidx ./${species}.fa ${ccc}:${sss}-${eee} | sed "s/>/>${id} /g" >./fas/${id}.fa; done

ls ./fas/*.fa | while read aaa; do echo -e "blastn -subject ${aaa} -query ${aaa} -outfmt 6 "; done >ERVK.genomic-self.blastn.sh

bash ERVK.genomic-self.blastn.sh >${species}.full-ERVK.genomic-self.blastn.out

awk '{FS=OFS="\t"}FILENAME==ARGV[1]{aaa[$1]=$2"\t"$3}FILENAME==ARGV[2]{if($1 in aaa) print $0,aaa[$1]}' <(less ./${concensus}.size.len_ltr.txt | sed "s/#/\t/g" | cut -f 1,3,4) <(less ${species}.full-ERVK.genomic-self.blastn.out | awk '{FS=OFS="\t"}{print $1,$0}' | sed "s/\./\t/1" | cut -f 1,3-) | awk '{if($8<$15 && $9<$15*1.5) print }' >${species}.full-ERVK.genomic-self.blastn.LTR_LR.txt

less ${species}.full-ERVK.genomic-self.blastn.LTR_LR.txt | cut -f 2,8,9,10,11 >${species}.full-ERVK.genomic-self.blastn.LTR_LR.v1.txt

less ${species}.full-ERVK.genomic-self.blastn.LTR_LR.v1.txt | cut -f 1 | sort | uniq -c | sort -k1,1nr | awk '{if($1 > 1) print}'
