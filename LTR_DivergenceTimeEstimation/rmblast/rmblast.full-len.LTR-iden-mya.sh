source ~/.bashrc

genome='ADL.v6b.fa'
concensus='ADL.earlGrey.families.ERVK-full.fa'
species='ADL.v6b'

mut_rate=0.000000015332
gen_time=12

mkdir -p LTR.fas

less ${species}.full-ERVK.genomic-self.blastn.LTR_LR.v1.txt | while read erv1 s1 e1 s2 e2; do
	samtools faidx ./fas/${erv1}.fa ${erv1}:${s1}-${e1} | sed "s/>/>L /" >./LTR.fas/${erv1}.fa
	samtools faidx ./fas/${erv1}.fa ${erv1}:${s2}-${e2} | sed "s/>/>R /" >>./LTR.fas/${erv1}.fa
done

cd LTR.fas

ls *.fa | grep -v "mafft" | while read aaa; do
	samtools faidx ${aaa} L >${aaa}.L
	samtools faidx ${aaa} R >${aaa}.R
	blastn -subject ${aaa}.L -query ${aaa}.R -outfmt 6 -out ${aaa}.L-R.blastn.out
done

ls *.fa | grep -v "mafft" | while read aaa; do
	bbb=$(echo ${aaa} | sed "s/\.fa//g")
	less ${bbb}.fa.L-R.blastn.out | sed "s/L/${bbb}.L/;s/R/${bbb}.R/" | sort -k12,12nr | head -1
done >../${species}.full-ERVK.genomic-self.blastn.LTR_LR.v1.self.out

cd ..

awk '{FS=OFS="\t"}FILENAME==ARGV[1]{aaa[$1]=$2}FILENAME==ARGV[2]{if($5 in aaa) print $1,$2,$3,$4,$5,$6,$7,aaa[$5]}' <(less ./${species}.full-ERVK.genomic-self.blastn.LTR_LR.v1.self.out | awk '{OFS="\t"}{print $1,$5/$4/(2*'${mut_rate}')/1000000*'${gen_time}'}' | sed "s/\.R//g") ${species}.full-ERVK.rmblastn.out.full-len-ltr.no-overlap.bed1.txt >${species}.full-ERVK.genomic-self.blastn.LTR_LR.v1.self.mya.bed1.txt
