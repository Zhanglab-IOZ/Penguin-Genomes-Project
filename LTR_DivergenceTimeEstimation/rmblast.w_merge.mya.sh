source ~/.bashrc

cat ./rmblast/ADL.v6b.full-ERVK.genomic-self.blastn.LTR_LR.v1.self.mya.bed1.txt <(cat rmblast_bed-merge/rnd*ltr1-ltr2.mya.txt | sed "s/_/\t/g" | awk '{OFS="\t"}{print $3,$4,$5,".",$1"_"$2"."$3"_"$4"_"$5,$1"_"$2,"LTR/ERVK",$6}') | sortBed >ADL.v6b.ERVK.full+merge.mya.bed1.txt

cat ADL.v6b.ERVK.full+merge.mya.bed1.txt | sort -k1,1V -k2,2n -k4,4Vr -k8,8nr | sort -k1,1V -k2,2n -u | sort -k1,1V -k3,3n -k4,4Vr -k8,8nr | sort -k1,1V -k3,3n -u | sort -k1,1V -k2,2n >ADL.v6b.ERVK.full+merge.mya.no-overlap.bed1.txt
