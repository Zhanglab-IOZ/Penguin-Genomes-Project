source ~/.bashrc

species="ADL.v6b"
rmblastn_out="ADL.v6b.full-ERVK.rmblastn.out"
rmblastn_bed1="ADL.v6b.full-ERVK.rmblastn.out.full-len-ltr.no-overlap.sort-u.bed1.txt"
len_ltr="ADL.earlGrey.families.ERVK-full.fa.size.len_ltr.txt"

ln -s ../rmblast/${rmblastn_out}
ln -s ../${len_ltr}
ln -s /public3/home/wangzn/yuanhao/penguin/genome/ADL.v6b/ADL.v6b.fa

bedtools intersect -v -a <(less ${rmblastn_out} | awk '{OFS="\t"}{if($9<$10) print $2,$9-1,$10,$1}' | sort -k1,1V -k2,2n) -b ../rmblast/${rmblastn_bed1} >${species}.ERVK.rmblastn.partial.pos.bed
bedtools intersect -v -a <(less ${rmblastn_out} | awk '{OFS="\t"}{if($9>$10) print $2,$10-1,$9,$1}' | sort -k1,1V -k2,2n) -b ../rmblast/${rmblastn_bed1} >${species}.ERVK.rmblastn.partial.neg.bed

less ${len_ltr} | while read erv len ltr; do
	ervf=$(echo ${erv} | cut -d "#" -f 1)
	awk '{if($4 == "'${erv}'") print }' ${species}.ERVK.rmblastn.partial.pos.bed | bedtools merge -i - -d 1000 | awk '{OFS="\t"}{print $0,"'${len}'","'${ltr}'","'${erv}'"}' >${ervf}.merge.bed
done
less ${len_ltr} | while read erv len ltr; do
	ervf=$(echo ${erv} | cut -d "#" -f 1)
	awk '{if($4 == "'${erv}'") print }' ${species}.ERVK.rmblastn.partial.neg.bed | bedtools merge -i - -d 1000 | awk '{OFS="\t"}{print $0,"'${len}'","'${ltr}'","'${erv}'"}' >>${ervf}.merge.bed
done

less ${len_ltr} | while read erv len ltr; do
	ervf=$(echo ${erv} | cut -d "#" -f 1)
	less ${ervf}.merge.bed | awk '($3-$2)/$4 > 0.5{print $1":"$2"-"$3}' | samtools faidx ./${species}.fa -r - >${ervf}.merge_0.5.fa
done
