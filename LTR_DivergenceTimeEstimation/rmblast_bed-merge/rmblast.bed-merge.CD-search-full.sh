source ~/.bashrc

ln -s ../ADL.earlGrey.families.ERVK-full.fa

len_ltr='./ADL.earlGrey.families.ERVK-full.fa.size.len_ltr.txt'
genome='./ADL.v6b.fa'
ERV1cons='./ADL.earlGrey.families.ERVK-full.fa'
mut_rate=0.000000015332
gen_time=12

#prepare a CD-search.local.ID.merge_0.5.fa.pbs.txt file

less ${len_ltr} | while read erv len ltr; do
	ervf=$(echo ${erv} | cut -d "#" -f 1)
	less ${ervf}.merge_0.5.cdd.out | grep "Query" | grep "Nucleotide" >${ervf}_merge_0.5.cdd.Query.site.txt
	less ${ervf}.merge_0.5.cdd.out | grep "Query_" | grep "Gag_" | awk '{print $2}' | awk -v FS="[" '{print $1}' | uniq >${ervf}_merge_0.5.cdd.gag.txt
	less ${ervf}.merge_0.5.cdd.out | grep "Query_" | grep -E "RNaseH|RNase_H" | awk '{print $2}' | awk -v FS="[" '{print $1}' | uniq >${ervf}_merge_0.5.cdd.RT.txt
	awk 'NR==FNR{a[$1]=$1;next}NR>FNR{if($1 in a)print $1}' ${ervf}_merge_0.5.cdd.gag.txt ${ervf}_merge_0.5.cdd.RT.txt | awk 'NR==FNR{a[$2]=$5;next}NR>FNR{if($1 in a) print a[$1]}' ${ervf}_merge_0.5.cdd.Query.site.txt - >${ervf}_merge_0.5.cdd.gag.RT.site.txt
	samtools faidx ${genome} -r ${ervf}_merge_0.5.cdd.gag.RT.site.txt >${ervf}_merge_0.5.cdd.gag.RT.site.fa
done

less ${len_ltr} | while read erv len ltr; do
	ervf=$(echo ${erv} | cut -d "#" -f 1)
	samtools faidx ${ERV1cons} ${erv}:1-${ltr} >${ervf}.LTR-L.fa
done

less ${len_ltr} | while read erv len ltr; do
	ervf=$(echo ${erv} | cut -d "#" -f 1)
	nucmer --maxmatch ${ervf}_merge_0.5.cdd.gag.RT.site.fa ${ervf}.LTR-L.fa --prefix ${ervf}.ltr_vs_cdd.gag.RT --threads 16
	show-coords ${ervf}.ltr_vs_cdd.gag.RT.delta -T >${ervf}.ltr_vs_cdd.gag.RT.coords.txt
	less ${ervf}.ltr_vs_cdd.gag.RT.coords.txt | awk '{if($5 > '${ltr}'/2 ) print $8}' | uniq -c | awk '$1==2{print $2}' | sort -k1,1V >${ervf}.cdd.gag.RT.2ltr.site.txt
done

less ${len_ltr} | while read erv len ltr; do
	ervf=$(echo ${erv} | cut -d "#" -f 1)
	qqq0=""
	i=0
	awk '{OFS="\t"}FILENAME==ARGV[1]{aaa[$1]=$1}FILENAME==ARGV[2]{if($8 in aaa) print }' ${ervf}.cdd.gag.RT.2ltr.site.txt ${ervf}.ltr_vs_cdd.gag.RT.coords.txt | awk '{if($5 > '${ltr}'/2) print }' | sort -k8,8V -k1,1n | while read qs qe ts te ql tl iden qqq ttt; do
		if [[ ${qqq} != ${qqq0} ]]; then
			qqq0=${qqq}
			i=1
			qs1=${qs}
		else
			i=$((${i} + 1))
			qe1=${qe}
			ccc=$(echo ${qqq} | sed "s/:/\t/g;s/-/\t/g" | cut -f 1)
			sss=$(echo ${qqq} | sed "s/:/\t/g;s/-/\t/g" | cut -f 2)
			sss1=$((${sss} + ${qs1} - 1))
			sss2=$((${sss} + ${qe1} - 1))
			sssl=$((sss2 - sss1 + 1))
			echo -e "${qqq}\t${ccc}:${sss1}-${sss2}\t${sssl}"
		fi
	done | awk '{if($3 > 0.5*'${len}' && $3 < 2*'${len}') print }' >${ervf}.cdd.gag.RT.2ltr.site.ltr-ltr.txt
done

less ${len_ltr} | while read erv len ltr; do
	ervf=$(echo ${erv} | cut -d "#" -f 1)
	awk '{OFS="\t"}FILENAME==ARGV[1]{aaa[$1]=$1}FILENAME==ARGV[2]{if($8 in aaa) print }' ${ervf}.cdd.gag.RT.2ltr.site.txt ${ervf}.ltr_vs_cdd.gag.RT.coords.txt | awk '{if($5 > '${ltr}'/2) print }' | sort -k8,8V -k1,1n | awk '{OFS="\t"}{print $1,$2,"|",$3,$4,"|",$5,$6,"|",$7,"|",$8,$9}' | column -t >${ervf}.cdd.gag.RT.2ltr.site.colum-t.txt
done

mkdir -p LTRs

awk '{print $1}' ${len_ltr} | cut -d "#" -f1 | while read i; do
	less ${i}.cdd.gag.RT.2ltr.site.colum-t.txt | awk 'NR%2 ==0{print $12,$1,$2}' | sed 's/:/ /g' | sed 's/-/ /g' | awk '{print $1":"$2+$4-1"-"$2+$5-1,$1"_"$2"_"$3}' | while read j; do
		m=($j)
		samtools faidx ${genome} ${m[0]} | sed "s/${m[0]}/${m[1]}/g" >./LTRs/${i}\.${m[1]}.cdd.gag.RT_ltr2.fa
	done
	less ${i}.cdd.gag.RT.2ltr.site.colum-t.txt | awk 'NR%2 ==1{print $12,$1,$2}' | sed 's/:/ /g' | sed 's/-/ /g' | awk '{print $1":"$2+$4-1"-"$2+$5-1,$1"_"$2"_"$3}' | while read j; do
		m=($j)
		samtools faidx ${genome} ${m[0]} | sed "s/${m[0]}/${m[1]}/g" >./LTRs/${i}.${m[1]}.cdd.gag.RT_ltr1.fa
	done
done

less ${len_ltr} | while read erv len ltr; do
	ervf=$(echo ${erv} | cut -d "#" -f 1)
	less ${ervf}.cdd.gag.RT.2ltr.site.ltr-ltr.txt | while read aaa bbb lll; do
		aaa1=$(echo ${aaa} | sed "s/:/_/g;s/-/_/g")
		echo -e "blastn -subject ./LTRs/${ervf}.${aaa1}.cdd.gag.RT_ltr1.fa -query ./LTRs/${ervf}.${aaa1}.cdd.gag.RT_ltr2.fa -outfmt 6 -out ./LTRs/${ervf}.${aaa1}.cdd.gag.RT_ltr1-ltr2.out"
	done >blastn.${ervf}.ltr1_ltr2_iden.sh
done

ls blastn.rnd*.ltr1_ltr2_iden.sh | while read aaa; do
	bash ${aaa}
done

less ${len_ltr} | while read erv len ltr; do
	ervf=$(echo ${erv} | cut -d "#" -f 1)
	cat ./LTRs/${ervf}.*.cdd.gag.RT_ltr1-ltr2.out | sort -k1,1V -k4,4nr | sort -k1,1 -u >${ervf}.merge_0.5.cdd.gag.RT_ltr1-ltr2.out
done

less ${len_ltr} | while read erv len ltr; do
	ervf=$(echo ${erv} | cut -d "#" -f 1)
	less ${ervf}.merge_0.5.cdd.gag.RT_ltr1-ltr2.out | awk '{OFS="\t"}{print "'${ervf}'",$1,$5/$4/(2*'${mut_rate}')/1000000*'${gen_time}'}' >${ervf}.merge_0.5.cdd.gag.RT_ltr1-ltr2.mya.txt
done

less ADL.earlGrey.families.ERVK-full.fa.size.len_ltr.txt | sed "s/#/\t/g" | cut -f 1 | while read aaa; do cat CD-search.local.ID.merge_0.5.fa.pbs.txt | sed "s/IDIDID/${aaa}/g" >CD-search.local.${aaa}.merge_0.5.fa.pbs; done
