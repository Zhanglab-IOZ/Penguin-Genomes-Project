### 🧬 Estimating Insertion Time of Full-Length Endogenous Retroviruses (ERVs) using Flanking LTRs: An Adelie ERVK Example

This repository contains the scripts and steps for estimating the insertion time of full-length Endogenous Retroviruses (ERVs) by comparing their flanking Long Terminal Repeats (LTRs), using **Adelie ERVK** as a case study.

---

#### 🛠️ Dependencies

The following resources and software packages are required to run the analysis:

* Manually corrected ERVK consensus sequences (`fasta` format)
* Length of the manually corrected ERVK consensus sequence and LTR length (`txt` format)
* NCBI Conserved Domain Database (CDD)
* NCBI `rmblast` package
* `KentUtils` package
* `bedtools`, `samtools`, etc.

---

#### 📜 Procedure: Step-by-Step

##### 1. Full-Length ERV Alignment and LTR Sequence Alignment (Complete Matches)


`./rmblast/rmblast.pbs`
`./rmblast/rmblast.full-len.sh`
`./rmblast/rmblast.full-len.LTR.sh`
`./rmblast/rmblast.full-len.LTR-iden-mya.sh`

##### 2. Full-Length ERV Alignment, Fragment Merging, and LTR Sequence Alignment (Fragmented Matches)


`./rmblast_bed-merge/rmblast.bed-merge.sh`
`./rmblast_bed-merge/rmblast.bed-merge.CD-search-full.sh`
`./rmblast_bed-merge/rmblast.bed-merge.CD-search.local.rnd-1_family-57f.merge_0.5.fa.pbs`
`./rmblast_bed-merge/rmblast.bed-merge.CD-search.local.rnd-5_family-10513f.merge_0.5.fa.pbs`

##### 3. Result Merging and Molecular Clock Estimation

The final step merges the results from the complete and fragmented alignments and uses the calculated LTR divergence to estimate the insertion time based on a molecular clock model.

`./rmblast.w_merge.mya.sh`

---

#### 💡 How It Works (Brief Summary)

The insertion time of a full-length ERV is typically estimated based on the assumption that its two LTRs were identical at the time of insertion (i.e., immediately after the retrotransposition event). Over evolutionary time, mutations accumulate independently in the two LTRs. The sequence divergence ($D$) between the two LTRs is calculated, and the insertion time ($T$) is estimated using the formula:

$$T = \frac{D}{2\mu}$$

where $\mu$ is the neutral substitution rate per site per year for the host species.
