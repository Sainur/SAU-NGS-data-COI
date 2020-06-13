### Analysis of COI data  
###### COI sequencing data were processed using USEARCH (Version 11)
#### Step-1: Merge forward and reverse sequences
`usearch -fastq_mergepairs data/*_R1.fastq -fastqout data/merged.fq  -relabel @`

#### Step-2: Remove primers
`usearch -fastx_truncate data/merged.fq -stripleft 26 -stripright 26 -fastqout data/stripped.fq`

#### Step-3: Quality filtering
`usearch -fastq_filter data/stripped.fq -fastq_trunclen 300 -fastq_maxee 1.0 \
 -fastaout final_result/Quality_filtering/reads_300.fasta`

#### Step-4:   Dereplication
`usearch --fastx_uniques final_result/Quality_filtering/reads_300.fasta -fastaout final_result/derep_300.fasta --sizeout` 

#### Step-5: Abundance sort and discard singletons 
   

`usearch -sortbysize final_result/derep_300.fasta -fastaout final_result/derep_sortbysize_300.fasta -minsize 2`

#### Step-6: Uclust clustering (Cluster sequences using UCLUST)
`usearch -cluster_fast final_result/derep_sortbysize_300.fasta -centroids final_result/cluster_otu_300.fasta \
 --id 0.97 -uc final_result/cluster_otu_300.uc -relabel Otu`
 
#### Step-7: Make  index file from database (SINTAX) 
`usearch -tax_stats Ref_database/MIDORI_LONGEST_20180221_COI.fasta -log MIDORI_LONGEST_20180221_COI_stats.txt
usearch -makeudb_usearch Ref_database/MIDORI_LONGEST_20180221_COI.fasta -output Ref_database/MIDORI_LONGEST_20180221_COI.udb`  

#### Step-8:Taxonomy
`usearch -sintax final_result/cluster_otu_300.fasta -db Ref_database/MIDORI_LONGEST_20180221_COI.udb \
 -tabbedout final_result/taxonomy_300_reads.sintax -strand both -sintax_cutoff 0.5`
 
#### Step-9: Make OTU table
`usearch -otutab final_result/Quality_filtering/reads_300.fasta -otus final_result/cluster_otu_300.fasta \
 -otutabout final_result/otutable_300.txt -biomout final_result/otutable_300.json -dbmatched -sizeout`
  

-------------------------------------------------------------------------------------------------------------------------------
 
###### Extra infomation 	 

```{r merged, echo=T}

Merging
  Fwd data/10_R1.fastq
  Rev data/10_R2.fastq
  Relabel reads as 10.#

00:00 2.4Mb  FASTQ base 33 for file data/10_R1.fastq
00:03 70Mb    100.0% 74.8% merged 
                                 
Merging
  Fwd data/11a_R1.fastq
  Rev data/11a_R2.fastq
  Relabel reads as 11a.#

00:06 71Mb    100.0% 73.0% merged
                                 
Merging
  Fwd data/12_R1.fastq
  Rev data/12_R2.fastq
  Relabel reads as 12.#

00:08 71Mb    100.0% 73.3% merged
                                 
Merging
  Fwd data/12b_R1.fastq
  Rev data/12b_R2.fastq
  Relabel reads as 12b.#

00:10 71Mb    100.0% 74.2% merged
                                 
Merging
  Fwd data/3_R1.fastq
  Rev data/3_R2.fastq
  Relabel reads as 3.#

00:12 42Mb    100.0% 73.6% merged
                                 
Merging
  Fwd data/5_R1.fastq
  Rev data/5_R2.fastq
  Relabel reads as 5.#

00:15 72Mb    100.0% 73.2% merged
                                 
Merging
  Fwd data/6b_R1.fastq
  Rev data/6b_R2.fastq
  Relabel reads as 6b.#

00:18 72Mb    100.0% 72.8% merged
                                 
Merging
  Fwd data/B10_R1.fastq
  Rev data/B10_R2.fastq
  Relabel reads as B10.#

00:19 28Mb    100.0% 72.8% merged
                                 
Merging
  Fwd data/B11_R1.fastq
  Rev data/B11_R2.fastq
  Relabel reads as B11.#

00:20 37Mb    100.0% 73.2% merged
                                 
Merging
  Fwd data/B12_R1.fastq
  Rev data/B12_R2.fastq
  Relabel reads as B12.#

00:21 34Mb    100.0% 73.4% merged
                                 
Merging
  Fwd data/B13_R1.fastq
  Rev data/B13_R2.fastq
  Relabel reads as B13.#

00:22 42Mb    100.0% 73.8% merged
                                 
Merging
  Fwd data/B14_R1.fastq
  Rev data/B14_R2.fastq
  Relabel reads as B14.#

00:24 74Mb    100.0% 74.4% merged
                                 
Merging
  Fwd data/B15_R1.fastq
  Rev data/B15_R2.fastq
  Relabel reads as B15.#

00:26 72Mb    100.0% 75.1% merged
                                 
Merging
  Fwd data/B16_R1.fastq
  Rev data/B16_R2.fastq
  Relabel reads as B16.#

00:28 74Mb    100.0% 75.8% merged
                                 
Merging
  Fwd data/B17_R1.fastq
  Rev data/B17_R2.fastq
  Relabel reads as B17.#

00:30 65Mb    100.0% 76.2% merged
                                 
Merging
  Fwd data/B18_R1.fastq
  Rev data/B18_R2.fastq
  Relabel reads as B18.#

00:31 43Mb    100.0% 76.4% merged
                                 
Merging
  Fwd data/B19_R1.fastq
  Rev data/B19_R2.fastq
  Relabel reads as B19.#

00:32 43Mb    100.0% 76.6% merged
                                 
Merging
  Fwd data/B1_R1.fastq
  Rev data/B1_R2.fastq
  Relabel reads as B1.#

00:33 56Mb    100.0% 76.7% merged
                                 
Merging
  Fwd data/B20_R1.fastq
  Rev data/B20_R2.fastq
  Relabel reads as B20.#

00:36 60Mb    100.0% 77.0% merged
                                 
Merging
  Fwd data/B3_R1.fastq
  Rev data/B3_R2.fastq
  Relabel reads as B3.#

00:37 64Mb    100.0% 77.1% merged
                                 
Merging
  Fwd data/B5_R1.fastq
  Rev data/B5_R2.fastq
  Relabel reads as B5.#

00:38 45Mb    100.0% 77.2% merged
                                 
Merging
  Fwd data/B6_R1.fastq
  Rev data/B6_R2.fastq
  Relabel reads as B6.#

00:42 77Mb    100.0% 77.0% merged
                                 
Merging
  Fwd data/B7_R1.fastq
  Rev data/B7_R2.fastq
  Relabel reads as B7.#

00:42 30Mb    100.0% 77.1% merged
                                 
Merging
  Fwd data/B8_R1.fastq
  Rev data/B8_R2.fastq
  Relabel reads as B8.#

00:45 78Mb    100.0% 77.2% merged
                                 
Merging
  Fwd data/B9_R1.fastq
  Rev data/B9_R2.fastq
  Relabel reads as B9.#

00:49 78Mb    100.0% 77.5% merged
                                 
Merging
  Fwd data/B_R1.fastq
  Rev data/B_R2.fastq
  Relabel reads as B.#

00:53 78Mb    100.0% 77.3% merged
                                 
Merging
  Fwd data/Ex1_R1.fastq
  Rev data/Ex1_R2.fastq
  Relabel reads as Ex1.#

00:55 78Mb    100.0% 77.1% merged
                                 
Merging
  Fwd data/Ex2_R1.fastq
  Rev data/Ex2_R2.fastq
  Relabel reads as Ex2.#

00:58 79Mb    100.0% 76.6% merged
                                 
Merging
  Fwd data/J_R1.fastq
  Rev data/J_R2.fastq
  Relabel reads as J.#

01:01 79Mb    100.0% 76.5% merged

Totals:
   1670570  Pairs (1.7M)
   1277674  Merged (1.3M, 76.48%)
    315824  Alignments with zero diffs (18.91%)
    331626  Too many diffs (> 5) (19.85%)
    292285  Fwd tails Q <= 2 trimmed (17.50%)
    491328  Rev tails Q <= 2 trimmed (29.41%)
     26518  Fwd too short (< 64) after tail trimming (1.59%)
      4713  Rev too short (< 64) after tail trimming (0.28%)
     30039  No alignment found (1.80%)
         0  Alignment too short (< 16) (0.00%)
     15251  Staggered pairs (0.91%) merged & trimmed
    213.38  Mean alignment length
    339.03  Mean merged length
      2.03  Mean fwd expected errors
      6.90  Mean rev expected errors
      0.09  Mean merged expected errors  """


-------------------------------------------------------------------------------------------------------------------------------
 
 Quality filtering  
00:00 1.8Mb  FASTQ base 33 for file data/stripped.fq
00:40 2.3Mb   100.0% Filtering, 89.7% passed
   1277644  Reads (1.3M)                    
    125273  Discarded reads length < 300
      6935  Discarded reads with expected errs > 1.00
   1145436  Filtered reads (1.1M, 89.7%)
-------------------------------------------------------------------------------------------------------------------------------
 
	Dreplication
	00:06 439Mb   100.0% Reading final_result/Quality_filtering/reads_300.fasta
00:10 457Mb   100.0% DF                                                    
00:10 468Mb  1145436 seqs, 244759 uniques, 194111 singletons (79.3%)
00:10 468Mb  Min size 1, median 1, max 43754, avg 4.68
00:16 454Mb   100.0% Writing final_result/derep_300.fasta
-------------------------------------------------------------------------------------------------------------------------------
 
	Clustering
	00:00 36Mb    100.0% Reading final_result/derep_sortbysize_300.fasta
00:00 22Mb    100.0% DF                                             
00:00 23Mb   50648 seqs, 50648 uniques, 50648 singletons (100.0%)
00:00 23Mb   Min size 1, median 1, max 1, avg 1.00
00:00 25Mb    100.0% DB
00:02 33Mb    100.0% 1791 clusters, max size 7018, avg 28.3
00:03 33Mb    100.0% Writing centroids to final_result/cluster_otu_300.fasta
                                                                            
      Seqs  50648 (50.6k)
  Clusters  1791
  Max size  7018
  Avg size  28.3
  Min size  1
Singletons  754, 1.5% of seqs, 42.1% of clusters
   Max mem  36Mb
      Time  3.00s
Throughput  16.9k seqs/sec.
```
