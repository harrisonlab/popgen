#!/bin/bash
input=/home/sobczm/popgen/renseq/input/transcriptomes/really_really_final_baits
scripts=/home/sobczm/bin/popgen/renseq

#Contig clustering and bait design for srlk genes
cd $input
usearch=/home/sobczm/bin/usearch/usearch9.0.2132_i86linux32
for fasta in onion_slrk.fasta
do
id=0.9
$usearch -cluster_fast $fasta -id $id -sort length -centroids ${fasta%.fasta}_${id}_centroids.fasta
done 

00:00 41Mb    100.0% Reading onion_slrk.fasta
00:00 149Mb   100.0% DF                      
00:00 149Mb  452 seqs, 452 uniques, 452 singletons (100.0%)
00:00 149Mb  Min size 1, median 1, max 1, avg 1.00
00:00 152Mb   100.0% DB
00:00 152Mb  Sort length... done.
00:04 164Mb   100.0% 175 clusters, max size 11, avg 2.6
00:04 164Mb   100.0% Writing centroids to onion_slrk_0.9_centroids.fasta
                                                                        
      Seqs  452
  Clusters  175
  Max size  11
  Avg size  2.6
  Min size  1
Singletons  74, 16.4% of seqs, 42.3% of clusters
   Max mem  166Mb
      Time  4.00s
Throughput  113.0 seqs/sec.


                                                                        
      Seqs  452
  Clusters  36
  Max size  75
  Avg size  12.6
  Min size  1
Singletons  4, 0.9% of seqs, 11.1% of clusters
   Max mem  163Mb
      Time  2.00s
Throughput  226.0 seqs/sec.

#RLK genes
for fasta in onion.RLK.LRR.fasta
do
id=0.5
$usearch -cluster_fast $fasta -id $id -sort length -centroids ${fasta%.fasta}_${id}_centroids.fasta
done 

00:00 45Mb    100.0% Reading onion.RLK.LRR.fasta
00:00 161Mb   100.0% DF                         
00:00 153Mb  1686 seqs, 1684 uniques, 1683 singletons (99.9%)
00:00 153Mb  Min size 1, median 1, max 3, avg 1.00
00:00 156Mb   100.0% DB
00:00 156Mb  Sort length... done.
00:25 170Mb   100.0% 196 clusters, max size 69, avg 8.6
00:25 170Mb   100.0% Writing centroids to onion.RLK.LRR_0.9_centroids.fasta
                                                                           
      Seqs  1684
  Clusters  196
  Max size  69
  Avg size  8.6
  Min size  1
Singletons  29, 1.7% of seqs, 14.8% of clusters
   Max mem  174Mb
      Time  26.0s
Throughput  64.8 seqs/sec.


for fasta in onion.RLK.LysM.fasta
do
id=0.5
$usearch -cluster_fast $fasta -id $id -sort length -centroids ${fasta%.fasta}_${id}_centroids.fasta
done 

00:00 40Mb    100.0% Reading onion.RLK.LysM.fasta
00:00 157Mb   100.0% DF                          
00:00 148Mb  57 seqs, 57 uniques, 57 singletons (100.0%)
00:00 148Mb  Min size 1, median 1, max 1, avg 1.00
00:00 151Mb   100.0% DB
00:00 151Mb  Sort length... done.
00:00 160Mb   100.0% 6 clusters, max size 21, avg 9.5
00:00 160Mb   100.0% Writing centroids to onion.RLK.LysM_0.9_centroids.fasta
                                                                            
      Seqs  57
  Clusters  6
  Max size  21
  Avg size  9.5
  Min size  1
Singletons  1, 1.8% of seqs, 16.7% of clusters
   Max mem  160Mb
      Time  1.00s
Throughput  57.0 seqs/sec.

#RLP genes
for fasta in onion.RLP.LRR.fasta
do
id=0.5
$usearch -cluster_fast $fasta -id $id -sort length -centroids ${fasta%.fasta}_${id}_centroids.fasta
done 

00:00 46Mb    100.0% Reading onion.RLP.LRR.fasta
00:00 155Mb   100.0% DF                         
00:00 155Mb  3265 seqs, 3264 uniques, 3263 singletons (100.0%)
00:00 155Mb  Min size 1, median 1, max 2, avg 1.00
00:00 158Mb   100.0% DB
00:00 158Mb  Sort length... done.
00:44 175Mb   100.0% 616 clusters, max size 59, avg 5.3
00:44 175Mb   100.0% Writing centroids to onion.RLP.LRR_0.5_centroids.fasta
                                                                           
      Seqs  3264
  Clusters  616
  Max size  59
  Avg size  5.3
  Min size  1
Singletons  179, 5.5% of seqs, 29.1% of clusters
   Max mem  175Mb
      Time  44.0s
Throughput  74.2 seqs/sec.


for fasta in onion.RLP.LysM.fasta
do
id=0.5
$usearch -cluster_fast $fasta -id $id -sort length -centroids ${fasta%.fasta}_${id}_centroids.fasta
done 

00:00 40Mb    100.0% Reading onion.RLP.LysM.fasta
00:00 144Mb   100.0% DF                          
00:00 148Mb  29 seqs, 29 uniques, 29 singletons (100.0%)
00:00 148Mb  Min size 1, median 1, max 1, avg 1.00
00:00 151Mb   100.0% DB
00:00 151Mb  Sort length... done.
00:00 158Mb   100.0% 6 clusters, max size 10, avg 4.8
00:00 158Mb   100.0% Writing centroids to onion.RLP.LysM_0.5_centroids.fasta
                                                                            
      Seqs  29
  Clusters  6
  Max size  10
  Avg size  4.8
  Min size  1
Singletons  1, 3.4% of seqs, 16.7% of clusters
   Max mem  158Mb
      Time  1.00s
Throughput  29.0 seqs/sec.


#First, hard-mask repetitive sequences in each centroid file
for fasta in *centroid*.fasta
do
$usearch -fastx_mask $fasta -qmask dust -fastaout "${fasta%.*}"_masked.fasta -hardmask
done

#Design baits
for fasta in *masked.fasta
do
n=5
python $scripts/create_baits.py --inp $fasta --coverage $n \
--out ${fasta%.fasta}_baits.fasta
done

#Remove baits containing any Ns
for baits in *baits.fasta
do
python $scripts/remove_N_fasta.py $baits
done

#Subtitute "_" and ":" with "-" in sequence headers as
#required by Mycoarray.
for baits in *noN.fasta
do
sed -i 's/_/-/g' $baits 
sed -i 's/:/-/g' $baits
done

#Put all RLP baits into one file
cat onion.RLK.LRR_0.5_centroids_masked_baits_noN.fasta >> onion.RLK.all_0.5_centroids_masked_baits_noN.fasta
cat onion.RLK.LysM_0.5_centroids_masked_baits_noN.fasta >> onion.RLK.all_0.5_centroids_masked_baits_noN.fasta
#Put all RLK baits into one file
cat onion.RLP.LRR_0.5_centroids_masked_baits_noN.fasta >> onion.RLP.all_0.5_centroids_masked_baits_noN.fasta
cat onion.RLP.LysM_0.5_centroids_masked_baits_noN.fasta >> onion.RLP.all_0.5_centroids_masked_baits_noN.fasta