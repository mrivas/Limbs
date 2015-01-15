ln -s /data/rivasas2/limbs/trinity/bat/St13_14_15.unique2.fastq .
ln -s /data2/rivasas2/limbs/trinity/bat/St13_14_15.unique.trinity/Trinity.fasta .
ln -s /data2/rivasas2/limbs/trinity/bat/funcAnnot/blastx.outfmt6 .
ln -s /data2/rivasas2/limbs/trinity/bat/St13_14_15.unique.trinity/Trinity.subset.fasta .
for file in /data/rivasas2/limbs/reads_clean/bat/*clipped.fastq; do ln -s $file . ; done
for file in /data2/rivasas2/limbs/trinity/bat/*.results; do ln -s $file .; done
for file in /data/rivasas2/limbs/alignment_clean/bat/*clipped/accepted_hits.bam; do name=$(echo $file | awk '{n=split($0,a,"/"); print a[n-1]"_accepted_hits.bam"}'); ln -s $file $name; done
for file in /data2/rivasas2/limbs/trinity/bat/*transcript.sorted.bam; do ln -s $file .; done
ln -s /data2/rivasas2/limbs/trinity/bat/St13_14_15.unique.trinity/map_file.txt .
