#for file in /data/rivasas2/limbs/reads_clean/pig/*trimmed.clipped.fastq*; do ln -s $file .;done
#for file in /data/rivasas2/limbs/alignment_clean/pig/*sorted.sam; do ln -s $file .; done

# FPKM
for feature in genes isoforms; do
for file in /data/rivasas2/limbs/cufflinks_clean/pig/cuffdiff.*/${feature}.fpkm_tracking; do
   prefix=$(echo $file | awk '{n=split($0,a,"/");split(a[n-1],b,".");print b[2]}')
   ln -s $file ${prefix}_${feature}.fpkm_tracking
done; done
# Differences
for feature in gene isoform; do
for file in /data/rivasas2/limbs/cufflinks_clean/pig/cuffdiff.*/${feature}_exp.diff; do
   prefix=$(echo $file | awk '{n=split($0,a,"/");split(a[n-1],b,".");print b[2]}')
   ln -s $file ${prefix}_${feature}_exp.diff
done; done

rm summary.txt
for stage in 22 26; do
	echo Pig_${stage} \
	",\`genes_FPKM <https://132.239.135.28/public/limbs/files/pig/"${stage}_genes.fpkm_tracking">\`_" \
	",\`isoforms_FPKM <https://132.239.135.28/public/limbs/files/pig/"${stage}_isoforms.fpkm_tracking">\`_" \
	",\`gene_diff <https://132.239.135.28/public/limbs/files/pig/"${stage}_gene_exp.diff">\`_" \
	",\`isoform_diff <https://132.239.135.28/public/limbs/files/pig/"${stage}_isoform_exp.diff">\`_" \
	>> summary.txt
done


#rm summary.txt
#for stage in 21 25; do
#for limb in FL HL; do
#for line in L004 L005 L006 L007; do
#for file in $(ls *sam | grep ${stage} | grep ${limb} | grep ${line} ); do
#	lane_number=$( echo $line | awk '{print substr($0,4)}')
#	name=${stage}${limb}_L${lane_number}
#	echo "\`"$name" <https://132.239.135.28/public/limbs/files/pig/"${file}">\`_" >> summary.txt
#
#done; done; done; done
#
