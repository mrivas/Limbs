# FASTQ files
#for file in /data/rivasas2/limbs/reads_clean/mouse/*trimmed.clipped.fastq*; do ln -s $file .;done
# SAM files
#for file in /data/rivasas2/limbs/alignment_clean/mouse/*sorted.sam; do ln -s $file .; done

specie="mouse"
for limb in FL HL; do
for stage in W2 W3_4 W6; do
rep=0
for fastq in $(ls /data/rivasas2/limbs/reads_clean/${specie}/*trimmed.clipped.fastq* | grep $limb | grep $stage);do
	rep=$(($rep+1))
	name=$(echo $fastq | awk '{n=split($0,a,"/");split(a[n],b,".trimmed");print b[1]}')
	alias=${specie}_${limb}_${stage}_${rep}
	
	# Link fastq files
	#ln -s ${fastq} .
	fastqLink="\`fastq <https://132.239.135.28/public/limbs/files/${specie}/${name}.trimmed.clipped.fastq.gz>\`_"
	
	# Link sam files
	# ln -s /data/rivasas2/limbs/alignment_clean/${specie}/${name}.trimmed.clipped.sorted.sam .
	samLink="\`sam <https://132.239.135.28/public/limbs/files/${specie}/${name}.trimmed.clipped.sorted.sam>\`_"
	aligRate=$(grep ${name} /data/rivasas2/limbs/alignment_clean/${specie}/mapped_reads.txt | awk '{print $3/$2}')
	
	# Link cufflinks gene expression files
	#ln -s /data2/rivasas2/limbs/cufflinks_time_series/${specie}/cufflinks_${name}/genes.fpkm_tracking cufflinks_${name}.genes.fpkm_tracking
	cufflinkLink="\`genes_fpkm <https://132.239.135.28/public/limbs/files/${specie}/cufflinks_${name}.genes.fpkm_tracking>\`_"


	# Table summary
	echo -e $name","${alias}","${fastqLink}","${samLink}","${aligRate}","${cufflinkLink}
done; done;done


# Cuffdiff files
#for feature in genes isoforms; do
#for file in /data/rivasas2/limbs/cufflinks_clean/mouse/cuffdiff.*/${feature}.fpkm_tracking; do
#   prefix=$(echo $file | awk '{n=split($0,a,"/");split(a[n-1],b,".");print b[2]}')
#   ln -s $file ${prefix}_${feature}.fpkm_tracking
#done; done
## Differences
#for feature in gene isoform; do
#for file in /data/rivasas2/limbs/cufflinks_clean/mouse/cuffdiff.*/${feature}_exp.diff; do
#   prefix=$(echo $file | awk '{n=split($0,a,"/");split(a[n-1],b,".");print b[2]}')
#   ln -s $file ${prefix}_${feature}_exp.diff
#done; done
#
#rm summary.txt
#for stage in W2 W3_4 W6; do
#	echo Mouse_${stage} \
#	",\`genes_FPKM <https://132.239.135.28/public/limbs/files/mouse/"${stage}_genes.fpkm_tracking">\`_" \
#	",\`isoforms_FPKM <https://132.239.135.28/public/limbs/files/mouse/"${stage}_isoforms.fpkm_tracking">\`_" \
#	",\`gene_diff <https://132.239.135.28/public/limbs/files/mouse/"${stage}_gene_exp.diff">\`_" \
#	",\`isoform_diff <https://132.239.135.28/public/limbs/files/mouse/"${stage}_isoform_exp.diff">\`_" \
#	>> summary.txt
#done

#rm summary.txt
#for stage in W2 W3_4 W6; do
#for limb in FL HL; do
#for line in L001 L002 L003 L004 L005 L006 L007; do
#for file in $(ls *sam | grep ${stage} | grep ${limb} | grep ${line} ); do
#	lane_number=$( echo $line | awk '{print substr($0,4)}')
#	name=${stage}${limb}_L${lane_number}
#	echo "\`"$name" <https://132.239.135.28/public/limbs/files/mouse/"${file}">\`_" >> summary.txt
#
#done; done; done; done

