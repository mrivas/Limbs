#for file in /data/rivasas2/limbs/reads_clean/opossum/*trimmed.clipped.fastq*; do ln -s $file .;done
#for file in /data/rivasas2/limbs/alignment_clean/opossum/*sorted.sam; do ln -s $file .; done
#for file in /data2/rivasas2/limbs/alignment_clean/opossum/*sorted.sam; do ln -s $file .; done

## FPKM
#for feature in genes isoforms; do
#for file in /data/rivasas2/limbs/cufflinks_clean/opossum/cuffdiff.*/${feature}.fpkm_tracking; do
#   prefix=$(echo $file | awk '{n=split($0,a,"/");split(a[n-1],b,".");print b[2]}')
#   ln -s $file ${prefix}_${feature}.fpkm_tracking
#done; done
## Differences
#for feature in gene isoform; do
#for file in /data/rivasas2/limbs/cufflinks_clean/opossum/cuffdiff.*/${feature}_exp.diff; do
#   prefix=$(echo $file | awk '{n=split($0,a,"/");split(a[n-1],b,".");print b[2]}')
#   ln -s $file ${prefix}_${feature}_exp.diff
#done; done

specie="opossum"
for limb in FL HL; do
for stage in 27 28 29 30 31 32; do
rep=0
for fastq in $(ls /data/rivasas2/limbs/reads_clean/${specie}/*trimmed.clipped.fastq* | grep $limb | grep $stage) \
             $(ls /data2/rivasas2/limbs/reads_clean/new_2014/${specie}/*trimmed.clipped.fastq* | grep $limb | grep $stage);do
    rep=$(($rep+1))
    name=$(echo $fastq | awk '{n=split($0,a,"/");split(a[n],b,".trimmed");print b[1]}')
    alias=${specie}_${limb}_${stage}_${rep}

    # Link fastq files
    #ln -s ${fastq} .
    fastqLink="\`fastq <https://132.239.135.28/public/limbs/files/${specie}/${name}.trimmed.clipped.fastq.gz>\`_"

    # Link sam files
    for file in $(ls /data/rivasas2/limbs/alignment_clean/${specie}/* \
                /data2/rivasas2/limbs/alignment_clean/${specie}/* \
                /data2/rivasas2/limbs/alignment_clean/new_2014/${specie}/* \
				| grep ${name}.trimmed.clipped.sorted.sam ); do
    	#ln -s $file .
		folder=$(echo $file | awk -v x=$name '{split($0,a,x);print a[1]}')
    	aligRate=$(grep ${name} ${folder}/mapped_reads.txt | awk '{print $3/$2}')
	done
    samLink="\`sam <https://132.239.135.28/public/limbs/files/${specie}/${name}.trimmed.clipped.sorted.sam>\`_"

#    # Link cufflinks gene expression files
#    for file in $(ls /data2/rivasas2/limbs/cufflinks_time_series/${specie}/cufflinks_${name}/genes.fpkm_tracking \
#                     /data2/rivasas2/limbs/cufflinks_time_series/new_2014/${specie}/cufflinks_${name}/genes.fpkm_tracking \
#				     | grep cufflinks_${name} ); do
#		ln -s $file cufflinks_${name}.genes.fpkm_tracking
#	done
    cufflinkLink="\`genes_fpkm <https://132.239.135.28/public/limbs/files/${specie}/cufflinks_${name}.genes.fpkm_tracking>\`_"


    # Table summary
    echo -e $name","${alias}","${fastqLink}","${samLink}","${aligRate}","${cufflinkLink}
done; done;done



#rm summary.txt
#for stage in 30 31; do
#	echo Opossum_${stage} \
#	",\`genes_FPKM <https://132.239.135.28/public/limbs/files/opossum/"${stage}_genes.fpkm_tracking">\`_" \
#	",\`isoforms_FPKM <https://132.239.135.28/public/limbs/files/opossum/"${stage}_isoforms.fpkm_tracking">\`_" \
#	",\`gene_diff <https://132.239.135.28/public/limbs/files/opossum/"${stage}_gene_exp.diff">\`_" \
#	",\`isoform_diff <https://132.239.135.28/public/limbs/files/opossum/"${stage}_isoform_exp.diff">\`_" \
#	>> summary.txt
#done

# Between stages comparisons
# FPKM
#for feature in genes isoforms; do
#for file in /data2/rivasas2/limbs/cufflinks_clean/opossum/cuffdiff.*/${feature}.fpkm_tracking; do
#   prefix=$(echo $file | awk '{n=split($0,a,"/");split(a[n-1],b,".");print b[2]}')
#   ln -s $file ${prefix}_${feature}.fpkm_tracking
#done; done
## Differences
#for feature in gene isoform; do
#for file in /data2/rivasas2/limbs/cufflinks_clean/opossum/cuffdiff.*/${feature}_exp.diff; do
#   prefix=$(echo $file | awk '{n=split($0,a,"/");split(a[n-1],b,".");print b[2]}')
#   ln -s $file ${prefix}_${feature}_exp.diff
#done; done
#rm summary.txt
#for stage in 28FL_31HL 29FL_32HL; do
#	echo Opossum_${stage} \
#	",\`genes_FPKM <https://132.239.135.28/public/limbs/files/opossum/"${stage}_genes.fpkm_tracking">\`_" \
#	",\`isoforms_FPKM <https://132.239.135.28/public/limbs/files/opossum/"${stage}_isoforms.fpkm_tracking">\`_" \
#	",\`gene_diff <https://132.239.135.28/public/limbs/files/opossum/"${stage}_gene_exp.diff">\`_" \
#	",\`isoform_diff <https://132.239.135.28/public/limbs/files/opossum/"${stage}_isoform_exp.diff">\`_" \
#	>> summary.txt
#done

#rm summary.txt
#for stage in 28 29 30 31 32; do
#for limb in FL HL; do
#for line in L001 L002 L003 L004 L005 L006 L007; do
#for file in $(ls *sam | grep ${stage} | grep ${limb} | grep ${line} ); do
#	lane_number=$( echo $line | awk '{print substr($0,4)}')
#	name=${stage}${limb}_L${lane_number}
#	echo "\`"${name}" <https://132.239.135.28/public/limbs/files/opossum/"${file}">\`_" >> summary.txt
#
#done; done; done; done
#
