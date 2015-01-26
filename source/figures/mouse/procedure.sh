#################################################################
# Fastq clean, SAM files, and cufflinks (genes|isoforms) files

# Table summary header
specie="mouse"
echo "Mouse" > ${specie}_results.table.rst
echo "=====" >> ${specie}_results.table.rst
echo "   " >> ${specie}_results.table.rst 
echo ".. csv-table:: Mouse: Summary of data and gene expression results." >> ${specie}_results.table.rst
echo "   :header: Library,Alias,Reads,Alignment,Alig rate,Genes expression,Isoforms expression" >> ${specie}_results.table.rst
echo "   " >> ${specie}_results.table.rst
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
	#ln -s /data2/rivasas2/limbs/cufflinks_time_series/${specie}/cufflinks_${name}/isoforms.fpkm_tracking cufflinks_${name}.isoforms.fpkm_tracking
    genesCufflinkLink="\`genes_fpkm <https://132.239.135.28/public/limbs/files/${specie}/cufflinks_${name}.genes.fpkm_tracking>\`_"
    isoformsCufflinkLink="\`isoforms_fpkm <https://132.239.135.28/public/limbs/files/${specie}/cufflinks_${name}.isoforms.fpkm_tracking>\`_"

    #Table summary
    echo -e "   "$name","${alias}","${fastqLink}","${samLink}","${aligRate}","${genesCufflinkLink}","${isoformsCufflinkLink} >> ${specie}_results.table.rst

done; done;done

########################################################################3
# Cuffdiff files


specie=mouse
echo "Mouse" > ${specie}_cuffdiff.table.rst
echo "=====" >> ${specie}_cuffdiff.table.rst
echo "   " >> ${specie}_cuffdiff.table.rst 
echo ".. csv-table:: Mouse: Gene expression differences between fore and hind limbs." >> ${specie}_cuffdiff.table.rst
echo "   :header: Stage,Genes expr,Genes diff, Isoforms expr,Isoforms diff" >> ${specie}_cuffdiff.table.rst
echo "   " >> ${specie}_cuffdiff.table.rst
for folder in $(ls -d /data/rivasas2/limbs/cufflinks_clean/${specie}/cuffdiff*); do
   stage=$(echo $folder | awk '{n=split($0,a,"/");split(a[n],b,".");print b[2]}')
#   ln -s ${folder}/genes.fpkm_tracking ${stage}_genes.fpkm_tracking
#   ln -s ${folder}/isoforms.fpkm_tracking ${stage}_isoforms.fpkm_tracking
#   ln -s ${folder}/gene_exp.diff ${stage}_genes_exp.diff
#   ln -s ${folder}/isoform_exp.diff ${stage}_isoform_exp.diff
   gfl="\`Gene expression <https://132.239.135.28/public/limbs/files/${specie}/${stage}_genes.fpkm_tracking>\`_"
   gdl="\`Isoform expression  <https://132.239.135.28/public/limbs/files/${specie}/${stage}_isoforms.fpkm_tracking>\`_"
   ifl="\`Gene differences <https://132.239.135.28/public/limbs/files/${specie}/${stage}_genes_exp.diff>\`_"
   idl="\`Isoform differences <https://132.239.135.28/public/limbs/files/${specie}/${stage}_isoform_exp.diff>\`_"
   echo -e "   "${stage}","${gfl}","${gdl}","${ifl}","${idl} >> ${specie}_cuffdiff.table.rst
done


