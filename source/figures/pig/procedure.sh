declare -A sname
sname["D20"]=20
sname["22_5"]=22
sname["21half"]=22
sname["24"]=26
sname["25half"]=26
dummy_stage="none"
specie=pig
echo "Pig" > ${specie}_results.table.rst
echo "===" >> ${specie}_results.table.rst
echo "   " >> ${specie}_results.table.rst
echo ".. csv-table:: Pig: Summary of data and gene expression results." >> ${specie}_results.table.rst
echo "   :header: Library,Alias,Reads,Alignment,Alig rate,Genes expression,Isoforms expression" >> ${specie}_results.table.rst
echo "   " >> ${specie}_results.table.rst
for limb in FL HL; do
for stage in D20 22_5 21half 24 25half; do
stageName=${sname[$stage]}
if [[ "$stageName" != "$dummy_stage" ]]; then
	rep=0
	dummy_stage=$stageName
fi
for fastq in $(ls /data/rivasas2/limbs/reads_clean/${specie}/*trimmed.clipped.fastq* | grep $limb | grep $stage) \
             $(ls /data2/rivasas2/limbs/reads_clean/new_2014/${specie}/*trimmed.clipped.fastq* | grep $limb | grep $stage);do
    rep=$(($rep+1))
    name=$(echo $fastq | awk '{n=split($0,a,"/");split(a[n],b,".trimmed");print b[1]}')
    alias=${specie}_${limb}_${stageName}_${rep}

    # Link fastq files
    #ln -s ${fastq} .
    fastqLink="\`fastq <https://132.239.135.28/public/limbs/files/${specie}/${name}.trimmed.clipped.fastq.gz>\`_"

    # Link sam files
    for file in $(ls /data/rivasas2/limbs/alignment_clean/${specie}/* \
                /data2/rivasas2/limbs/alignment_clean/new_2014/${specie}/* \
                | grep ${name}.trimmed.clipped.sorted.sam ); do
        #ln -s $file .
        folder=$(echo $file | awk -v x=$name '{split($0,a,x);print a[1]}')
        aligRate=$(grep ${name} ${folder}/mapped_reads.txt | awk '{print $3/$2}')
    done
    samLink="\`sam <https://132.239.135.28/public/limbs/files/${specie}/${name}.trimmed.clipped.sorted.sam>\`_"

    # Link cufflinks gene expression files
#    for file in $(ls /data2/rivasas2/limbs/cufflinks_time_series/${specie}/cufflinks_*/genes.fpkm_tracking \
#                     /data2/rivasas2/limbs/cufflinks_time_series/new_2014/${specie}/cufflinks_*/genes.fpkm_tracking \
#                    | grep cufflinks_${name} ); do
#       ln -s $file cufflinks_${name}.genes.fpkm_tracking
#   done
#    for file in $(ls /data2/rivasas2/limbs/cufflinks_time_series/${specie}/cufflinks_*/isoforms.fpkm_tracking \
#                     /data2/rivasas2/limbs/cufflinks_time_series/new_2014/${specie}/cufflinks_*/isoforms.fpkm_tracking \
#                    | grep cufflinks_${name} ); do
#       ln -s $file cufflinks_${name}.isoforms.fpkm_tracking
#   done
    genesCufflinkLink="\`genes_fpkm <https://132.239.135.28/public/limbs/files/${specie}/cufflinks_${name}.genes.fpkm_tracking>\`_"
	isoformsCufflinkLink="\`isoforms_fpkm <https://132.239.135.28/public/limbs/files/${specie}/cufflinks_${name}.isoforms.fpkm_tracking>\`_"

	#Table summary
	echo -e "   "$name","${alias}","${fastqLink}","${samLink}","${aligRate}","${genesCufflinkLink}","${isoformsCufflinkLink} >> ${specie}_results.table.rst

done; done;done

###########################################################################
# Cuffdiff

specie=pig
echo "Pig" > ${specie}_cuffdiff.table.rst
echo "===" >> ${specie}_cuffdiff.table.rst
echo "   " >> ${specie}_cuffdiff.table.rst
echo ".. csv-table:: Pig: Gene expression differences between fore and hind limbs." >> ${specie}_cuffdiff.table.rst
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
