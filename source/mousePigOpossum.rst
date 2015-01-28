.. _mousePigOpossum:

Mouse, pig, and opossum analyzes
================================

Data cleaning and alignment
---------------------------

As with bat, the first step was to pre-processed all RNA-seq libraries (from mouse, pig, and opossum) to eliminate adaptors, and low quality bases ( trimming qual<20 at the reads’ 3 prime end).

.. code-block:: bash

   # Illumina adaptors
   # from page 6 on:
   # http://supportres.illumina.com/documents/myillumina/6378de81-c0cc-47d0-9281-724878bb1c30/2012-09-18_illuminacustomersequenceletter.pdf
   declare -A file_adaptor
   file_adaptor=(
   [ACAGTG]=GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG
   [CGATGT]=GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG
   [CTTGTA]=GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG
   [GCCAAT]=GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG
   [GTGAAA]=GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACGATCTCGTATGCCGTCTTCTGCTTG
   [CAGATC]=GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG
   )
   
   for file in ${reads_raw}/*.fastq; do
      length=$((${#reads_raw}+1))
      file=${file:$length:(-6)}
      # TRIMMS LOW QUALITIES
      echo "Trimming "${file}" =================================================="
      fastq_quality_trimmer -Q 33 -v -t 20 -l 20 -i ${reads_raw}/${file}.fastq -o ${file}.trimmed.fastq
      # CLIPPS SEQUENCING ADAPTORS
      sequence=$(awk -v x=$file 'END{n=split(x,a,"_");print a[n-3]}' /dev/null)
      adaptor=${file_adaptor[$sequence]}
      echo "Clipping: "${file}" =================================================="
      fastx_clipper -Q 33 -a ${adaptor} -l 20 -n -v -i ${file}.trimmed.fastq -o ${file}.trimmed.clipped.fastq
   done
   cd $base

Once cleaned, the libraries were aligned agains their corresponding reference genomes (Table 1).

.. csv-table:: Table 1: Reference genome and annotations
   :header: Mouse, Opossum, Pig
   
   `Mus_musculus.GRCm38.73.fa.gz  <ftp://ftp.ensembl.org/pub/release-73/fasta/mus_musculus/dna/Mus_musculus.GRCm38.73.dna_sm.toplevel.fa.gz>`_, `Monodelphis_domestica.BROADO5.73.fa.gz <ftp://ftp.ensembl.org/pub/release-73/fasta/monodelphis_domestica/dna/Monodelphis_domestica.BROADO5.73.dna_sm.toplevel.fa.gz>`_, `Sus_scrofa.Sscrofa10.2.73.fa.gz <ftp://ftp.ensembl.org/pub/release-73/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa10.2.73.dna_sm.toplevel.fa.gz>`_
   `Mus_musculus.GRCm38.73.gtf.gz <ftp://ftp.ensembl.org/pub/release-73/gtf/mus_musculus/Mus_musculus.GRCm38.73.gtf.gz>`_, `Monodelphis_domestica.BROADO5.73.gtf.gz <ftp://ftp.ensembl.org/pub/release-73/gtf/monodelphis_domestica/Monodelphis_domestica.BROADO5.73.gtf.gz>`_, `Sus_scrofa.Sscrofa10.2.73.gtf.gz <ftp://ftp.ensembl.org/pub/release-73/gtf/sus_scrofa/Sus_scrofa.Sscrofa10.2.73.gtf.gz>`_

.. code-block:: bash

   $pathToStarDir/STAR \
       --runMode genomeGenerate \
       --genomeDir $path_to_genome_index_star \
       --genomeFastaFiles $path_to_genome_fa \
       --runThreadN 8 \
       --sjdbGTFfile $gtf

   cd $alignments
   
   for pathFile in ${reads_clean}/*trimmed.clipped.fastq
   do
      length=$((${#reads_clean}+1))
      oFile=${pathFile:$length:(-6)}
      echo $oFile
   
      $pathToStarDir/STAR \
          --genomeDir $path_to_index \
          --readFilesIn $pathFile \
          --runThreadN 10 \
          --genomeLoad NoSharedMemory  \
          --outFilterScoreMin 0 \
          --outFilterMultimapNmax 10  \
          --outFilterMismatchNmax 3 \
          --clip3pNbases 0 \
          --clip5pNbases 0 \
          --outFileNamePrefix ${oFile}. \
          --outSAMmode Full \
          --outSAMattributes Standard \
          --outSAMstrandField intronMotif \
          --outReadsUnmapped Fastx
      java -jar /home/rivasas2/tools/picard-tools-1.52/SortSam.jar \
          INPUT=${oFile}.Aligned.out.sam \
          OUTPUT=${oFile}.sorted.sam \
          SORT_ORDER=coordinate
   done
   #
   cd $base:

We used cufflinks to compute gene expression.

.. code-block:: bash

   gtf=/home/rivasas2/tools/genomes/$species/*.gtf
   fasta=/home/rivasas2/tools/genomes/$species/*.fa
   alignments=/data2/rivasas2/limbs/alignment_clean/new_2014/$species
   cufflinks_time_series=/data2/rivasas2/limbs/cufflinks_time_series/new_2014/$species
   
   cd $cufflinks_time_series
   
   for sam in ${alignments}/*.sorted.sam; do
       length=$((${#alignments}+1))
       name=${sam:$length:(-27)}
       echo "######################################################################"
       echo Gene expression $species $name
       echo "######################################################################"
   
       cufflinks \
           -G ${gtf} \
           -b $fasta \
           -u \
           -p 7 \
           -N \
           --frag-len-mean 350 --frag-len-std-dev 100 \
           -o cufflinks_$name \
           $sam
   
   done
   
   cd -


Results of alignment and gene expression
----------------------------------------

The cleaned libraries, alignments, and expression values are presented in the following links.

.. toctree::

   mouse_results.table
   opossum_results.table
   pig_results.table

.. note::

   Library names follows the convention: species_stage_limb_replicate. For instance, mouse library of stage W2, forward limb, on the first replicate is  named: mouse_W2_FL_1.


Gene expression analyses
------------------------

Gene expression differences between fore vs hind limbs were computed at each developmental stage. For each comparison, all replicates were used at once.


.. code-block:: bash

   ##################################################################
   # 3. CUFFDIFF
   #################################################################

   # Between limbs comparisons at same stage
   cd $cufflinks_clean
   
   for stage in ${stage_all[$species]}; do
      alignment_FL=""
      alignment_HL=""
      for file in ${alignments}/*.sorted.sam; do
          if [[ $file =~ $stage && $file =~ "FL" ]]; then
              if [ "$alignment_FL" = "" ]; then
                  alignment_FL=$file
              else 
                  alignment_FL=$alignment_FL,$file
              fi
          elif [[ $file =~ $stage && $file =~ "HL" ]]; then
              if [ "$alignment_HL" = "" ]; then
                  alignment_HL=$file
              else
                  alignment_HL=$alignment_HL,$file
              fi
          fi
      done
      

   
      echo CUFFDIFF $stage ================================================
      echo FL files ------------------------------------------------
      echo $alignment_FL
      echo HL files -------------------------------------------------
      echo $alignment_HL
      cuffdiff \
          ${gtf} \
          -p 10 --frag-len-mean 350 --frag-len-std-dev 100 \
          --multi-read-correct \
          --frag-bias-correct ${refGenome} \
          -o cuffdiff.${stage} \
          -L St${stage}_FL,St${stage}_HL \
          $alignment_FL $alignment_HL
   done
   cd $base

Since in opossum fore and hind limbs have similar morphology development at different stages, we compared the equivalent stages

.. code-block:: bash

   # Between limbs comparisons at diferent stages

   ####################################################################
    Comparison: 28 FL - 31 HL
   ####################################################################
   cd $cufflinks_clean
   
   alignment_FL=""
   alignment_HL=""
   alignments=/data2/rivasas2/limbs/alignment_clean/opossum
   for file in ${alignments}/*.sorted.sam; do
      if [[ $file =~ "28" && $file =~ "FL" ]]; then
          if [ "$alignment_FL" = "" ]; then
              alignment_FL=$file
          else 
              alignment_FL=$alignment_FL,$file
          fi
      fi
   done
   alignments=/data/rivasas2/limbs/alignment_clean/opossum
   for file in ${alignments}/*.sorted.sam; do
      if [[ $file =~ "31" && $file =~ "HL" ]]; then
          if [ "$alignment_HL" = "" ]; then
              alignment_HL=$file
          else
              alignment_HL=$alignment_HL,$file
          fi
      fi
   done
   

   echo CUFFDIFF $stage ================================================
   echo FL files ------------------------------------------------
   echo $alignment_FL
   echo HL files -------------------------------------------------
   echo $alignment_HL
   cuffdiff \
      ${gtf} \
      -p 10 --frag-len-mean 350 --frag-len-std-dev 100 \
      --multi-read-correct \
      --frag-bias-correct ${refGenome} \
      -o cuffdiff.28FL_31HL \
      -L St28_FL,St31_HL \
      $alignment_FL $alignment_HL
   
   cd $base

   #####################################################################
   # Comparison: 29 FL - 32 HL
   #####################################################################
   
   cd $cufflinks_clean
   
   alignment_FL=""
   alignment_HL=""
   alignments=/data2/rivasas2/limbs/alignment_clean/opossum
   for file in ${alignments}/*.sorted.sam; do
      if [[ $file =~ "29" && $file =~ "FL" ]]; then
          if [ "$alignment_FL" = "" ]; then
              alignment_FL=$file
          else 
              alignment_FL=$alignment_FL,$file
          fi
      fi
   done
   alignments=/data2/rivasas2/limbs/alignment_clean/opossum
   for file in ${alignments}/*.sorted.sam; do
      if [[ $file =~ "32" && $file =~ "HL" ]]; then
          if [ "$alignment_HL" = "" ]; then
              alignment_HL=$file
          else
              alignment_HL=$alignment_HL,$file
          fi
      fi
   done

   
   echo CUFFDIFF $stage ================================================
   echo FL files ------------------------------------------------
   echo $alignment_FL
   echo HL files -------------------------------------------------
   echo $alignment_HL
   cuffdiff \
      ${gtf} \
      -p 10 --frag-len-mean 350 --frag-len-std-dev 100 \
      --multi-read-correct \
      --frag-bias-correct ${refGenome} \
      -o cuffdiff.29FL_32HL \
      -L St29_FL,St32_HL \
      $alignment_FL $alignment_HL
   
   cd $base


Results of differential gene expression analyses
------------------------------------------------

Results per specie in the following links.

.. toctree::

   mouse_cuffdiff.table
   opossum_cuffdiff.table
   pig_cuffdiff.table
