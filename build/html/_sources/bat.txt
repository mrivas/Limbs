.. _bat:

Bat analysis
============

*De novo* transcriptome assembly
--------------------------------

We pre-processed all *Carollia perspicillata*'s RNA-seq libraries to eliminate duplicates, adaptors, and low quality bases ( trimming qual<20 at the reads' 3 prime end).

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
   reads_raw="/path/to/raw/fastq/folder" 
   reads_clean="/path/to/outpu/folder"

   for file in ${reads_raw}/*; do
   
      # TRIMMS LOW QUALITIES
      length=$((${#reads_raw}+1))
      #file=${file:$length:(-8)}
      file=${file:$length:(-6)}
   
      echo "Trimming "${file}" =================================================="
      fastq_quality_trimmer -Q 33 -v -t 20 -l 20 -i ${reads_raw}/${file}.fastq -o ${file}.trimmed.fastq
   
      # CLIPPS SEQUENCING ADAPTORS
      sequence=$(awk -v x=$file 'END{n=split(x,a,"_");print a[n-3]}' /dev/null)
      adaptor=${file_adaptor[$sequence]}
   
      echo "Clipping: "${file} $adaptor" =================================================="
      fastx_clipper -Q 33 -a ${adaptor} -l 20 -n -v -i ${file}.trimmed.fastq -o ${file}.trimmed.clipped.fastq
       
   done


.. table:: Output cleaned RNA-seq libraries. 

   ======= =======================================================================================================================================================================
   Library Cleaned fastq files
   ======= =======================================================================================================================================================================
   13FL_L5 `bat_2011_109_13HL_ACAGTG_L005_R1_001.trimmed.clipped.fastq <https://132.239.135.28/public/limbs/files/bat/bat_2011_109_13HL_ACAGTG_L005_R1_001.trimmed.clipped.fastq>`_
   13HL_L7 `bat_2013_105_13HL_CTTGTA_L007_R1_001.trimmed.clipped.fastq <https://132.239.135.28/public/limbs/files/bat/bat_2013_105_13HL_CTTGTA_L007_R1_001.trimmed.clipped.fastq>`_
   13FL_L4 `Bat_2013_161_13FL_CGATGT_L004_R1_001.trimmed.clipped.fastq <https://132.239.135.28/public/limbs/files/bat/Bat_2013_161_13FL_CGATGT_L004_R1_001.trimmed.clipped.fastq>`_
   13HL_L6 `bat_2013_161_13HL_GCCAAT_L006_R1_001.trimmed.clipped.fastq <https://132.239.135.28/public/limbs/files/bat/bat_2013_161_13HL_GCCAAT_L006_R1_001.trimmed.clipped.fastq>`_
   13FL_L1 `St13FL162_CGATGT_L001_R1_001.trimmed.clipped.fastq <https://132.239.135.28/public/limbs/files/bat/St13FL162_CGATGT_L001_R1_001.trimmed.clipped.fastq>`_
   13FL_L2 `St13FL178_CGATGT_L002_R1_001.trimmed.clipped.fastq <https://132.239.135.28/public/limbs/files/bat/St13FL178_CGATGT_L002_R1_001.trimmed.clipped.fastq>`_
   14FL_L1 `St14FL104_ACAGTG_L001_R1_001.trimmed.clipped.fastq <https://132.239.135.28/public/limbs/files/bat/St14FL104_ACAGTG_L001_R1_001.trimmed.clipped.fastq>`_
   14FL_L2 `St14FL122_ACAGTG_L002_R1_001.trimmed.clipped.fastq <https://132.239.135.28/public/limbs/files/bat/St14FL122_ACAGTG_L002_R1_001.trimmed.clipped.fastq>`_
   14HL_L1 `St14HL104_GCCAAT_L001_R1_001.trimmed.clipped.fastq <https://132.239.135.28/public/limbs/files/bat/St14HL104_GCCAAT_L001_R1_001.trimmed.clipped.fastq>`_
   14HL_L2 `St14HL122_GCCAAT_L002_R1_001.trimmed.clipped.fastq <https://132.239.135.28/public/limbs/files/bat/St14HL122_GCCAAT_L002_R1_001.trimmed.clipped.fastq>`_
   15FL_L1 `St15FL105_CTTGTA_L001_R1_001.trimmed.clipped.fastq <https://132.239.135.28/public/limbs/files/bat/St15FL105_CTTGTA_L001_R1_001.trimmed.clipped.fastq>`_
   15FL_L2 `St15FL118_CTTGTA_L002_R1_001.trimmed.clipped.fastq <https://132.239.135.28/public/limbs/files/bat/St15FL118_CTTGTA_L002_R1_001.trimmed.clipped.fastq>`_
   15HL_L1 `St15HL105_GTGAAA_L001_R1_001.trimmed.clipped.fastq <https://132.239.135.28/public/limbs/files/bat/St15HL105_GTGAAA_L001_R1_001.trimmed.clipped.fastq>`_
   15HL_L2 `St15HL118_GTGAAA_L002_R1_001.trimmed.clipped.fastq <https://132.239.135.28/public/limbs/files/bat/St15HL118_GTGAAA_L002_R1_001.trimmed.clipped.fastq>`_
   ======= =======================================================================================================================================================================

.. note::

   Library names follows the convention: ${stage}${limb}_${line}. For instance, library of stage 13 on forward limb, sequenced on line 5 is named: 13FL_L5

At first, we attempted to align the cleaned libraries against *Myotis lucifugus*, a closely-related bat species with a readily available genome sequence and gene annotation at ENSEMBL:

| `myoLuc.fa <ftp://ftp.ensembl.org/pub/release-71/fasta/myotis_lucifugus/dna/Myotis_lucifugus.Myoluc2.0.71.dna_sm.toplevel.fa.gz>`_
| `myoLuc.gtf <ftp://ftp.ensembl.org/pub/release-71/gtf/myotis_lucifugus/Myotis_lucifugus.Myoluc2.0.71.gtf.gz>`_
|

However, the mapping (using `Tophat <http://ccb.jhu.edu/software/tophat/>`_) of libraries at stages 14 and 15 resulted in only ~9% alignment rates.


.. code-block:: bash

   # Index genome with Bowtie2
   bowtie2-build myoLuc2.fa myoLuc2

   for pathFile in *trimmed.clipped.fastq
   do
      file=${pathFile:33:(-6)}
      echo ${file}"=================================================="
      tophat -g 1 -p 10 --GTF myoLuc.gtf -o ${file} myoLuc2 ${file}.fastq
   done

.. table:: Alignment rates against *Myotis lucifugus*

   ========================================================================================================================= ==================
   Alignment file                                                                                                            Alignment rate (%)
   ========================================================================================================================= ==================
   `14FL_L1 <https://132.239.135.28/public/limbs/files/bat/St14FL104_ACAGTG_L001_R1_001.trimmed.clipped_accepted_hits.bam>`_ 8.23
   `14FL_L2 <https://132.239.135.28/public/limbs/files/bat/St14FL122_ACAGTG_L002_R1_001.trimmed.clipped_accepted_hits.bam>`_ 8.67
   `14HL_L1 <https://132.239.135.28/public/limbs/files/bat/St14HL104_GCCAAT_L001_R1_001.trimmed.clipped_accepted_hits.bam>`_ 9.11
   `14HL_L2 <https://132.239.135.28/public/limbs/files/bat/St14HL122_GCCAAT_L002_R1_001.trimmed.clipped_accepted_hits.bam>`_ 9.00
   `15FL_L1 <https://132.239.135.28/public/limbs/files/bat/St15FL105_CTTGTA_L001_R1_001.trimmed.clipped_accepted_hits.bam>`_ 7.72
   `15FL_L2 <https://132.239.135.28/public/limbs/files/bat/St15FL118_CTTGTA_L002_R1_001.trimmed.clipped_accepted_hits.bam>`_ 7.99
   `15HL_L1 <https://132.239.135.28/public/limbs/files/bat/St15HL105_GTGAAA_L001_R1_001.trimmed.clipped_accepted_hits.bam>`_ 8.69
   `15HL_L2 <https://132.239.135.28/public/limbs/files/bat/St15HL118_GTGAAA_L002_R1_001.trimmed.clipped_accepted_hits.bam>`_ 8.17
   ========================================================================================================================= ==================

As a result, we resorted to a *de novo* transcriptome assembly strategy. Here, the idea was to deduce *Carollia perspicillata* trancriptome by stitching together overlapping RNA-seq reads. To do so, we pooled all cleaned libraries into a single file, and removed duplicated sequences to reduce the computation footprint.

`St13_14_15.unique2.fastq <https://132.239.135.28/public/limbs/files/bat/St13_14_15.unique2.fastq>`_.
   
These reads were then assembled into transcripts using `Trinity <http://trinityrnaseq.sourceforge.net/>`_.

.. sidebar:: Output

   `Trinity.fasta <https://132.239.135.28/public/limbs/files/bat/Trinity.fasta>`_

.. code-block:: bash

   $ $TRINITY_HOME/Trinity.pl \
      --seqType fq --JM 80G --single St13_14_15.unique2.fastq \
      --output /data2/rivasas2/limbs/trinity/bat/St13_14_15.unique.trinity \
      --CPU 12

obtaining **350,733** gene transcripts. 

To eliminate spurious transcripts, we kept only those which sequences matched 
(`blastx <http://www.ncbi.nlm.nih.gov/books/NBK52640/>`_) 
the protein sequences of the SwissProt database 
(`uniprot_sprot.fasta <http://sourceforge.net/projects/trinotate/files/TRINOTATE_RESOURCES/TrinotateResources-20130704/uniprot_sprot.fasta.gz>`_). 

.. sidebar:: Output

   `blastx.outfmt6 <https://132.239.135.28/public/limbs/files/bat/blastx.outfmt6>`_

.. code-block:: bash

   $ blastx -query Trinity.fasta \
      -db uniprot_sprot.fasta \
      -out blastx.outfmt6 -evalue 1e-20 -num_threads 10 -max_target_seqs 1 -outfmt 6
|
|

Only **88,930** found a match. We used this information to subset the transcriptome using `samtools <http://samtools.sourceforge.net/>`_.

.. sidebar:: Output

   `Trinity.subset.fasta <https://132.239.135.28/public/limbs/files/bat/Trinity.subset.fasta>`_

.. code-block:: bash

   $ cut -f1 blastx.outfmt6 > isoforms_subset_names.txt
   $ samtools faidx Trinity.fasta
   $ xargs samtools faidx Trinity.fasta < isoforms_subset_names.txt > Trinity.subset.fasta
|
|

The final set of transcripts then was reduced to **88,930** sequences, which were used as references to map the RNA-seq libraries and compute gene expression.

Gene expression
---------------

We used `RSEM <http://deweylab.biostat.wisc.edu/rsem/>`_ to align the reads ( like `Tophat <http://ccb.jhu.edu/software/tophat/>`_, this tool internally uses `bowtie <https://http://bowtie-bio.sourceforge.net/index.shtml>` ) and to compute the expression values.

First, we used `Trinity's <http://trinityrnaseq.sourceforge.net/>`_ script to extract the gene names from our reference transcriptome.

.. sidebar:: Output

   `map_file.txt <https://132.239.135.28/public/limbs/files/bat/map_file.txt>`_

.. code-block:: bash

   extract-transcript-to-gene-map-from-trinity Trinity.subset.fasta map_file.txt
|
|
|

At this point we were ready to run `RSEM <http://deweylab.biostat.wisc.edu/rsem/>`_ on our cleaned RNA-seq libraries. We didn't remove duplicated reads from the libraries since they may be important when distinguishing differences in expression.


.. code-block:: bash

   rsem-prepare-reference \
      --transcript-to-gene-map map_file.txt \
      --quiet \
      Trinity.subset.fasta \
      reference_name.subset
   
   declare -A lanes
   lanes["13_FL"]="L001 L002 L004"
   lanes["13_HL"]="L005 L006 L007"
   lanes["14_FL"]="L001 L002"
   lanes["14_HL"]="L001 L002"
   lanes["15_FL"]="L001 L002"
   lanes["15_HL"]="L001 L002"
   
   for sample in 13 14 15; do
   for limb in FL HL; do
   for lane in ${lanes[${sample}_${limb}]}; do
      
      for file in $(ls *clipped.fastq | grep $sample | grep $limb | grep $lane ); do
          
          echo $file L${lane:3} ================================
   
          rsem-calculate-expression \
          -p 10 --fragment-length-mean 425.0 --fragment-length-sd 150.0 --quiet \
          --calc-ci \
          ${file} \
          reference_name.subset \
          ${sample}${limb}_L${lane:3}.RSEM > ${sample}${limb}_L${lane:3}.RSEM.log 2>&1
   
      done
   
   done; done; done

.. _alignment_rates_deNovo_table:

.. table:: Aligned cleaned RNAA-seq against *de novo* transcriptome, and expression at gene and transcripts levels

   ============================================================================================= ================== ============================================================================================= ================================================================================================
   Alignment file                                                                                 Alignment rate(%) Gene level expression                                                                         Transcript level expression
   ============================================================================================= ================== ============================================================================================= ================================================================================================
   `13FL_L1 <https://132.239.135.28/public/limbs/files/bat/13FL_L1.RSEM.transcript.sorted.bam>`_ 82.7               `13FL_L1_gene <https://132.239.135.28/public/limbs/files/bat/13FL_L1.RSEM.isoforms.results>`_ `13FL_L1_isoform <https://132.239.135.28/public/limbs/files/bat/13FL_L1.RSEM.isoforms.results>`_
   `13FL_L2 <https://132.239.135.28/public/limbs/files/bat/13FL_L2.RSEM.transcript.sorted.bam>`_ 83.8               `13FL_L2_gene <https://132.239.135.28/public/limbs/files/bat/13FL_L2.RSEM.isoforms.results>`_ `13FL_L2_isoform <https://132.239.135.28/public/limbs/files/bat/13FL_L2.RSEM.isoforms.results>`_
   `13FL_L4 <https://132.239.135.28/public/limbs/files/bat/13FL_L4.RSEM.transcript.sorted.bam>`_ 82.5               `13FL_L4_gene <https://132.239.135.28/public/limbs/files/bat/13FL_L4.RSEM.isoforms.results>`_ `13FL_L4_isoform <https://132.239.135.28/public/limbs/files/bat/13FL_L4.RSEM.isoforms.results>`_
   `13HL_L5 <https://132.239.135.28/public/limbs/files/bat/13HL_L5.RSEM.transcript.sorted.bam>`_ 81.9               `13FL_L4_gene <https://132.239.135.28/public/limbs/files/bat/13FL_L4.RSEM.isoforms.results>`_ `13FL_L4_isoform <https://132.239.135.28/public/limbs/files/bat/13FL_L4.RSEM.isoforms.results>`_
   `13HL_L6 <https://132.239.135.28/public/limbs/files/bat/13HL_L6.RSEM.transcript.sorted.bam>`_ 85.6               `13HL_L6_gene <https://132.239.135.28/public/limbs/files/bat/13HL_L6.RSEM.isoforms.results>`_ `13HL_L6_isoform <https://132.239.135.28/public/limbs/files/bat/13HL_L6.RSEM.isoforms.results>`_
   `13HL_L7 <https://132.239.135.28/public/limbs/files/bat/13HL_L7.RSEM.transcript.sorted.bam>`_ 80.8               `13HL_L7_gene <https://132.239.135.28/public/limbs/files/bat/13HL_L7.RSEM.isoforms.results>`_ `13HL_L7_isoform <https://132.239.135.28/public/limbs/files/bat/13HL_L7.RSEM.isoforms.results>`_
   `14FL_L1 <https://132.239.135.28/public/limbs/files/bat/14FL_L1.RSEM.transcript.sorted.bam>`_ 78.6               `14FL_L1_gene <https://132.239.135.28/public/limbs/files/bat/14FL_L1.RSEM.isoforms.results>`_ `14FL_L1_isoform <https://132.239.135.28/public/limbs/files/bat/14FL_L1.RSEM.isoforms.results>`_
   `14FL_L2 <https://132.239.135.28/public/limbs/files/bat/14FL_L2.RSEM.transcript.sorted.bam>`_ 80.5               `14FL_L2_gene <https://132.239.135.28/public/limbs/files/bat/14FL_L2.RSEM.isoforms.results>`_ `14FL_L2_isoform <https://132.239.135.28/public/limbs/files/bat/14FL_L2.RSEM.isoforms.results>`_ 
   `14HL_L1 <https://132.239.135.28/public/limbs/files/bat/14HL_L1.RSEM.transcript.sorted.bam>`_ 80.3               `14HL_L1_gene <https://132.239.135.28/public/limbs/files/bat/14HL_L1.RSEM.isoforms.results>`_ `14HL_L1_isoform <https://132.239.135.28/public/limbs/files/bat/14HL_L1.RSEM.isoforms.results>`_
   `14HL_L2 <https://132.239.135.28/public/limbs/files/bat/14HL_L2.RSEM.transcript.sorted.bam>`_ 80.1               `14HL_L2_gene <https://132.239.135.28/public/limbs/files/bat/14HL_L2.RSEM.isoforms.results>`_ `14HL_L2_isoform <https://132.239.135.28/public/limbs/files/bat/14HL_L2.RSEM.isoforms.results>`_
   `15FL_L1 <https://132.239.135.28/public/limbs/files/bat/15FL_L1.RSEM.transcript.sorted.bam>`_ 76.6               `15FL_L1_gene <https://132.239.135.28/public/limbs/files/bat/15FL_L1.RSEM.isoforms.results>`_ `15FL_L1_isoform <https://132.239.135.28/public/limbs/files/bat/15FL_L1.RSEM.isoforms.results>`_
   `15FL_L2 <https://132.239.135.28/public/limbs/files/bat/15FL_L2.RSEM.transcript.sorted.bam>`_ 77.3               `15FL_L2_gene <https://132.239.135.28/public/limbs/files/bat/15FL_L2.RSEM.isoforms.results>`_ `15FL_L2_isoform <https://132.239.135.28/public/limbs/files/bat/15FL_L2.RSEM.isoforms.results>`_
   `15HL_L1 <https://132.239.135.28/public/limbs/files/bat/15HL_L1.RSEM.transcript.sorted.bam>`_ 79.2               `15HL_L1_gene <https://132.239.135.28/public/limbs/files/bat/15HL_L1.RSEM.isoforms.results>`_ `15HL_L1_isoform <https://132.239.135.28/public/limbs/files/bat/15HL_L1.RSEM.isoforms.results>`_
   `15HL_L2 <https://132.239.135.28/public/limbs/files/bat/15HL_L2.RSEM.transcript.sorted.bam>`_ 79.2               `15HL_L2_gene <https://132.239.135.28/public/limbs/files/bat/15HL_L2.RSEM.isoforms.results>`_ `15HL_L2_isoform <https://132.239.135.28/public/limbs/files/bat/15HL_L2.RSEM.isoforms.results>`_
   ============================================================================================= ================== ============================================================================================= ================================================================================================

From the alignment rate column it can be seen that the *de novo* assembly improved the mapping in one order of magnitud.
