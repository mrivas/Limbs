Mouse, pig, and opossum analyzes
================================

Alignment
---------

As with bat, the first step was to pre-processed all RNA-seq libraries (from mouse, pig, and opossum) to eliminate duplicates, adaptors, and low quality bases ( trimming qual<20 at the reads’ 3 prime end).

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


.. csv-table:: Mouse results.
   :header: "Lib-name","Alias","Clean Reads","Alignment","Alignment rate","Expression"

   Mus_W2_FL_1_2_CGATGT_L002_R1_001,mouse_FL_W2_1,`fastq <https://132.239.135.28/public/limbs/files/mouse/Mus_W2_FL_1_2_CGATGT_L002_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/mouse/Mus_W2_FL_1_2_CGATGT_L002_R1_001.trimmed.clipped.sorted.sam>`_,0.96979,`genes_fpkm <https://132.239.135.28/public/limbs/files/mouse/cufflinks_Mus_W2_FL_1_2_CGATGT_L002_R1_001.genes.fpkm_tracking>`_
   Mus_W2_FL3_4_ACAGTG_L003_R1_001,mouse_FL_W2_2,`fastq <https://132.239.135.28/public/limbs/files/mouse/Mus_W2_FL3_4_ACAGTG_L003_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/mouse/Mus_W2_FL3_4_ACAGTG_L003_R1_001.trimmed.clipped.sorted.sam>`_,0.986523,`genes_fpkm <https://132.239.135.28/public/limbs/files/mouse/cufflinks_Mus_W2_FL3_4_ACAGTG_L003_R1_001.genes.fpkm_tracking>`_
   Mus_W2_FL_GCCAAT_L004_R1_001,mouse_FL_W2_3,`fastq <https://132.239.135.28/public/limbs/files/mouse/Mus_W2_FL_GCCAAT_L004_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/mouse/Mus_W2_FL_GCCAAT_L004_R1_001.trimmed.clipped.sorted.sam>`_,0.974812,`genes_fpkm <https://132.239.135.28/public/limbs/files/mouse/cufflinks_Mus_W2_FL_GCCAAT_L004_R1_001.genes.fpkm_tracking>`_
   mouse_W3_4_FL2R_ACAGTG_L006_R1_001,mouse_FL_W3_4_1,`fastq <https://132.239.135.28/public/limbs/files/mouse/mouse_W3_4_FL2R_ACAGTG_L006_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/mouse/mouse_W3_4_FL2R_ACAGTG_L006_R1_001.trimmed.clipped.sorted.sam>`_,0.983271,`genes_fpkm <https://132.239.135.28/public/limbs/files/mouse/cufflinks_mouse_W3_4_FL2R_ACAGTG_L006_R1_001.genes.fpkm_tracking>`_
   Mus_W3_4_FL3_GTGAAA_L002_R1_001,mouse_FL_W3_4_2,`fastq <https://132.239.135.28/public/limbs/files/mouse/Mus_W3_4_FL3_GTGAAA_L002_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/mouse/Mus_W3_4_FL3_GTGAAA_L002_R1_001.trimmed.clipped.sorted.sam>`_,0.978826,`genes_fpkm <https://132.239.135.28/public/limbs/files/mouse/cufflinks_Mus_W3_4_FL3_GTGAAA_L002_R1_001.genes.fpkm_tracking>`_
   Mouse_W6_FL2_CGATGT_L001_R1_001,mouse_FL_W6_1,`fastq <https://132.239.135.28/public/limbs/files/mouse/Mouse_W6_FL2_CGATGT_L001_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/mouse/Mouse_W6_FL2_CGATGT_L001_R1_001.trimmed.clipped.sorted.sam>`_,0.987349,`genes_fpkm <https://132.239.135.28/public/limbs/files/mouse/cufflinks_Mouse_W6_FL2_CGATGT_L001_R1_001.genes.fpkm_tracking>`_
   Mus_W6_FL1_CGATGT_L003_R1_001,mouse_FL_W6_2,`fastq <https://132.239.135.28/public/limbs/files/mouse/Mus_W6_FL1_CGATGT_L003_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/mouse/Mus_W6_FL1_CGATGT_L003_R1_001.trimmed.clipped.sorted.sam>`_,0.986541,`genes_fpkm <https://132.239.135.28/public/limbs/files/mouse/cufflinks_Mus_W6_FL1_CGATGT_L003_R1_001.genes.fpkm_tracking>`_
   Mus_W6_FL3_ACAGTG_L004_R1_001,mouse_FL_W6_3,`fastq <https://132.239.135.28/public/limbs/files/mouse/Mus_W6_FL3_ACAGTG_L004_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/mouse/Mus_W6_FL3_ACAGTG_L004_R1_001.trimmed.clipped.sorted.sam>`_,0.9837,`genes_fpkm <https://132.239.135.28/public/limbs/files/mouse/cufflinks_Mus_W6_FL3_ACAGTG_L004_R1_001.genes.fpkm_tracking>`_
   mouse_W2_HL3R_CGATGT_L005_R1_001,mouse_HL_W2_1,`fastq <https://132.239.135.28/public/limbs/files/mouse/mouse_W2_HL3R_CGATGT_L005_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/mouse/mouse_W2_HL3R_CGATGT_L005_R1_001.trimmed.clipped.sorted.sam>`_,0.987322,`genes_fpkm <https://132.239.135.28/public/limbs/files/mouse/cufflinks_mouse_W2_HL3R_CGATGT_L005_R1_001.genes.fpkm_tracking>`_
   Mus_W2_HL1_2_CTTGTA_L005_R1_001,mouse_HL_W2_2,`fastq <https://132.239.135.28/public/limbs/files/mouse/Mus_W2_HL1_2_CTTGTA_L005_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/mouse/Mus_W2_HL1_2_CTTGTA_L005_R1_001.trimmed.clipped.sorted.sam>`_,0.983087,`genes_fpkm <https://132.239.135.28/public/limbs/files/mouse/cufflinks_Mus_W2_HL1_2_CTTGTA_L005_R1_001.genes.fpkm_tracking>`_
   mouse_W3_4_HL1R_GCCAAT_L007_R1_001,mouse_HL_W3_4_1,`fastq <https://132.239.135.28/public/limbs/files/mouse/mouse_W3_4_HL1R_GCCAAT_L007_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/mouse/mouse_W3_4_HL1R_GCCAAT_L007_R1_001.trimmed.clipped.sorted.sam>`_,0.989685,`genes_fpkm <https://132.239.135.28/public/limbs/files/mouse/cufflinks_mouse_W3_4_HL1R_GCCAAT_L007_R1_001.genes.fpkm_tracking>`_
   mouse_W3_4_HL3R_CTTGTA_L004_R1_001,mouse_HL_W3_4_2,`fastq <https://132.239.135.28/public/limbs/files/mouse/mouse_W3_4_HL3R_CTTGTA_L004_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/mouse/mouse_W3_4_HL3R_CTTGTA_L004_R1_001.trimmed.clipped.sorted.sam>`_,0.982345,`genes_fpkm <https://132.239.135.28/public/limbs/files/mouse/cufflinks_mouse_W3_4_HL3R_CTTGTA_L004_R1_001.genes.fpkm_tracking>`_
   Mus_W3_4_HL2_GCCAAT_L003_R1_001,mouse_HL_W3_4_3,`fastq <https://132.239.135.28/public/limbs/files/mouse/Mus_W3_4_HL2_GCCAAT_L003_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/mouse/Mus_W3_4_HL2_GCCAAT_L003_R1_001.trimmed.clipped.sorted.sam>`_,0.987093,`genes_fpkm <https://132.239.135.28/public/limbs/files/mouse/cufflinks_Mus_W3_4_HL2_GCCAAT_L003_R1_001.genes.fpkm_tracking>`_
   Mouse_W6_HL1_GTGAAA_L001_R1_001,mouse_HL_W6_1,`fastq <https://132.239.135.28/public/limbs/files/mouse/Mouse_W6_HL1_GTGAAA_L001_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/mouse/Mouse_W6_HL1_GTGAAA_L001_R1_001.trimmed.clipped.sorted.sam>`_,0.983112,`genes_fpkm <https://132.239.135.28/public/limbs/files/mouse/cufflinks_Mouse_W6_HL1_GTGAAA_L001_R1_001.genes.fpkm_tracking>`_
   mouse_W6_HL3_GTGAAA_L005_R1_001,mouse_HL_W6_2,`fastq <https://132.239.135.28/public/limbs/files/mouse/mouse_W6_HL3_GTGAAA_L005_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/mouse/mouse_W6_HL3_GTGAAA_L005_R1_001.trimmed.clipped.sorted.sam>`_,0.980267,`genes_fpkm <https://132.239.135.28/public/limbs/files/mouse/cufflinks_mouse_W6_HL3_GTGAAA_L005_R1_001.genes.fpkm_tracking>`_
   Mus_W6_HL2_GCCAAT_L005_R1_001,mouse_HL_W6_3,`fastq <https://132.239.135.28/public/limbs/files/mouse/Mus_W6_HL2_GCCAAT_L005_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/mouse/Mus_W6_HL2_GCCAAT_L005_R1_001.trimmed.clipped.sorted.sam>`_,0.986401,`genes_fpkm <https://132.239.135.28/public/limbs/files/mouse/cufflinks_Mus_W6_HL2_GCCAAT_L005_R1_001.genes.fpkm_tracking>`_
   Mus_W6_HL3_CTTGTA_L002_R1_001,mouse_HL_W6_4,`fastq <https://132.239.135.28/public/limbs/files/mouse/Mus_W6_HL3_CTTGTA_L002_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/mouse/Mus_W6_HL3_CTTGTA_L002_R1_001.trimmed.clipped.sorted.sam>`_,0.987156,`genes_fpkm <https://132.239.135.28/public/limbs/files/mouse/cufflinks_Mus_W6_HL3_CTTGTA_L002_R1_001.genes.fpkm_tracking>`_

.. csv-table:: Opossum results.
   :header: "Lib-name","Alias","Clean Reads","Alignment","Alignment rate","Expression"

   opossumSt27FL1_CGATGT_L005_R1_001,opossum_FL_27_1,`fastq <https://132.239.135.28/public/limbs/files/opossum/opossumSt27FL1_CGATGT_L005_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/opossumSt27FL1_CGATGT_L005_R1_001.trimmed.clipped.sorted.sam>`_,0.931492,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_opossumSt27FL1_CGATGT_L005_R1_001.genes.fpkm_tracking>`_
   opossumSt27FL2_ACAGTG_L003_R1_001,opossum_FL_27_2,`fastq <https://132.239.135.28/public/limbs/files/opossum/opossumSt27FL2_ACAGTG_L003_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/opossumSt27FL2_ACAGTG_L003_R1_001.trimmed.clipped.sorted.sam>`_,0.930409,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_opossumSt27FL2_ACAGTG_L003_R1_001.genes.fpkm_tracking>`_
   opossumSt27FL3_GCCAAT_L002_R1_001,opossum_FL_27_3,`fastq <https://132.239.135.28/public/limbs/files/opossum/opossumSt27FL3_GCCAAT_L002_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/opossumSt27FL3_GCCAAT_L002_R1_001.trimmed.clipped.sorted.sam>`_,0.913644,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_opossumSt27FL3_GCCAAT_L002_R1_001.genes.fpkm_tracking>`_
   Mono_St28_FL1_CGATGT_L004_R1_001,opossum_FL_28_1,`fastq <https://132.239.135.28/public/limbs/files/opossum/Mono_St28_FL1_CGATGT_L004_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/Mono_St28_FL1_CGATGT_L004_R1_001.trimmed.clipped.sorted.sam>`_,0.909821,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_Mono_St28_FL1_CGATGT_L004_R1_001.genes.fpkm_tracking>`_
   Mono_St28_FL2_ACAGTG_L005_R1_001,opossum_FL_28_2,`fastq <https://132.239.135.28/public/limbs/files/opossum/Mono_St28_FL2_ACAGTG_L005_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/Mono_St28_FL2_ACAGTG_L005_R1_001.trimmed.clipped.sorted.sam>`_,0.908274,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_Mono_St28_FL2_ACAGTG_L005_R1_001.genes.fpkm_tracking>`_
   Mono_St28_FL3_GCCAAT_L002_R1_001,opossum_FL_28_3,`fastq <https://132.239.135.28/public/limbs/files/opossum/Mono_St28_FL3_GCCAAT_L002_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/Mono_St28_FL3_GCCAAT_L002_R1_001.trimmed.clipped.sorted.sam>`_,0.870133,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_Mono_St28_FL3_GCCAAT_L002_R1_001.genes.fpkm_tracking>`_
   opossum_St28_FL1_CGATGT_L006_R1_001,opossum_FL_28_4,`fastq <https://132.239.135.28/public/limbs/files/opossum/opossum_St28_FL1_CGATGT_L006_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/opossum_St28_FL1_CGATGT_L006_R1_001.trimmed.clipped.sorted.sam>`_,0.908772,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_opossum_St28_FL1_CGATGT_L006_R1_001.genes.fpkm_tracking>`_
   opossum_St28_FL2_ACAGTG_L007_R1_001,opossum_FL_28_5,`fastq <https://132.239.135.28/public/limbs/files/opossum/opossum_St28_FL2_ACAGTG_L007_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/opossum_St28_FL2_ACAGTG_L007_R1_001.trimmed.clipped.sorted.sam>`_,0.913338,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_opossum_St28_FL2_ACAGTG_L007_R1_001.genes.fpkm_tracking>`_
   opossum_St28_FL3_GCCAAT_L004_R1_001,opossum_FL_28_6,`fastq <https://132.239.135.28/public/limbs/files/opossum/opossum_St28_FL3_GCCAAT_L004_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/opossum_St28_FL3_GCCAAT_L004_R1_001.trimmed.clipped.sorted.sam>`_,0.866846,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_opossum_St28_FL3_GCCAAT_L004_R1_001.genes.fpkm_tracking>`_
   St29_Control_FL_index19_GTGAAA_L001_R1_001,opossum_FL_29_1,`fastq <https://132.239.135.28/public/limbs/files/opossum/St29_Control_FL_index19_GTGAAA_L001_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/St29_Control_FL_index19_GTGAAA_L001_R1_001.trimmed.clipped.sorted.sam>`_,0.920155,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_St29_Control_FL_index19_GTGAAA_L001_R1_001.genes.fpkm_tracking>`_
   St29_Control_FL_index2_CGATGT_L001_R1_001,opossum_FL_29_2,`fastq <https://132.239.135.28/public/limbs/files/opossum/St29_Control_FL_index2_CGATGT_L001_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/St29_Control_FL_index2_CGATGT_L001_R1_001.trimmed.clipped.sorted.sam>`_,0.902657,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_St29_Control_FL_index2_CGATGT_L001_R1_001.genes.fpkm_tracking>`_
   St29_Control_FL_index7_CAGATC_L001_R1_001,opossum_FL_29_3,`fastq <https://132.239.135.28/public/limbs/files/opossum/St29_Control_FL_index7_CAGATC_L001_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/St29_Control_FL_index7_CAGATC_L001_R1_001.trimmed.clipped.sorted.sam>`_,0.922021,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_St29_Control_FL_index7_CAGATC_L001_R1_001.genes.fpkm_tracking>`_
   St30FLControl1_GCCAAT_L004_R1_001,opossum_FL_30_1,`fastq <https://132.239.135.28/public/limbs/files/opossum/St30FLControl1_GCCAAT_L004_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/St30FLControl1_GCCAAT_L004_R1_001.trimmed.clipped.sorted.sam>`_,0.92152,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_St30FLControl1_GCCAAT_L004_R1_001.genes.fpkm_tracking>`_
   St30FLControl2_GCCAAT_L005_R1_001,opossum_FL_30_2,`fastq <https://132.239.135.28/public/limbs/files/opossum/St30FLControl2_GCCAAT_L005_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/St30FLControl2_GCCAAT_L005_R1_001.trimmed.clipped.sorted.sam>`_,0.870099,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_St30FLControl2_GCCAAT_L005_R1_001.genes.fpkm_tracking>`_
   St30FLControl3_GCCAAT_L006_R1_001,opossum_FL_30_3,`fastq <https://132.239.135.28/public/limbs/files/opossum/St30FLControl3_GCCAAT_L006_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/St30FLControl3_GCCAAT_L006_R1_001.trimmed.clipped.sorted.sam>`_,0.929081,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_St30FLControl3_GCCAAT_L006_R1_001.genes.fpkm_tracking>`_
   St31FLControl1_GTGAAA_L004_R1_001,opossum_FL_31_1,`fastq <https://132.239.135.28/public/limbs/files/opossum/St31FLControl1_GTGAAA_L004_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/St31FLControl1_GTGAAA_L004_R1_001.trimmed.clipped.sorted.sam>`_,0.91526,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_St31FLControl1_GTGAAA_L004_R1_001.genes.fpkm_tracking>`_
   St31FLControl2_GTGAAA_L005_R1_001,opossum_FL_31_2,`fastq <https://132.239.135.28/public/limbs/files/opossum/St31FLControl2_GTGAAA_L005_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/St31FLControl2_GTGAAA_L005_R1_001.trimmed.clipped.sorted.sam>`_,0.94274,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_St31FLControl2_GTGAAA_L005_R1_001.genes.fpkm_tracking>`_
   St31FLControl3_GTGAAA_L006_R1_001,opossum_FL_31_3,`fastq <https://132.239.135.28/public/limbs/files/opossum/St31FLControl3_GTGAAA_L006_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/St31FLControl3_GTGAAA_L006_R1_001.trimmed.clipped.sorted.sam>`_,0.923467,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_St31FLControl3_GTGAAA_L006_R1_001.genes.fpkm_tracking>`_
   Mono_St30_HL3_4_GTGAAA_L003_R1_001,opossum_HL_30_1,`fastq <https://132.239.135.28/public/limbs/files/opossum/Mono_St30_HL3_4_GTGAAA_L003_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/Mono_St30_HL3_4_GTGAAA_L003_R1_001.trimmed.clipped.sorted.sam>`_,0.918085,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_Mono_St30_HL3_4_GTGAAA_L003_R1_001.genes.fpkm_tracking>`_
   opossum_St30_HL1R_CTTGTA_L005_R1_001,opossum_HL_30_2,`fastq <https://132.239.135.28/public/limbs/files/opossum/opossum_St30_HL1R_CTTGTA_L005_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/opossum_St30_HL1R_CTTGTA_L005_R1_001.trimmed.clipped.sorted.sam>`_,0.917669,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_opossum_St30_HL1R_CTTGTA_L005_R1_001.genes.fpkm_tracking>`_
   opossum_St30_HL3R_GTGAAA_L006_R1_001,opossum_HL_30_3,`fastq <https://132.239.135.28/public/limbs/files/opossum/opossum_St30_HL3R_GTGAAA_L006_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/opossum_St30_HL3R_GTGAAA_L006_R1_001.trimmed.clipped.sorted.sam>`_,0.910475,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_opossum_St30_HL3R_GTGAAA_L006_R1_001.genes.fpkm_tracking>`_
   Mono_St31_HL1_CGATGT_L005_R1_001,opossum_HL_31_1,`fastq <https://132.239.135.28/public/limbs/files/opossum/Mono_St31_HL1_CGATGT_L005_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/Mono_St31_HL1_CGATGT_L005_R1_001.trimmed.clipped.sorted.sam>`_,0.920572,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_Mono_St31_HL1_CGATGT_L005_R1_001.genes.fpkm_tracking>`_
   Mono_St31_HL2_ACAGTG_L002_R1_001,opossum_HL_31_2,`fastq <https://132.239.135.28/public/limbs/files/opossum/Mono_St31_HL2_ACAGTG_L002_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/Mono_St31_HL2_ACAGTG_L002_R1_001.trimmed.clipped.sorted.sam>`_,0.918843,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_Mono_St31_HL2_ACAGTG_L002_R1_001.genes.fpkm_tracking>`_
   St31_Control_HL_index19_GTGAAA_L003_R1_001,opossum_HL_31_3,`fastq <https://132.239.135.28/public/limbs/files/opossum/St31_Control_HL_index19_GTGAAA_L003_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/St31_Control_HL_index19_GTGAAA_L003_R1_001.trimmed.clipped.sorted.sam>`_,0.923268,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_St31_Control_HL_index19_GTGAAA_L003_R1_001.genes.fpkm_tracking>`_
   St31_Control_HL_index2_CGATGT_L003_R1_001,opossum_HL_31_4,`fastq <https://132.239.135.28/public/limbs/files/opossum/St31_Control_HL_index2_CGATGT_L003_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/St31_Control_HL_index2_CGATGT_L003_R1_001.trimmed.clipped.sorted.sam>`_,0.920553,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_St31_Control_HL_index2_CGATGT_L003_R1_001.genes.fpkm_tracking>`_
   St31_Control_HL_index7_CAGATC_L003_R1_001,opossum_HL_31_5,`fastq <https://132.239.135.28/public/limbs/files/opossum/St31_Control_HL_index7_CAGATC_L003_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/St31_Control_HL_index7_CAGATC_L003_R1_001.trimmed.clipped.sorted.sam>`_,0.921626,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_St31_Control_HL_index7_CAGATC_L003_R1_001.genes.fpkm_tracking>`_
   Mono_St_32_HL1_CTTGTA_L004_R1_001,opossum_HL_32_1,`fastq <https://132.239.135.28/public/limbs/files/opossum/Mono_St_32_HL1_CTTGTA_L004_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/Mono_St_32_HL1_CTTGTA_L004_R1_001.trimmed.clipped.sorted.sam>`_,0.918143,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_Mono_St_32_HL1_CTTGTA_L004_R1_001.genes.fpkm_tracking>`_
   Mono_St32_HL2_GTGAAA_L005_R1_001,opossum_HL_32_2,`fastq <https://132.239.135.28/public/limbs/files/opossum/Mono_St32_HL2_GTGAAA_L005_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/Mono_St32_HL2_GTGAAA_L005_R1_001.trimmed.clipped.sorted.sam>`_,0.913137,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_Mono_St32_HL2_GTGAAA_L005_R1_001.genes.fpkm_tracking>`_
   opossum_St32_HL3R_GTGAAA_L007_R1_001,opossum_HL_32_3,`fastq <https://132.239.135.28/public/limbs/files/opossum/opossum_St32_HL3R_GTGAAA_L007_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/opossum/opossum_St32_HL3R_GTGAAA_L007_R1_001.trimmed.clipped.sorted.sam>`_,0.923339,`genes_fpkm <https://132.239.135.28/public/limbs/files/opossum/cufflinks_opossum_St32_HL3R_GTGAAA_L007_R1_001.genes.fpkm_tracking>`_

.. csv-table:: Pig results.
   :header: "Lib-name","Alias","Clean Reads","Alignment","Alignment rate","Expression"

   pigD20FL1_CGATGT_L002_R1_001,pig_FL_20_1,`fastq <https://132.239.135.28/public/limbs/files/pig/pigD20FL1_CGATGT_L002_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/pig/pigD20FL1_CGATGT_L002_R1_001.trimmed.clipped.sorted.sam>`_,0.914928,`genes_fpkm <https://132.239.135.28/public/limbs/files/pig/cufflinks_pigD20FL1_CGATGT_L002_R1_001.genes.fpkm_tracking>`_
   pigD22_5_23_5FL1_CGATGT_L003_R1_001,pig_FL_22_1,`fastq <https://132.239.135.28/public/limbs/files/pig/pigD22_5_23_5FL1_CGATGT_L003_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/pig/pigD22_5_23_5FL1_CGATGT_L003_R1_001.trimmed.clipped.sorted.sam>`_,0.924877,`genes_fpkm <https://132.239.135.28/public/limbs/files/pig/cufflinks_pigD22_5_23_5FL1_CGATGT_L003_R1_001.genes.fpkm_tracking>`_
   pigD22_5_23_5FL2_ACAGTG_L004_R1_001,pig_FL_22_2,`fastq <https://132.239.135.28/public/limbs/files/pig/pigD22_5_23_5FL2_ACAGTG_L004_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/pig/pigD22_5_23_5FL2_ACAGTG_L004_R1_001.trimmed.clipped.sorted.sam>`_,0.924981,`genes_fpkm <https://132.239.135.28/public/limbs/files/pig/cufflinks_pigD22_5_23_5FL2_ACAGTG_L004_R1_001.genes.fpkm_tracking>`_
   pigD24FL1_CGATGT_L004_R1_001,pig_FL_26_1,`fastq <https://132.239.135.28/public/limbs/files/pig/pigD24FL1_CGATGT_L004_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/pig/pigD24FL1_CGATGT_L004_R1_001.trimmed.clipped.sorted.sam>`_,0.917696,`genes_fpkm <https://132.239.135.28/public/limbs/files/pig/cufflinks_pigD24FL1_CGATGT_L004_R1_001.genes.fpkm_tracking>`_
   pigD24FL2_ACAGTG_L002_R1_001,pig_FL_26_2,`fastq <https://132.239.135.28/public/limbs/files/pig/pigD24FL2_ACAGTG_L002_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/pig/pigD24FL2_ACAGTG_L002_R1_001.trimmed.clipped.sorted.sam>`_,0.91807,`genes_fpkm <https://132.239.135.28/public/limbs/files/pig/cufflinks_pigD24FL2_ACAGTG_L002_R1_001.genes.fpkm_tracking>`_
   pigD24FL3_GCCAAT_L003_R1_001,pig_FL_26_3,`fastq <https://132.239.135.28/public/limbs/files/pig/pigD24FL3_GCCAAT_L003_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/pig/pigD24FL3_GCCAAT_L003_R1_001.trimmed.clipped.sorted.sam>`_,0.911923,`genes_fpkm <https://132.239.135.28/public/limbs/files/pig/cufflinks_pigD24FL3_GCCAAT_L003_R1_001.genes.fpkm_tracking>`_
   pig_21half_22half_FL_CGATGT_L007_R1_001,pig_FL_22_1,`fastq <https://132.239.135.28/public/limbs/files/pig/pig_21half_22half_FL_CGATGT_L007_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/pig/pig_21half_22half_FL_CGATGT_L007_R1_001.trimmed.clipped.sorted.sam>`_,0.91307,`genes_fpkm <https://132.239.135.28/public/limbs/files/pig/cufflinks_pig_21half_22half_FL_CGATGT_L007_R1_001.genes.fpkm_tracking>`_
   pig_25half_26half_FL_GCCAAT_L005_R1_001,pig_FL_26_1,`fastq <https://132.239.135.28/public/limbs/files/pig/pig_25half_26half_FL_GCCAAT_L005_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/pig/pig_25half_26half_FL_GCCAAT_L005_R1_001.trimmed.clipped.sorted.sam>`_,0.906529,`genes_fpkm <https://132.239.135.28/public/limbs/files/pig/cufflinks_pig_25half_26half_FL_GCCAAT_L005_R1_001.genes.fpkm_tracking>`_
   pigD20HL1_GCCAAT_L004_R1_001,pig_HL_20_1,`fastq <https://132.239.135.28/public/limbs/files/pig/pigD20HL1_GCCAAT_L004_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/pig/pigD20HL1_GCCAAT_L004_R1_001.trimmed.clipped.sorted.sam>`_,0.87476,`genes_fpkm <https://132.239.135.28/public/limbs/files/pig/cufflinks_pigD20HL1_GCCAAT_L004_R1_001.genes.fpkm_tracking>`_
   pigD20HL2_CTTGTA_L005_R1_001,pig_HL_20_2,`fastq <https://132.239.135.28/public/limbs/files/pig/pigD20HL2_CTTGTA_L005_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/pig/pigD20HL2_CTTGTA_L005_R1_001.trimmed.clipped.sorted.sam>`_,0.918537,`genes_fpkm <https://132.239.135.28/public/limbs/files/pig/cufflinks_pigD20HL2_CTTGTA_L005_R1_001.genes.fpkm_tracking>`_
   pigD22_5_23_5HL1_GCCAAT_L005_R1_001,pig_HL_22_1,`fastq <https://132.239.135.28/public/limbs/files/pig/pigD22_5_23_5HL1_GCCAAT_L005_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/pig/pigD22_5_23_5HL1_GCCAAT_L005_R1_001.trimmed.clipped.sorted.sam>`_,0.913773,`genes_fpkm <https://132.239.135.28/public/limbs/files/pig/cufflinks_pigD22_5_23_5HL1_GCCAAT_L005_R1_001.genes.fpkm_tracking>`_
   pigD22_5_23_5HL2_CTTGTA_L002_R1_001,pig_HL_22_2,`fastq <https://132.239.135.28/public/limbs/files/pig/pigD22_5_23_5HL2_CTTGTA_L002_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/pig/pigD22_5_23_5HL2_CTTGTA_L002_R1_001.trimmed.clipped.sorted.sam>`_,0.930277,`genes_fpkm <https://132.239.135.28/public/limbs/files/pig/cufflinks_pigD22_5_23_5HL2_CTTGTA_L002_R1_001.genes.fpkm_tracking>`_
   pigD24HL1_CTTGTA_L004_R1_001,pig_HL_26_1,`fastq <https://132.239.135.28/public/limbs/files/pig/pigD24HL1_CTTGTA_L004_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/pig/pigD24HL1_CTTGTA_L004_R1_001.trimmed.clipped.sorted.sam>`_,0.919365,`genes_fpkm <https://132.239.135.28/public/limbs/files/pig/cufflinks_pigD24HL1_CTTGTA_L004_R1_001.genes.fpkm_tracking>`_
   pigD24HL2_CTTGTA_L003_R1_001,pig_HL_26_2,`fastq <https://132.239.135.28/public/limbs/files/pig/pigD24HL2_CTTGTA_L003_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/pig/pigD24HL2_CTTGTA_L003_R1_001.trimmed.clipped.sorted.sam>`_,0.920233,`genes_fpkm <https://132.239.135.28/public/limbs/files/pig/cufflinks_pigD24HL2_CTTGTA_L003_R1_001.genes.fpkm_tracking>`_
   pigD24HL3_ACAGTG_L005_R1_001,pig_HL_26_3,`fastq <https://132.239.135.28/public/limbs/files/pig/pigD24HL3_ACAGTG_L005_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/pig/pigD24HL3_ACAGTG_L005_R1_001.trimmed.clipped.sorted.sam>`_,0.917207,`genes_fpkm <https://132.239.135.28/public/limbs/files/pig/cufflinks_pigD24HL3_ACAGTG_L005_R1_001.genes.fpkm_tracking>`_
   pig_21half_22half_HL_ACAGTG_L004_R1_001,pig_HL_22_1,`fastq <https://132.239.135.28/public/limbs/files/pig/pig_21half_22half_HL_ACAGTG_L004_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/pig/pig_21half_22half_HL_ACAGTG_L004_R1_001.trimmed.clipped.sorted.sam>`_,0.911548,`genes_fpkm <https://132.239.135.28/public/limbs/files/pig/cufflinks_pig_21half_22half_HL_ACAGTG_L004_R1_001.genes.fpkm_tracking>`_
   pig_25half_26half_HL_CTTGTA_L006_R1_001,pig_HL_26_1,`fastq <https://132.239.135.28/public/limbs/files/pig/pig_25half_26half_HL_CTTGTA_L006_R1_001.trimmed.clipped.fastq.gz>`_,`sam <https://132.239.135.28/public/limbs/files/pig/pig_25half_26half_HL_CTTGTA_L006_R1_001.trimmed.clipped.sorted.sam>`_,0.919551,`genes_fpkm <https://132.239.135.28/public/limbs/files/pig/cufflinks_pig_25half_26half_HL_CTTGTA_L006_R1_001.genes.fpkm_tracking>`_

.. note::

   Library names follows the convention: ${stage}${limb}_${line}. For instance, mouse library of stage W2 on forward limb, sequenced on line 2 is named: W2FL_L2

Alignment

.. csv-table:: Reference genome and annotations
   :header: "Mouse (Mus_musculus.GRCm38.73)", "Opossum (Monodelphis_domestica.BROADO5.73)", "Pig (Sus_scrofa.Sscrofa10.2.73)"
   :widths: 4, 4, 4
   
   `Genome  <ftp://ftp.ensembl.org/pub/release-73/fasta/mus_musculus/dna/Mus_musculus.GRCm38.73.dna_sm.toplevel.fa.gz>`_, `Genome <ftp://ftp.ensembl.org/pub/release-73/fasta/monodelphis_domestica/dna/Monodelphis_domestica.BROADO5.73.dna_sm.toplevel.fa.gz>`_, `Genome <ftp://ftp.ensembl.org/pub/release-73/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa10.2.73.dna_sm.toplevel.fa.gz>`_
   `Annotations <ftp://ftp.ensembl.org/pub/release-73/gtf/mus_musculus/Mus_musculus.GRCm38.73.gtf.gz>`_, `Annotations <ftp://ftp.ensembl.org/pub/release-73/gtf/monodelphis_domestica/Monodelphis_domestica.BROADO5.73.gtf.gz>`_, `Annotations <ftp://ftp.ensembl.org/pub/release-73/gtf/sus_scrofa/Sus_scrofa.Sscrofa10.2.73.gtf.gz>`_



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


Gene expression and fore vs hind limbs differences
--------------------------------------------------

To improve the statistical inference of differential expression we used all replicates


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

.. csv-table:: Expression and difference files
   :header: "Stage", "Gene exp", "Isoforms exp", "Gene diff", "Isoform exp"

   Mouse_W2 ,`genes_FPKM <https://132.239.135.28/public/limbs/files/mouse/W2_genes.fpkm_tracking>`_ ,`isoforms_FPKM <https://132.239.135.28/public/limbs/files/mouse/W2_isoforms.fpkm_tracking>`_ ,`gene_diff <https://132.239.135.28/public/limbs/files/mouse/W2_gene_exp.diff>`_ ,`isoform_diff <https://132.239.135.28/public/limbs/files/mouse/W2_isoform_exp.diff>`_
   Mouse_W3_4 ,`genes_FPKM <https://132.239.135.28/public/limbs/files/mouse/W3_4_genes.fpkm_tracking>`_ ,`isoforms_FPKM <https://132.239.135.28/public/limbs/files/mouse/W3_4_isoforms.fpkm_tracking>`_ ,`gene_diff <https://132.239.135.28/public/limbs/files/mouse/W3_4_gene_exp.diff>`_ ,`isoform_diff <https://132.239.135.28/public/limbs/files/mouse/W3_4_isoform_exp.diff>`_
   Mouse_W6 ,`genes_FPKM <https://132.239.135.28/public/limbs/files/mouse/W6_genes.fpkm_tracking>`_ ,`isoforms_FPKM <https://132.239.135.28/public/limbs/files/mouse/W6_isoforms.fpkm_tracking>`_ ,`gene_diff <https://132.239.135.28/public/limbs/files/mouse/W6_gene_exp.diff>`_ ,`isoform_diff <https://132.239.135.28/public/limbs/files/mouse/W6_isoform_exp.diff>`_
   Opossum_30 ,`genes_FPKM <https://132.239.135.28/public/limbs/files/opossum/30_genes.fpkm_tracking>`_ ,`isoforms_FPKM <https://132.239.135.28/public/limbs/files/opossum/30_isoforms.fpkm_tracking>`_ ,`gene_diff <https://132.239.135.28/public/limbs/files/opossum/30_gene_exp.diff>`_ ,`isoform_diff <https://132.239.135.28/public/limbs/files/opossum/30_isoform_exp.diff>`_
   Opossum_31 ,`genes_FPKM <https://132.239.135.28/public/limbs/files/opossum/31_genes.fpkm_tracking>`_ ,`isoforms_FPKM <https://132.239.135.28/public/limbs/files/opossum/31_isoforms.fpkm_tracking>`_ ,`gene_diff <https://132.239.135.28/public/limbs/files/opossum/31_gene_exp.diff>`_ ,`isoform_diff <https://132.239.135.28/public/limbs/files/opossum/31_isoform_exp.diff>`_
   Pig_22 ,`genes_FPKM <https://132.239.135.28/public/limbs/files/pig/22_genes.fpkm_tracking>`_ ,`isoforms_FPKM <https://132.239.135.28/public/limbs/files/pig/22_isoforms.fpkm_tracking>`_ ,`gene_diff <https://132.239.135.28/public/limbs/files/pig/22_gene_exp.diff>`_ ,`isoform_diff <https://132.239.135.28/public/limbs/files/pig/22_isoform_exp.diff>`_
   Pig_26 ,`genes_FPKM <https://132.239.135.28/public/limbs/files/pig/26_genes.fpkm_tracking>`_ ,`isoforms_FPKM <https://132.239.135.28/public/limbs/files/pig/26_isoforms.fpkm_tracking>`_ ,`gene_diff <https://132.239.135.28/public/limbs/files/pig/26_gene_exp.diff>`_ ,`isoform_diff <https://132.239.135.28/public/limbs/files/pig/26_isoform_exp.diff>`_


Between limbs comparisons at diferent stages in opossum
-------------------------------------------------------

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

.. csv-table:: Expression and difference files
   :header: "Stage", "Gene exp", "Isoforms exp", "Gene diff", "Isoform exp"

   Opossum_28FL_31HL ,`genes_FPKM <https://132.239.135.28/public/limbs/files/opossum/28FL_31HL_genes.fpkm_tracking>`_ ,`isoforms_FPKM <https://132.239.135.28/public/limbs/files/opossum/28FL_31HL_isoforms.fpkm_tracking>`_ ,`gene_diff <https://132.239.135.28/public/limbs/files/opossum/28FL_31HL_gene_exp.diff>`_ ,`isoform_diff <https://132.239.135.28/public/limbs/files/opossum/28FL_31HL_isoform_exp.diff>`_
   Opossum_29FL_32HL ,`genes_FPKM <https://132.239.135.28/public/limbs/files/opossum/29FL_32HL_genes.fpkm_tracking>`_ ,`isoforms_FPKM <https://132.239.135.28/public/limbs/files/opossum/29FL_32HL_isoforms.fpkm_tracking>`_ ,`gene_diff <https://132.239.135.28/public/limbs/files/opossum/29FL_32HL_gene_exp.diff>`_ ,`isoform_diff <https://132.239.135.28/public/limbs/files/opossum/29FL_32HL_isoform_exp.diff>`_




