.. _betweenSpeciesNewData:

Between species analyzes NEW DATA
=================================

In this chapter, we analyzed the gene expression differences between species across development. The pipeline is divided into two big tasks. The first, is to determine the orthologs genes between all four species. The second, to determine the set of genes with the most divergent/conserved expression profiles among species.

Determine orthologous genes
---------------------------

First, we got the DNA sequences (fasta format) of mouse, opossum, and pig from `ENSEMBL <https://www.ensembl.org/biomart>`_, using Perl scrips: 
`mouse.pl <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/mouse.pl>`_, 
`opossum.pl <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/opossum.pl>`_, and 
`pig.pl <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/pig.pl>`_

.. sidebar:: Output

   | `mouse.genes.fasta <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/mouse.genes.fasta>`_ 
   | `opossum.genes.fasta <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/opossum.genes.fasta>`_ 
   | `pig.genes.fasta <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/pig.genes.fasta>`_

.. code-block:: bash

   #############################
   # Get FASTA sequence of genes
   #############################

   for specie in mouse opossum pig; do
      perl $specie.pl > $specie.genes.fasta
   done
           

Then, using bat's transcriptome as reference 
(`Trinity.subset.fasta <https://132.239.135.28/public/limbs/files/bat/Trinity.subset.fasta>`_)
we aligned (blastn, E-value 1e-20), one by one, the transcriptomes of mouse, opossum, and pig
. 

.. sidebar:: Output

   | `mouse.blastn.outfmt6 <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/mouse.blastn.outfmt6>`_ 
   | `opossum.blastn.outfmt6 <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/opossum.blastn.outfmt6>`_ 
   | `pig.blastn.outfmt6 <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/pig.blastn.outfmt6>`_ 

.. code-block:: bash

   #######################################
   # Blast BAT genes against other species
   #######################################
   
   for specie in mouse opossum pig; do
      echo makdeblastdb $species =====================
      makeblastdb -in $specie.genes.fasta -dbtype nucl
      echo blastn bat-$specie ------------------------
      blastn \
          -query Trinity.subset.fasta \
          -db $specie.genes.fasta \
          -out $specie.blastn.outfmt6 \
          -evalue 1e-20 \
          -num_threads 10 \
          -max_target_seqs 1 \
          -outfmt 6
   done



To reduce ambiguity, we filter out bat genes matching more than one orthologs gene in any other species. We extracted the gene-id columns (1 and 2) of the output files of blastn,

.. sidebar:: Output

   `bat_mouse.orthologs.txt <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/bat_mouse.orthologs.txt>`_ 
   `bat_opossum.orthologs.txt <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/bat_opossum.orthologs.txt>`_ 
   `bat_pig..orthologs.txt <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/bat_pig.orthologs.txt>`_ 

.. code-block:: bash
   
   ############################################
   # Filter out the BAT genes with 
   # isoforms matching different ortholog genes
   ############################################
   
   for specie in mouse opossum pig; do
   
      echo $specie ============================   
      
      awk 'BEGIN{FS=OFS="\t"} {
	     split($1,a,"_seq")
		 print a[1],$2}' ${specie}.blastn.outfmt6 \ 
		 | sort -u -k1,1 \
		 | awk 'BEGIN{FS=OFS="\t"} {
		    split($2,a,"|") 
			print $1,a[1]}' > bat_${specie}.orthologs.txt
   done


removed duplicates,

.. sidebar:: Output

   `bat_mouse.orthologs.uniq.txt <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/bat_mouse.orthologs.uniq.txt>`_ 
   `bat_opossum.orthologs.uniq.txt <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/bat_opossum.orthologs.uniq.txt>`_ 
   `bat_pig..orthologs.uniq.txt <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/bat_pig.orthologs.uniq.txt>`_ 

.. code-block:: bash

   #######################################
   # Remove non-unique bat-species matches
   #######################################

   for specie in mouse opossum pig; do
      sort -k2 bat_${species}.orthologs.txt \
	  | uniq -f1 -u \
	  > bat_${species}.orthologs.uniq.txt
   done


and finally, we used a python script 
`(getAllOrthologs2.py)  <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/getAllOrthologs2.py>`_
to determine the bat genes that have orthologs sequences in all the other three species.

.. sidebar:: Output

   `all.orthologs.uniq.txt <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/all.orthologs.uniq.txt>`_ 

.. code-block:: bash

   ###############################################################
   # Find bat genes with orthologous genes in all the other species
   ################################################################

   python getAllOrthologs2.py > all.orthologs.uniq.txt



Conservation of ortholog's gene expression across species
---------------------------------------------------------

Conservation metric
*******************

We analyzed how conserved are the gene expression profiles of bat, mouse, opossum, and pig across embryonic limb development on both: hind and fore limbs. The equivalent developmental stages between species is presented in Table1. 

.. sidebar:: Table 1: Equivalent stages between species. 

   \*In opossum fore and hind limbs have different stage correspondances, equivalente stages are presetned as: fore-limb / hind-limb.

.. _table_correspondance:

   ====== ===== === ======== ===
   Stage  Mouse Bat Opossum* Pig
   ====== ===== === ======== ===
   ridge  W2    13  27 / 30  20 
   bud    W3_4  14  28 / 31  22
   paddle W6    15  29 / 32  26
   ====== ===== === ======== ===

|
|
|
|

To quantify conservation at each development stage, we use the mean of all species pairwise Spearman coefficients ( :math:`c` ):

.. math:: 
   
   c = \frac{ 1 }{ \binom{n}{k} } \sum_{i=1}^{k-1} \sum_{j>i}^{k} r_{i,j} 

Where :math:`r_{i,j}` is the Spearman coefficient between species :math:`i` and :math:`j` at a given stage, and :math:`k` is the total number of species under study in a particular experimental setting (3 for setting 1 and 2, and 4 for setting 3). We selected Spearman rather than Pearson coefficient, as the former is robust against out-layers. 

Results
*******

In what follows, all analyzes were done using `R <http://www.r-project.org/>`_, a free, open-source, data-analysis software. The R scripts for hind and fore limbs can be download here:

| `analysis.all.FL.R <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/analysis.all.FL.R>`_
| `analysis.all.HL.R <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/analysis.all.HL.R>`_
|

Normalization and scaling of the RNA-seq libraries
..................................................

The gene expression values of each RNA-seq library were normalized by gene length and library size. By each stage gene expression were scaled using the `DESeq <http://genomebiology.com/2010/11/10/r106>`_ method ( see `Figure 1`_).

.. _Figure 1:

.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/raw.allFL.svg
   :width: 40 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/raw.allHL.svg
   :width: 40 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/normalized.allFL.svg
   :width: 40 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/normalized.allHL.svg
   :width: 40 %
Figure 1: Normalization and scaling of RNA-seq libraries. Gene expression was normalized by gene length and library size (first row). Then libraries were scaled to reduce between library variance using DESeq method (second row). Columns 1, 2, and 3 correspond to experiment settings 1, 2, and 3, respectively. FL: fore limbs; HL: hind limbs; b\_: bat; m\_: mouse; o\_: opposum; p\_: pig

Gene expression profiles
........................

Based on 6,583 orthologous (that are  common to all 4 species), we use equation (1) to measure between-species gene expression  conservation at each developmental stage. All pairwise Spearman coefficients are presented  on Figure 2, where it can be seen that all species are positively correlated ( all Spearman  coefficients were above 0.5 ). On experimental settings 1 and 2, the conservation level  (Figure 3A, and 3B) decreases from stages 3-4 to 6, and the same tendency is observed on  experimental setting 3, where there is a constant decrease in conservation from stage 2 to 6  (Figure 3C).   To estimate how robust this conservation measurements are against different sets or  orthologous genes, we created gene subsamples with sizes ranging from 50 to 100% of all  orthologous genes (Figure 3). For each sample size, we measured conservation on 500 sets  of genes randomly selected. Based on the resulting distributions (shown as boxplots on  Figure 3) it can be observed that for experimental setting 3, only 70% of the genes are  necessary to find a statistically significant difference between conservation levels 2 and 3-4  (not overlapping 95% confidence intervals; Figure 3C). However on all experimental settings,  the conservation differences between stages 3-4 and 6 were highly dependent on the chosen  genes. In this cases, using 90% or less orthologous genes produces not significant  conservation differences between stages (the 95% confidence intervals overlap in all cases;  Figures 3A, 3B, and 3C).  

.. _Figure 2:

.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/heatmap2.allFL.svg
   :width: 30 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/heatmap34.allFL.svg
   :width: 30 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/heatmap6.allFL.svg
   :width: 30 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/heatmap2.allHL.svg
   :width: 30 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/heatmap34.allHL.svg
   :width: 30 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/heatmap6.allHL.svg
   :width: 30 % 
Figure 2: Pairwise Spearman coefficients values. Rows a, b, and c correspond to   experimental settings 1, 2, and 3. FL: fore limbs; HL: hind limbs; b\_: bat; m\_: mouse; o\_:   opposum; p\_: pig 

.. _Figure 3:

.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/divergency.allFL.svg
   :width: 30 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/divergency.allHL.svg
   :width: 30 % 
Figure 3: Conservation of gene expression between species. Figures A, B, and C correspond to experimental settings 1, 2, and 3 respectively. FL: fore limbs; HL: hind limbs.

Discussion
**********

We found that gene expression conservation between bat, mouse, and opossum decreases   from stage 2 to 6. This trend is the opposite of the one observed for morphological   conservation at equivalent limb developmental stages ( species: mouse, opossum, pig, and   horse; Ross et al., 2013). This lack of correlation between genetic and morphological   conservation patterns may spring from the fact that we only used orthologous genes in our   analysis. We hypothesized that between-species morphological differences may be driven by   genes that have been under the influence of divergent selective pressures on bat, mouse,   and opossum and therefore unlikely to be orthologous


List of divergent and conserved genes
*************************************

List of those genes that are conserved (or are different) among species at the different stages

I've cluster the orthologous genes according to their cross-species standard deviation (SD). High cross-species SD means that a gene is divergent across species, and low SD that the gene is conserved.

I'm attaching lists with the lowest SD (lower than 5% and 25% quantiles ) and greatest SD (greater than 75% and %95 quantiles) for all three scenarios: FL, HL, and noPig (HL but without pig species).  

The first columns of each list are the gene ID of each species, followed by their expression values, and with the last column containing the cross-species SD. The distributions figures of the SD for each scenario are presented in Figure 4.


.. toctree::

   conserved_divergent.table.rst


.. _Figure 4:

.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/sd_distributions_allFL.svg
   :width: 70 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/sd_distributions_allHL.svg
   :width: 70 % 
Figure 4: Distribution of standard deviation. Rows 1 and 2 correspond to fore and hind limbs, respectively. FL: fore limbs; HL: hind limbs.


Phylogenetic analysis
=====================


:cite:`Ross2013` reported that the morphology disparity of forelimbs among mammals (mice, opossums, horses, and pigs) decreases from ridge to bud stage, but start increasing again from paddle stage. A possible gene expression explanation for this hourglass-like morphological development, is that during early stages of limb development the transcriptome is dominated by species-specific genes, then control moves towards genes that are common to all mammals, and finally goes back to species-specific genes.

To test this hypothesis, for each species we computed the evolutionary age of their genes. Young genes are the ones that have evolved recently and therefore are more species-specific. Conversely, old genes have been inherited from common ancestors and therefore are common among species. We computed the evolutionary age of each gene by determining their most distant phylogenetic node among bat, pig, mouse, and opossum containing a detectable (blast hit, E-value 1e-5) homologue :cite:`Grosse2012`. Then, for each stage, :math:`s`, we calculated the transcriptome age index (TAI) as the sum of each gene age weighted by its expression :cite:`Domazet2010`: 

.. math::

   TAI_s = \sum_i{ a_i \left( \frac { e_{is}  }{ \sum_i e_{is} } \right) }

where, :math:`a_i` and :math:`e_{is}` are the age and expression value of gene :math:`i`.

To analyze the level of specialization of the genome, for each specie we computed their evolution age. This procedure, called phylostratography,   

.. _Figure 5:

.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/orthology_heatmap.svg
   :width: 45 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/phylogenetic_tree.svg
   :width: 45 % 
Figure 5: 

.. _Figure 6:

.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/tai.FL-HL.svg
   :width: 90 % 
Figure 6: TAI values for fore and hind limbs.   

.. _Figure 7:

.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/bat_vtai.FL.svg
   :width: 30 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/pig_vtai.FL.svg
   :width: 30 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/mouse_vtai.FL.svg
   :width: 30 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/bat_vtai.HL.svg
   :width: 30 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/pig_vtai.HL.svg
   :width: 30 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/mouse_vtai.HL.svg
   :width: 30 % 
Figure 7: Fore and Hind limbs distribution of VTAI.
