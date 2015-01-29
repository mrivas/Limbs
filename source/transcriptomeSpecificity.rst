.. _transcriptomeSpecificity:

Lineage-specificity of the transcriptomes across development
============================================================

Ross et al :cite:`Ross2013` reported that morphology disparity of forelimbs among mammals (mice, opossums, horses, and pigs) decreases from ridge to bud stage, but increases again from paddle stage. This hourglass shape of the compared limb development processes may be explained at the gene expression level by the species-specificity of the genes dominating each developmental stage. Here, we hypothesized that at early limb developmental stages the transcriptome is dominated by species-specific genes, then --at the bottleneck of the hourglass-- control moves towards genes that are common to all mammals, and finally goes back to species-specific genes.

To test this hypothesis, for all four species we computed the evolutionary age of each of their genes by determining their most distant phylogenetic node among bat, pig, mouse, and opossum containing a detectable (blast hit, E-value 1e-5) homologue :cite:`Quint2012`. In this frame, youth correspond to genes recently evolved and, therefore, more likely to be specific of each species. Conversely, old genes have been inherited from common ancestors and are present among all descendant species. The specificity of each species transcriptome was then quantified across limb development using the transcriptome age index (:math:`TAI`) :cite:`Domazet-Loso2010`, which is the sum of each gene evolutionary age weighted by its expression :cite:`Domazet-Loso2010`. For a given species at a given stage, :math:`s`, :math:`TAI_s` is mathematically defined as:

.. math::

   TAI_s = \sum_{i=1}^n{ ps_i \left( \frac { e_{is}  }{ \sum_{i=1}^n e_{is} } \right) }

where, :math:`ps_i` and :math:`e_{is}` are the age and expression values of gene :math:`i`, and :math:`n` is the total number of genes.

Genes evolutionary age
**********************

To determine genes homologies among any two species, we pooled the fasta sequences of all species genes to form a single database. Then, one by one  we aligned the genome of each species against this database using an E-value of 1e-5.

.. code-block:: bash

   ########################################################
   # Create Blast database with all species genes
   ########################################################
   
   echo "Create database with all species genes sequences"
   cp /data2/rivasas2/limbs/trinity/bat/St13_14_15.unique.trinity/Trinity.subset.fasta allSpecies.genes.fasta
   cat /data2/rivasas2/limbs/trinity/bat/blast/mouse.genes.fasta >> allSpecies.genes.fasta
   cat /data2/rivasas2/limbs/trinity/bat/blast/pig.genes.fasta >> allSpecies.genes.fasta
   cat /data2/rivasas2/limbs/trinity/bat/blast/opossum.genes.fasta >> allSpecies.genes.fasta

   makeblastdb -in allSpecies.genes.fasta -dbtype nucl
   
   ########################################################
   # Blast BAT genes against other species
   ########################################################

   declare -A genes
   genes["bat"]=/data2/rivasas2/limbs/trinity/bat/St13_14_15.unique.trinity/Trinity.subset.fasta
   genes["mouse"]=/data2/rivasas2/limbs/trinity/bat/blast/mouse.genes.fasta
   genes["opossum"]=/data2/rivasas2/limbs/trinity/bat/blast/opossum.genes.fasta
   genes["pig"]=/data2/rivasas2/limbs/trinity/bat/blast/pig.genes.fasta
   
   for specie in bat mouse opossum pig; do
       echo "blastn" $specie "---------------------------------------"
       blastn \
           -query ${genes[$specie]} \
           -db allSpecies.genes.fasta \
           -out ${specie}.blastn.outfmt6 \
           -evalue 1e-5 \
           -num_threads 10 \
           -max_hsps 1 \
           -outfmt 6
   done

Then, we used a python script, `phylostratum.py <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/phylostratum.py>`_,  to find among each gene alignments its most distant node on the phylogenetic tree (Figure 1, left).

.. code-block:: bash

   ##########################################################################
   # Extract most distant phylogenetic node among bat, pig, mouse and opossum
   ##########################################################################
   for specie in bat pig mouse opossum; do
      echo "Finding phylogenetic nodes for " $specie
      python phylostratum.py $specie
   done


The results, plotted in Figure 1 (left), show that bat, pig, mouse, and opossum have 4k, 11k, 21k, and 15k genes that are specific to each one of these species (they don't have and homologue in any other species). Therefore, they were assigned to phylostratum (*ps*) 1, and their evolutionary age was assigned a value of 1. On the other end, the number of genes that are common among all four species (*ps* 4) ranges from 6,902 (pig) to 11,296 (bat). These genes were assigned an evolutionary age of 4. Between these youngest and oldest genes there are *ps* 2 and *ps* 3 corresponding to the branching points of the phylogenetic tree. The ages of each set of genes was assigned concordantly to its :math:`ps`. To test if the assignment of homologue genes between any two species was coherent with the expected phylogeny, we did a hierachichal clustering based on the pair-wise number of homologue genes. The results, Figure 1 (right), show that the clustering mirrors exactly the phylogeny of this species. Thus, supporting the method to call homologue genes.  


Since *ps* 4 comprehend a set of further evolutionary ages 

.. _Figure 1:

.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/phylogenetic_tree.svg
   :width: 45 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/orthology_heatmap.svg
   :width: 45 % 
Figure 1: Phylogenetic tree and clustering of species bases of number of pair-wise orthologs genes. **On the left**, the phylogenetic tree of bat, pig, mouse, and opossum. For each of the four phylogenetic strata (*ps* 1, *ps* 2, *ps* 3, and *ps* 4) the number of genes assigned to it are presented in a color-coded manner. At each branch of the tree are presented, in black, the estimated time of evolutionary departure of any two species :cite:`Hedges2009`. **On the right**, the hierarchical clustering of species based on the number of pair-wise homologue genes. The number of homologue genes of species `j` (columns) found on species :math:`i` (rows) was normalized by the total number of genes on species :math:`j`.

Once obtained the evolutionary age for each gene on each species, we computed :math:`TAI_s` for each species across ridge, bud, and paddle stages. The results, Figure 2, show that on fore-limbs, the :math:`TAI_s` values of pig and mouse are coherent with our hypothesis. On both species, the ridge and paddle stages are dominated by young genes whereas older genes predominates on the intermediate budge stage. Opossum, on the other hand, had the opposite trend. Bat showed a path unlike the previous two types. For bat, its transcriptome tended to gain dominance of young genes as limb development progresses. 


.. _Figure 2:

.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/tai.FL-HL.svg
   :width: 90 % 
Figure 2: TAI values for fore and hind limbs.   

Since opossum, mouse, and pig were all reported to have limb developmental process whose morphological divergences resemble a hourglass shape, we expected to observe a similar trend among the three of them. To understand opossum discrepancy, we analyzed the statistical significant of each species :math:`TAI_s` trend using the procedure proposed by :cite:`Quint2012`, where the variance of :math:`TAI_s` across stages (ridge, budge, and paddle), :math:`VTAI`, is use as test statistic. The null distribution was obtained by sampling 1000 surrogates of :math:`VTAI`. Each surrogates was generated by randomly permuting the *ps* assignations. The null distribution was modeled as a gamma distribution, and it parameters estimated using the `MASS <http://cran.r-project.org/web/packages/MASS/index.html>`_ library in `R <http://www.r-project.org/>`_ (R scripts for p-value computations on `forelimbs <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/analysis.FLgeneAge.R>`_ and `hindlimbs <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/analysis.HLgeneAge.R>`_).

The statistical results, Figure 3 first row, show that mouse and pig :math:`TAI_s` trends are highly significant (p-values: 1.36e-4 and 3.12e-4), whereas opossum is not (p-value 0.256). Interestingly, we found that bat's results are also statistically significat (p-value 1.11e-4). Based on this, we concluded that our hypothesis is supported only by the results on pig and mouse. The lack of statisticall significance of opossum results may stem from the reduced number of evolutionary ages that could be assigned to its genes (*ps* 1 and *ps* 4). As for bat, the widening use of young genes across across development may be the result of its highly specialzed fore-limbs which requires a growing use of bat-specific genes.

Regardless of the fact that, according to the best of our knowledge, there isn't cross species limb morphological divergence information at hindlimbs, we repeated the previous analysis on hindlimbs. The results, Figure 2 right, show that bat, and mouse transcriptomes are dominated at early and late stages by old genes, with young genes dominating the middle stage. Pig, has the opossite trend. All three species have statistically significant changes on :math:`TAI_s` values (p-values: bat 5.08e-4, pig 3.12e-4, and mouse 2.29e-10 ). As before, opossum didn't yield statistically significant results (p-value 0.151).


.. _Figure 3:

.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/bat_vtai.FL.svg
   :width: 23 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/pig_vtai.FL.svg
   :width: 23 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/mouse_vtai.FL.svg
   :width: 23 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/opossum_vtai.FL.svg
   :width: 23 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/bat_vtai.HL.svg
   :width: 23 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/pig_vtai.HL.svg
   :width: 23 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/mouse_vtai.HL.svg
   :width: 23 % 
.. image:: https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/opossum_vtai.HL.svg
   :width: 23 % 
Figure 3: Fore and Hind limbs distribution of :math:`VTAI` surrogates. The **first** and **second** rows correspond to the results of each species at fore and hind-limbs. When the value of :math:`VTAI` was not plot if its value was to large compared to the typical values of its surrogates. 

Bibliography
------------

.. bibliography:: Mendeley.bib
   :style: plain

