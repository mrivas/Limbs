#for file in /data2/rivasas2/limbs/cufflinks_time_series/allSpecies/new_2014/*/*txt; do ln -s $file . ; done
#for file in /data2/rivasas2/limbs/trinity/bat/blast/*pl; do ln -s $file .; done
#for file in /data2/rivasas2/limbs/trinity/bat/blast/*genes.fasta; do ln -s $file .; done
#for file in /data2/rivasas2/limbs/trinity/bat/blast/*outfmt6; do ln -s $file .; done
#for file in /data2/rivasas2/limbs/trinity/bat/blast/bat*.orthologs.txt; do ln -s $file .; done
# for file in /data2/rivasas2/limbs/trinity/bat/blast/bat*.orthologs.uniq.txt; do ln -s $file .; done
#ln -s /data2/rivasas2/limbs/trinity/bat/blast/getAllOrthologs2.py .
#ln -s /data2/rivasas2/limbs/trinity/bat/blast/all.orthologs.uniq.txt .
#for file in /data2/rivasas2/limbs/cufflinks_time_series/allSpecies/new_2014/*/*svg; do ln -s $file; done
#for file in /data2/rivasas2/limbs/cufflinks_time_series/allSpecies/new_2014/*/analysis.*; do ln -s $file .; done
#
###########################################################################
## List of conserved and divergent ortholog genes
#
#declare -A levels
#levels["25"]=75
#levels["5"]=95
#echo "List of divergent and conserved genes" >  conserved_divergent.table.rst
#echo "=====================================" >> conserved_divergent.table.rst
#echo "   " >> conserved_divergent.table.rst
#echo ".. csv-table:: Conserved and divergent list of genes across species." >> conserved_divergent.table.rst
#echo "   :header: Conserved, Divergent" >> conserved_divergent.table.rst
#echo "   " >> conserved_divergent.table.rst
#for limb in FL HL; do
#for stage in beginning early late; do
#for level in 25 5; do
#	fileCon=$( ls * | grep $level | grep $stage | grep $limb | grep conserved | awk '{split($0,a,"%.txt");print a[1]"%25.txt"}' )
#	fileDiv=$( ls * | grep ${levels[$level]} | grep $stage | grep $limb | grep divergent | awk '{split($0,a,"%.txt");print a[1]"%25.txt"}' )
#	
#	nameCon=conserved_${limb}_${stage}${level}"%"
#	nameDiv=divergent_${limb}_${stage}${levels[$level]}"%"
#	
#	linkCon="\`$nameCon <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/${fileCon}>\`_"
#	linkDiv="\`$nameDiv <https://132.239.135.28/public/limbs/files/betweenSpeciesNewData/${fileDiv}>\`_"
#	
#	echo -e "   "${linkCon}","${linkDiv} >> conserved_divergent.table.rst
#
#done; done; done


##################################################################
# Phylography figures
#for file in /data2/rivasas2/limbs/phylostratum/*svg; do
#	ln -s $file .
#done

##################################################################
# Phylography TAI figures forward limbs
#for file in /data2/rivasas2/limbs/phylostratum/FL/*svg; do
#	ln -s $file .
#done
ln -s /data2/rivasas2/limbs/phylostratum/HL/*R .
# Phylography TAI figures hind limbs
#for file in /data2/rivasas2/limbs/phylostratum/HL/*svg; do
#	ln -s $file .
#done
ln -s /data2/rivasas2/limbs/phylostratum/FL/*R .
#ln -s /data2/rivasas2/limbs/docs/source/figures/betweenSpeciesNewData/tai.FL-HL.svg . 
#ln -s /data2/rivasas2/limbs/docs/source/figures/betweenSpeciesNewData/phylogenetic_tree.svg . 

#ln -s /data2/rivasas2/limbs/phylostratum/phylostratum.py .


