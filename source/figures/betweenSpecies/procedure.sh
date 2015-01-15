for file in /data2/rivasas2/limbs/cufflinks_time_series/allSpecies/noUniq/*txt; do ln -s $file . ; done
for file in /data2/rivasas2/limbs/trinity/bat/blast/*pl; do ln -s $file .; done
for file in /data2/rivasas2/limbs/trinity/bat/blast/*genes.fasta; do ln -s $file .; done
for file in /data2/rivasas2/limbs/trinity/bat/blast/*outfmt6; do ln -s $file .; done
for file in /data2/rivasas2/limbs/trinity/bat/blast/bat*.orthologs.txt; do ln -s $file .; done
 for file in /data2/rivasas2/limbs/trinity/bat/blast/bat*.orthologs.uniq.txt; do ln -s $file .; done
ln -s /data2/rivasas2/limbs/trinity/bat/blast/getAllOrthologs2.py .
ln -s /data2/rivasas2/limbs/trinity/bat/blast/all.orthologs.uniq.txt .
for file in /data2/rivasas2/limbs/cufflinks_time_series/allSpecies/noUniq/*svg; do ln -s $file; done
for file in /data2/rivasas2/limbs/cufflinks_time_series/allSpecies/analysis.*; do ln -s $file .; done
