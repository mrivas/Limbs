for specie in mouse pig opossum; do
	for file in /var/www/public/limbs/files/${specie}/*table.rst; do 
		ln -s $file .
	done
done

ln -s /var/www/public/limbs/files/betweenSpeciesNewData/conserved_divergent.table.rst
