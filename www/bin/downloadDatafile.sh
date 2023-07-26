#!/bin/sh

cd /var/www/data/new/

/usr/bin/wget http://snapshot.geneontology.org/annotations/sgd.gaf.gz
/usr/bin/wget http://snapshot.geneontology.org/ontology/go-basic.obo

/bin/gunzip -f sgd.gaf.gz

/bin/cp ../gene_association.sgd ../gene_association.sgd_old
/bin/cp ../gene_ontology.obo ../gene_ontology.obo_old

/usr/bin/grep -v "$(printf '\t')IBA$(printf '\t')" sgd.gaf | /usr/bin/grep -v "$(printf '\t')IEA$(printf '\t')" | /usr/bin/grep -v "$(printf '\t')CPX-" > ../gene_association_mapper.sgd 
/usr/bin/grep -v -v "$(printf '\t')CPX-" sgd.gaf > ../gene_association.sgd

/bin/mv go-basic.obo ../gene_ontology.obo
/bin/cp ../gene_association.sgd ../unparsed/
/bin/cp ../gene_ontology.obo ../unparsed/

echo "creating slim component gaf file..."

/var/www/bin/map2slim /var/www/data/slim_component.lst /var/www/data/gene_ontology.obo /var/www/data/gene_association_mapper.sgd -o /var/www/data/new/slim_component_gene_association.sgd

echo "creating slim process gaf file..."

/var/www/bin/map2slim /var/www/data/slim_process.lst /var/www/data/gene_ontology.obo /var/www/data/gene_association_mapper.sgd -o /var/www/data/new/slim_process_gene_association.sgd

echo "creating slim function gaf file..."

/var/www/bin/map2slim /var/www/data/slim_function.lst /var/www/data/gene_ontology.obo /var/www/data/gene_association_mapper.sgd -o /var/www/data/new/slim_function_gene_association.sgd

/bin/cp ../slim_component_gene_association.sgd ../slim_component_gene_association.sgd_old
/bin/cp ../slim_process_gene_association.sgd ../slim_process_gene_association.sgd_old
/bin/cp ../slim_function_gene_association.sgd ../slim_function_gene_association.sgd_old

/bin/mv slim_component_gene_association.sgd ../
/bin/mv slim_process_gene_association.sgd ../
/bin/mv slim_function_gene_association.sgd ../

