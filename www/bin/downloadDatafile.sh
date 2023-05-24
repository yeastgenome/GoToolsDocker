#!/bin/sh -x

echo "downloadDatafile.sh start: $(/bin/date)

if [ ! -d /var/www/data/new ]; then
	/bin/mkdir /var/www/data/new
	if [ $? -ne 0 ]; then
		echo "Cannot mkdir /var/www/data/new"
		exit $?
	fi
fi

cd /var/www/data/new/

/usr/bin/wget http://snapshot.geneontology.org/annotations/sgd.gaf.gz
if [ $? -ne 0 ]; then
	echo "Error: wget http://snapshot.geneontology.org/annotations/sgd.gaf.gz"
	exit $?
fi

/usr/bin/wget http://snapshot.geneontology.org/ontology/go-basic.obo
if [ $? -ne 0 ]; then
        echo "Error: wget http://snapshot.geneontology.org/ontology/go-basic.obo"
        exit $?
fi

/bin/gunzip -f sgd.gaf.gz

/bin/cp -p ../gene_association.sgd ../gene_association.sgd_old
/bin/cp -p ../gene_ontology.obo ../gene_ontology.obo_old
/bin/mv sgd.gaf ../gene_association.sgd
/bin/mv go-basic.obo ../gene_ontology.obo

echo "creating slim component gaf file..."

/var/www/bin/map2slim /var/www/data/slim_component.lst /var/www/data/gene_ontology.obo /var/www/data/gene_association.sgd -o /var/www/data/new/slim_component_gene_association.sgd

echo "creating slim process gaf file..."

/var/www/bin/map2slim /var/www/data/slim_process.lst /var/www/data/gene_ontology.obo /var/www/data/gene_association.sgd -o /var/www/data/new/slim_process_gene_association.sgd

echo "creating slim function gaf file..."

/var/www/bin/map2slim /var/www/data/slim_function.lst /var/www/data/gene_ontology.obo /var/www/data/gene_association.sgd -o /var/www/data/new/slim_function_gene_association.sgd

/bin/cp -p ../slim_component_gene_association.sgd ../slim_component_gene_association.sgd_old
/bin/cp -p ../slim_process_gene_association.sgd ../slim_process_gene_association.sgd_old
/bin/cp -p ../slim_function_gene_association.sgd ../slim_function_gene_association.sgd_old

/bin/mv slim_component_gene_association.sgd ../
/bin/mv slim_process_gene_association.sgd ../
/bin/mv slim_function_gene_association.sgd ../

echo "downloadDatafile.sh end: $(/bin/date)

exit 0
