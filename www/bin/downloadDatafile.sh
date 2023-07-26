#!/bin/sh -x

OUTPUT_FILE=/tmp/output.log            # output written here
OUTPUT2_FILE=/tmp/output2.log          # processed output to be included in JSON email message
MESSAGE_JSON_FILE=/tmp/message.json    # pre-processed JSON email message
MESSAGE2_JSON_FILE=/tmp/message2.json  # post-processed JSON email message

echo "downloadDatafile.sh start:  $(/bin/date)" | /bin/tee $OUTPUT_FILE

if [ ! -d /var/www/data/new ]; then
	/bin/mkdir /var/www/data/new
	if [ $? -ne 0 ]; then
		echo "Cannot mkdir /var/www/data/new" | /bin/tee -a $OUTPUT_FILE
		exit $?
	fi
fi

cd /var/www/data/new/
/usr/bin/wget http://snapshot.geneontology.org/annotations/sgd.gaf.gz 2>&1 | /bin/tee -a $OUTPUT_FILE
if [ $? -ne 0 ]; then
	echo "Error: wget http://snapshot.geneontology.org/annotations/sgd.gaf.gz" | /bin/tee -a $OUTPUT_FILE
	exit $?
fi

/usr/bin/wget http://snapshot.geneontology.org/ontology/go-basic.obo 2>&1 | /bin/tee -a $OUTPUT_FILE
if [ $? -ne 0 ]; then
        echo "Error: wget http://snapshot.geneontology.org/ontology/go-basic.obo" | /bin/tee -a $OUTPUT_FILE
        exit $?
fi

/bin/gunzip -f sgd.gaf.gz 2>&1 | /bin/tee -a $OUTPUT_FILE

/bin/cp -p ../gene_association.sgd ../gene_association.sgd_old
/bin/cp -p ../gene_ontology.obo ../gene_ontology.obo_old

/usr/bin/grep -v "$(printf '\t')IBA$(printf '\t')" sgd.gaf | /usr/bin/grep -v "$(printf '\t')IEA$(printf '\t')" | /usr/bin/grep -v "$(printf '\t')CPX-" > ../gene_association_mapper.sgd 
/usr/bin/grep -v -v "$(printf '\t')CPX-" sgd.gaf > ../gene_association.sgd

/bin/mv go-basic.obo ../gene_ontology.obo
/bin/cp -p ../gene_association.sgd ../unparsed/
/bin/cp -p ../gene_ontology.obo ../unparsed/

echo "creating slim component gaf file..." | /bin/tee -a $OUTPUT_FILE

/var/www/bin/map2slim /var/www/data/slim_component.lst /var/www/data/gene_ontology.obo /var/www/data/gene_association_mapper.sgd -o /var/www/data/new/slim_component_gene_association.sgd 2>&1 | /bin/tee -a $OUTPUT_FILE

echo "creating slim process gaf file..." | /bin/tee -a $OUTPUT_FILE

/var/www/bin/map2slim /var/www/data/slim_process.lst /var/www/data/gene_ontology.obo /var/www/data/gene_association_mapper.sgd -o /var/www/data/new/slim_process_gene_association.sgd 2>&1 | /bin/tee -a $OUTPUT_FILE

echo "creating slim function gaf file..." | /bin/tee -a $OUTPUT_FILE

/var/www/bin/map2slim /var/www/data/slim_function.lst /var/www/data/gene_ontology.obo /var/www/data/gene_association_mapper.sgd -o /var/www/data/new/slim_function_gene_association.sgd 2>&1 | /bin/tee -a $OUTPUT_FILE

/bin/cp -p ../slim_component_gene_association.sgd ../slim_component_gene_association.sgd_old
/bin/cp -p ../slim_process_gene_association.sgd ../slim_process_gene_association.sgd_old
/bin/cp -p ../slim_function_gene_association.sgd ../slim_function_gene_association.sgd_old

/bin/mv slim_component_gene_association.sgd ../
/bin/mv slim_process_gene_association.sgd ../
/bin/mv slim_function_gene_association.sgd ../

echo "downloadDatafile.sh finished:  $(/bin/date)" | /bin/tee -a $OUTPUT_FILE

* create and send email report

# add \n characters to end of each line in OUTPUT_FILE for JSON message
/usr/bin/touch $OUTPUT2_FILE
/usr/bin/awk '{printf "%s\\n", $0}' $OUTPUT_FILE > $OUTPUT2_FILE

# create JSON email message
echo '{"Data": "From: '$(echo $EMAIL_FROM)'\nTo: '$(echo $EMAIL_TO)'\nSubject: downloadDatafile.sh report\nMIME-Version: 1.0\nContent-type: Multipart/Mixed; boundary=\"NextPart\"\n\n--NextPart\nContent-Type: text/plain\n\ndownloadDatafile.sh completed successfully\n\n--NextPart\nContent-Type: text/plain;\nContent-Disposition: attachment; filename=\"downloadDatafile_report.txt\"\n\n'$(cat $OUTPUT2_FILE)'\n--NextPart--"}' > $MESSAGE_JSON_FILE

# replace literal newline characters with literal '\n' characters in JSON message
/usr/bin/sed -i 's/$/\\n/' $MESSAGE_JSON_FILE
/usr/bin/touch $MESSAGE2_JSON_FILE
/usr/bin/tr -d '\n' < $MESSAGE_JSON_FILE > $MESSAGE2_JSON_FILE
/usr/bin/sed -i 's/}\\n/}\n/' $MESSAGE2_JSON_FILE  # add final trailing newline

/usr/local/bin/aws ses send-raw-email --cli-binary-format raw-in-base64-out --raw-message file://${MESSAGE2_JSON_FILE} --region $AWS_SES_REGION

exit 0
