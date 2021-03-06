import json
import os
import hashlib
from pathlib import Path
import boto3
from flask import send_from_directory, Response

dataDir = '/var/www/data/'
binDir = '/var/www/bin/'
tmpDir = '/var/www/tmp/'

gaf = dataDir + 'gene_association.sgd'
gtfScript = binDir + 'GOTermFinder.pl'

S3_BUCKET = os.environ['S3_BUCKET']
s3_root_url = "https://" + S3_BUCKET + ".s3.amazonaws.com/"

def set_download_file(filename):

    if filename.endswith('.ps'): # or filename.endswith('.svg'):
        return send_from_directory(tmpDir, filename, as_attachment=True, mimetype='application/text', attachment_filename=(str(filename)))

    if filename.endswith('.png') or filename.endswith('.svg'):
        return send_from_directory(tmpDir, filename, as_attachment=True, mimetype='image/png+svg', attachment_filename=(str(filename)))
        
    f = open(tmpDir + filename, encoding="utf-8")
    content = f.read()
    f.close()
    return "<pre>" + content + "</pre>"

def upload_file_to_s3(file, filename):

    filename = 'gotermfinder/' + filename

    s3 = boto3.client('s3')
    file.seek(0)
    s3.upload_fileobj(file, S3_BUCKET, filename, ExtraArgs={'ACL': 'public-read'})

    return s3_root_url + filename

def get_download_url(tmpFile):

    downloadFile = tmpDir + tmpFile
    if not tmpFile.endswith('.png'):
        thisFile = Path(str(downloadFile))
        md5sum = None
        with thisFile.open(mode="rb") as fh:
            md5sum = hashlib.md5(fh.read()).hexdigest()
        newFileName = downloadFile
        file_suffix = tmpFile.split('.')[-1]
        if file_suffix not in ['txt', 'svg', 'png', 'html', 'ps']:
            file_suffix = 'txt'
        if md5sum:
            tmpFile = md5sum + "." + file_suffix
            newFileName = tmpDir + tmpFile
        os.rename(downloadFile, newFileName)
        downloadFile = newFileName
    file = open(downloadFile, "rb")
       
    s3_url = upload_file_to_s3(file, tmpFile)
        
    return s3_url

def get_param(request, name):

    p = request.args
    f = request.form
    return f.get(name) if f.get(name) else p.get(name)

def parse_gaf_file(gaf_file):
    
    f = open(gaf_file, encoding="utf-8")
    namemapping = {}
    aliasmapping = {}
    found = {}
    for line in f:
        pieces = line.strip().split('\t')
        if pieces[0].startswith('!'):
            continue
        if pieces[1] in found:
            continue
        found[pieces[1]] = 1
        names = pieces[10].split('|')
        i = 0
        for name in names:
            namemapping[name] = pieces[2]
            if i > 0:
                aliasmapping[name] = pieces[2]
            
            i = i + 1
        namemapping[pieces[1]] = pieces[2]
        aliasmapping[pieces[1]] = pieces[2]
    
    f.close()

    return (namemapping, aliasmapping)

def create_gene_list(genelist_file, genes, namemapping, aliasmapping):

    fw = open(genelist_file, 'w')
    
    genes = genes.split('|')
    found = {}
    for gene in genes:
        name = gene
        if gene in namemapping:
            name = namemapping[gene]
        if name in found:
            continue
        found[name] = 1
        if gene in aliasmapping:
            gene = aliasmapping[gene]
        fw.write(gene + "\n")
    fw.close()
    
def get_html_content(id):

    htmlFile = tmpDir + id + '.html'
    imageHtmlFile = tmpDir + id + '_ImageHtml.html'
    imageFile = tmpDir + id + '_Image.html'

    ## get rid of html and body tags in html file
    f = open(htmlFile, encoding="utf-8")
    html = f.read()
    f.close()
    
    html = html.replace("<html><body>", "").replace("</body></html>", "")
    html = html.replace("color=red", "color=maroon")
    html = html.replace('<a name="table" />', '')
    html = html.replace("infowin", "_extwin")

    ## fix image URL in image file and get rid of html and body tags
    f = open(imageFile, encoding="utf-8")
    image = f.read()
    f.close()

    image = image.replace("<img src='./", "<img src='" + s3_root_url + 'gotermfinder/')

    fw = open(imageFile, "w")
    fw.write(image)
    fw.close()
    
    image = image.replace("<html><body>", "").replace("</body></html>", "")
    image = "<b>Nodes" + image.split("</font><br><br><b>Nodes")[1]

    ## fix image URL in imageHtml file
    f = open(imageHtmlFile, encoding="utf-8")
    imageHtml = f.read()
    f.close()

    imageHtml = imageHtml.replace("<img src='./", "<img src='" + s3_root_url)

    fw = open(imageHtmlFile, "w")
    fw.write(imageHtml)
    fw.close()
    
    return (html, image)

def enrichment_search(request, id):

    geneList = tmpDir + id + '.txt'
    tmpTab = tmpDir + id + '_tab.txt'
    
    genes = get_param(request, 'genes')
    genes = genes.upper().replace('SGD:', '')
    aspect = get_param(request, 'aspect')
    if aspect is None:
        aspect = 'P'
    (namemapping, aliasmapping) = parse_gaf_file(gaf)
    create_gene_list(geneList, genes, namemapping, aliasmapping)

    cmd = gtfScript + ' -a ' + aspect + ' -g ' + gaf + ' -t ' + tmpDir + ' ' + geneList + ' -F'
    output = os.popen(cmd).read()
    
    f = open(tmpTab)
    data = []
    for line in f:
        if line.startswith('GOID'):
            continue
        pieces = line.strip().split('\t')
        pvalue = ''
        if 'e-' in pieces[2]:
            ## if it is in scientific notation, we want up to two of the decimal places 
            # 6.036693772627e-26
            pvalue_raw = pieces[2].split('.')
            pvalue = pvalue_raw[0] + '.' + pvalue_raw[1][0:2] + 'e-' + pvalue_raw[1].split('e-')[1]
        elif '.' in pieces[2]:
            # otherwise, we'll take up to five decimal places
            pvalue_raw = pieces[2].split('.')
            pvalue = pvalue_raw[0] + '.' + pvalue_raw[1][0:5]
        else:
            pvalue = pieces[2]
        data.append({ 'goid': pieces[0],
                      'term': pieces[1],
                      'pvalue': pvalue,
                      'num_gene_annotated': pieces[4] })
    f.close()

    return data

def gtf_search(request, id):

    geneList = tmpDir + id + '.txt'
    gene4bgList = tmpDir + id + '_4bg.txt'
    tmpTab = tmpDir + id + '_tab.txt'
    
    # 'COR5|CYT1|Q0105|QCR2|S000001929|S000002937|S000003809|YEL024W|YEL039C|YGR183C|YHR001W-A'

    genes = get_param(request, 'genes')
    if genes is None:
        return { " ERROR": "NO GENE NAME PASSED IN" }
    genes = genes.upper().replace('SGD:', '')
    
    aspect = get_param(request, 'aspect')
    if aspect is None:
        aspect = 'F'        
    (namemapping, aliasmapping) = parse_gaf_file(gaf)
    create_gene_list(geneList, genes, namemapping, aliasmapping)

    genes4bg = get_param(request, 'genes4bg')
    evidence = get_param(request, 'evidence')
    pvalue = get_param(request, 'pvalue')
    if pvalue is None:
        pvalue = 0.01
        
    option = ""
    if pvalue == "" and pvalue != 0.01:
        option = " -p " + pvalue

    if genes4bg:
        genes4bg = genes4bg.replace('|', '\n')
        fw = open(gene4bgList, 'w')
        fw.write(genes4bg + "\n")
        fw.close()
        option = option + " -b " + gene4bgList
    
    if evidence:
        if evidence.endswith('|'):
            evidence = evidence[0:-1]
        if evidence.startswith('|'):
            evidence = evidence[1:]
        evidence = evidence.replace('|', ',')
        option = option + " -e " + evidence

    p = request.args
    f = request.form
    FDR = p.get('FDR') if p.get('FDR') else f.get('FDR')

    if FDR:
        option = option + " -F"
        
    cmd = gtfScript + ' -a ' + aspect + ' -g ' + gaf + ' -t ' + tmpDir + ' ' + geneList + ' -v'
    if option != '':
    	cmd = cmd + option

    output = os.popen(cmd).read()
    
    if 'No significant GO terms' in output:
        if 'were found for your input list of genes.' in output:
            output = output.split('were found for your input list of genes')[0] + ' were found for your input list of genes.'
        else:
            output = "No significant GO terms were found for your input list of genes."
        return { "output": output }
    # elif os.path.exists(tmpTab):
    #    return  { "output": "<pre>" + output + "</pre>" }
    else:
        png_url = get_download_url(id+'.png')
        (html, imageHtml) = get_html_content(id)        
        return { "html": html,
                 "image_html": imageHtml,
		 "image_page": get_download_url(id+'_Image.html'),
		 "tab_page": get_download_url(id+'_tab.txt'),
		 "term_page": get_download_url(id+'_terms.txt'),
                 "table_page": get_download_url(id+'.html'),
		 "png_page": png_url,
		 "svg_page": get_download_url(id+'.svg'),
		 "ps_page": get_download_url(id+'.ps'),
		 "input_page": get_download_url(id+'.txt') }
    
    




