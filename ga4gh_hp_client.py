import sys
import argparse
import requests
import json

global projects, chrs, variantSetIds, url, outfile
projects = {'cervix': ['CESC-US'], 'stomach': ['GACA-CN', 'STAD-US'], 'bone': ['BOCA-FR', 'BOCA-UK'], 'breast': ['BRCA-UK', 'BRCA-US'], 'bladder': ['BLCA-CN', 'BLCA-US'], 'brain': ['PBCA-DE', 'LGG-US', 'NBL-US', 'GBM-US'], 'head and neck': ['ORCA-IN', 'THCA-SA', 'THCA-US'], 'lung': ['LUSC-CN', 'LUSC-KR', 'LUSC-US'], 'blood': ['MALY-DE', 'CLLE-ES', 'CMDI-UK', 'LAML-KR', 'ALL-US'], 'colorectal': ['COCA-CN', 'COAD-US', 'READ-US'], 'skin': ['SKCM-US'], 'uterus': ['UCEC-US'], 'esophagus': ['ESCA-CN', 'ESAD-UK'], 'pancreas': ['PACA-CA', 'PACA-IT', 'PACA-AU', 'PAEN-AU'], 'prostate': ['EOPC-DE', 'PRAD-UK', 'PRAD-CA', 'PRAD-US'], 'ovary': ['OV-AU', 'OV-US'], 'kidney': ['RECA-EU', 'RECA-CN', 'KIRC-US', 'KIRP-US'], 'liver': ['LIAD-FR', 'LICA-FR', 'LIHM-FR', 'LINC-JP', 'LIRI-JP', 'LIHC-US'], 'all': ['MALY-DE', 'BOCA-FR', 'ESCA-CN', 'RECA-EU', 'LIAD-FR', 'LICA-FR', 'LIHM-FR', 'EOPC-DE', 'BLCA-CN', 'CLLE-ES', 'PBCA-DE', 'COCA-CN', 'ORCA-IN', 'LUSC-CN', 'PACA-CA', 'PACA-IT', 'GACA-CN', 'CMDI-UK', 'BOCA-UK', 'ESAD-UK', 'RECA-CN', 'LINC-JP', 'LIRI-JP', 'PRAD-UK', 'LAML-KR', 'BRCA-UK', 'THCA-SA', 'LUSC-KR', 'OV-AU', 'PACA-AU', 'PRAD-CA', 'LGG-US', 'NBL-US', 'COAD-US', 'KIRC-US', 'BLCA-US', 'ALL-US', 'GBM-US', 'BRCA-US', 'CESC-US', 'READ-US', 'KIRP-US', 'LIHC-US', 'LUSC-US', 'OV-US', 'PRAD-US', 'SKCM-US', 'STAD-US', 'UCEC-US', 'THCA-US','PAEN-AU']}
chrs = {'1':249250621,'2':243199373,'3':198022430,'4':191154276,'5':180915260,'6':171115067,'7':159138663,'8':146364022,'9':141213431,'10':135534747,'11':135006516,'12':133851895,'13':115169878,'14':107349540,'15':102531392,'16':90354753,'17':81195210,'18':78077248,'19':59128983,'20':63025520,'21':48129895,'22':51304566,'X':155270560,'Y':59373566, 'M':16571}

# for testing
#chrs = {'1':249250621,'2':243199373,'3':198022430}
variantSetIds = {
    "cervix": 6, 
    "stomach": 16, 
    "liver": 10, 
    "skin": 15, 
    "bladder": 1, 
    "brain": 2, 
    "head and neck": 5, 
    "lung": 12, 
    "breast": 4, 
    "colorectal": 7, 
    "kidney": 9, 
    "uterus": 17, 
    "esophagus": 8, 
    "pancreas": 13, 
    "prostate": 14, 
    "ovary": 11, 
    "bone": 3, 
    "blood": 0
  }
outfile = open('oncotator_input.txt', 'w')
header = 'tumor_name\tchr\tstart_position\tend_position\tbuild\tref_allele\talt_allele\ttumor_f\tt_ref_count\tt_alt_count\n'

def GASearchVariantsRequest(query, url):
    request = {'start':1, 'end': None, 'referenceName': None, 'pageSize': None, 'pageToken': None, 'callSetIds': None, 'variantName': None, 'variantSetIds': variantSetIds[query.lower()]}
    job_list = []
    for key, val in chrs.iteritems():
        request['end'] = val
        request['referenceName'] = key
        oncotator_build(request, url)

def GASearchVariantsRequest1(query, url):
    request = {'start':1, 'end': None, 'referenceName': None, 'pageSize': None, 'pageToken': None, 'callSetIds': None, 'variantName': None, 'variantSetIds': variantSetIds[query.lower()]}
    job_list = []
    for key, val in chrs.iteritems():
        request['end'] = val
        request['referenceName'] = key
        oncotator_build(request, url)

def variants_search(this_request, this_url):
    r = requests.post(this_url + "/variants/search", data=json.dumps(this_request), headers={'content-type':'application/json'})
    return json.loads(r.text)

def vs_recurse(this_request, this_url):
    r_json = variants_search(this_request, this_url)
    responses = r_json['variants']
    nextPageToken = r_json['nextPageToken']

    while nextPageToken != None:
        new_request = this_request
        new_request['pageToken'] = nextPageToken
        new_response = variants_search(new_request, this_url)
        responses = responses + new_response['variants']
        nextPageToken = new_response['nextPageToken']
    return responses

def oncotator_build(request, query_url):
    build = 'hg19'

    for url in query_url:
        # request to url
        variants = vs_recurse(request, url)

        # write results to oncotator input file
        for variant in variants:
            ref_allele = variant['referenceBases']
            alt = variant['alternateBases']
            contig = variant['referenceName']
            start_position = variant['start']
            end_position = variant['end']
            # TODO add better checking
            for call in variant['calls']:
                for allele in call['genotype']:
                    if allele != '0':
                        line = []
                        if 'callSetId' != None:
                            line.append(call['callSetId'])
                        else:
                            line.append("None")
                        line.append(contig)
                        line.append(start_position)
                        line.append(end_position)
                        line.append('37')
                        line.append(ref_allele)
                        line.append(alt[int(allele)-1])
                        if ('AF' in call['info']):
                            line.append(call['info']['AF'][0])
                        else:
                            line.append(" ")
                        if ('AN' in call['info']):
                            line.append(call['info']['AN'][0])
                        else:
                            line.append(" ")
                        if ('AC' in call['info']):
                            line.append(call['info']['AC'][0])
                        else:
                            line.append(" ")
                        outfile.write('\t'.join(map(str, line)))
                        outfile.write('\n')

if __name__ == "__main__":

    # example query:
    # python ga4gh_hp_client.py -u http://<host1>:<port_on_host1>/variants/search http://<host2>:<port_on_host2> http://<host3>:<port_on_host3>/variants/search -q liver
    parser = argparse.ArgumentParser(description='Variant DB GA4GH Communicator')
    parser.add_argument('-u', "--urls", dest="urls", nargs = "+", type=str, required=True, default=None, help="Path to databases")
    parser.add_argument('-q', '--query', dest="query", nargs='?', default=None,
        help='Query by cancer type.')

    args = parser.parse_args()

    url = args.urls
    outfile.write(header)

    GASearchVariantsRequest(args.query, args.urls)
    outfile.close()
