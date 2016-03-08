import sys, os, argparse, requests, json
from garequests import GASearchVariantsRequest
from threading import Thread, active_count

def all_chrs():
	# chrs = {'16':90354753}
	chrs = {'1':249250621,'2':243199373,'3':198022430,'4':191154276,'5':180915260,'6':171115067,'7':159138663,'8':146364022,'9':141213431,'10':135534747,'11':135006516,'12':133851895,'13':115169878,'14':107349540,'15':102531392,'16':90354753,'17':81195210,'18':78077248,'19':59128983,'20':63025520,'21':48129895,'22':51304566,'X':155270560,'Y':59373566}
	all_chrs_reqs = [GASearchVariantsRequest(start=1, end=v, referenceName=k) for k, v in chrs.items()]
	return all_chrs_reqs

def searchVariants(url, GARequest, pass_through, outfile):
	"""
		POST GASearchVariantsRequest to /variants/search
		return empty GASearchVariantsResponse if pass through option set
	"""
	payload = GARequest.getJson()
	try:
		r = requests.post(url + "/variants/search", data=payload, headers={'content-type':'application/json'})
		response = r.content
	except requests.exceptions.ConnectionError as e:
		# connection issue
		if pass_through:
			print "\nIssue sending request" + str(payload) + '\n\n'+ str(e)
			print "Proceeding with pass through."

		response = json.dumps({'variants':[], 'nextPageToken': None})
		
		# return "\nIssue sending request" + str(payload) + '\n\n'+ str(e)

	outf = open(outfile, "w")
	oncotator_build(json.loads(response), outf)

def oncotator_build(response, outfile):
		build = 'hg19'

		# write results to oncotator input file
		for variant in response['variants']:
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

def parallelRequests(url, all_requests, outdir, pass_through, combinedOutputFile):
	r_threads = list()
	count = 0

	for request in all_requests:
		count+=1
		print 'Starting request process', request
		outfile = outdir + "/" + str(count) + "_interm.csv"
		r_thread = Thread(target=searchVariants,
											name="Thread " + str(count) + " :" + outfile,
											args=(url, request, pass_through, outfile))
		r_thread.start()
		r_threads.append((r_thread, outfile))

	combinedOutput = outdir + "/" + combinedOutputFile
	countCombined = 0
	isCombined = [False] * len(r_threads)

	while( active_count() > 0 and countCombined != len(r_threads) ) :
		index = -1
		for (r_thread, outFile) in r_threads:
			index += 1
			if( isCombined[index] ):
				continue
			if( r_thread.is_alive() ):
				r_thread.join(120)
			else:
				os.system('cat %s >> %s'%(outFile, combinedOutput))
				countCombined += 1
				isCombined[index] = True


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='GA4GH Client')
	parser.add_argument('-u', '--urls', dest='urls', nargs = "+", type=str, required=True, default=None, help="Path to databases")
	parser.add_argument('-d', '--outdir', dest='outdir', type=str, required=True, default=None, help='Output Destination')
	parser.add_argument('-o', '--output', dest="output", type=str, required=True, default=None, help='Final output file')

	parser.add_argument('-a', '--all', dest="all", action='store_true', help='Flag to query across all chromosomes.')

	parser.add_argument('-vs', '--variantSetIds', nargs='+', type=str, required=False, default=[""],
		help='The IDs of the variant sets to search over.')
	parser.add_argument('-c', '--callSetIds', nargs='+', type=str, required=False, default=None,
		help='Only return variant calls which belong to call sets with these IDs. Leaving this blank returns all variant calls.')

	parser.add_argument('-r', '--referenceName', type=str, required=False, default=None,
		help='Only return variants on this reference.')
	parser.add_argument('-s', '--start', type=long, required=False, default=None,
		help='Required. The end of the window (0-based, exclusive) for which overlapping variants should be returned.')
	parser.add_argument('-e', '--end', type=long, required=False, default=None,
		help='The beginning of the window (0-based, inclusive) for which overlapping variants should be returned. Genomic positions are non-negative integers less than reference length. Requests spanning the join of circular genomes are represented as two requests one on each side of the join (position 0).')
	
	parser.add_argument('-p', '--pageSize', type=int, required=False, default=None,
		help='Specifies the maximum number of results to return in a single page. If unspecified, a system default will be used.')
	parser.add_argument('-t', '--pageToken', type=int, required=False, default=None,
		help='The continuation token, which is used to page through large result sets. To get the next page of results, set this parameter to the value of nextPageToken from the previous response.')
	
	parser.add_argument('-P', '--pass_through', action='store_true', help='If set, issues with requests will return empty GASearchVariantsResponse.')

	args = parser.parse_args()

	if not args.all and (not args.referenceName or not args.start or not args.end):
		parser.error('-r/--referenceName and -s/--start and -e/--end required if -a/--all not specified.')
	elif args.all:
		print 'querying whole genome'
		parallelRequests(args.urls[0], all_chrs(), args.outdir, args.pass_through, args.output)