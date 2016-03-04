import sys, argparse, requests, json
from garequests import GASearchVariantsRequest

def searchVariants(url, GARequest, pass_through=False):
	"""
		POST GASearchVariantsRequest to /variants/search
		return empty GASearchVariantsResponse if pass through option set
	"""
	payload = GARequest.getJson()
	try:
		r = requests.post(url + "/variants/search", data=payload, headers={'content-type':'application/json'})
	except requests.exceptions.ConnectionError as e:
		# connection issue
		if pass_through:
			print "\nIssue sending request" + str(payload) + '\n\n'+ str(e)
			print "Proceeding with pass through."

			return json.dumps({'variants':[], 'nextPageToken': None})
		
		return "\nIssue sending request" + str(payload) + '\n\n'+ str(e)

	return r.content

def parallelRequest(urls, request_pool, pass_through=False, outdir, output):
	"""
		Consolidate search variant requests, and submit in parallel
	"""
	m_threads = list()
	count = 0

	for request in request_pool:
		count +=1

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='GA4GH Client')
	parser.add_argument('-u', '--urls', dest='urls', nargs = "+", type=str, required=True, default=None, help="Path to databases")
	parser.add_argument('-d', '--out_dir', dest='out_dir', type=str, required=True, default=None, help='Output Destination')
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
	else:
		# minimum set of information for now
		request = GASearchVariantsRequest(start=args.start, end=args.end, referenceName=args.referenceName, callSetIds=args.callSetIds)
		response = searchVariants(args.urls[0], request, pass_through=args.pass_through)
		print response