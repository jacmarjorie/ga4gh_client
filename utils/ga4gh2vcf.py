import sys, json, argparse
import vcf
import collections

# CallData is a class inherited by PyVCF
class CallData(collections.namedtuple('calldata', ['GT'])):
    __slots__ = ()

    _types = []
    _nums = []

    def __str__(self):
        dat = ", ".join(["%s=%s" % (x, y)
            for (x, y) in zip(self._fields, self)])
        return "CallData(" + dat + ')'

    def __reduce__(self):
        args = super(CallData, self).__reduce__()
        return make_calldata_tuple, (fields, )

class ga4gh2vcf():

	def __init__(self, garesponse):
		self.variants = garesponse['variants']
		self.callsets = self.getAllCallSets(self.variants)

	def getAllCallSets(self, variants):
		"""
		Get all calls
		"""
		allcs = []
		for v in variants:
			allcs += [c['callSetName'] for c in v['calls']]
		return sorted(set(allcs))

	def writeRecords(self, vcf_out):

		# prepare file
		self.writeVCFHeader(vcf_out)
		vcf_reader = vcf.Reader(filename=vcf_out)
		vcf_writer = vcf.Writer(open(vcf_out, 'w'), vcf_reader)

		for v in self.variants:
			# make record object, call object will refer back to this
			# omitting info list for now
			r = vcf.model._Record(v['referenceName'], v['start'], str(v['id']), v['referenceBases'], v['alternateBases'], 100, "PASS", v['info'], "GT", {})

			# set up samples objects
			variant_calls = {c['callSetName']:CallData("|".join(map(str,c['genotype']))) for c in v['calls']}
			sample_indexes = dict(zip(self.callsets, [CallData("0|0")]*len(self.callsets)))
			sample_indexes.update(variant_calls)

			# record.samples
			samples = []
			for c, g in sample_indexes.items():
				vcf_call = vcf.model._Call(r, c, g)
				samples.append(vcf_call)
				sample_indexes[c] = samples.index(vcf_call)

			# update sample portion
			r.samples = samples
			r._sample_indexes = sample_indexes

			# write vcf multisample object to file
			vcf_writer.write_record(r)

		vcf_writer.close()

	def writeVCFHeader(self, vcf_out, close=False):
		with open(vcf_out, 'w') as new_vcf:
			main_header = """##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=2016316
##reference=hg19
##source=fromGA4GH
##contig=<ID=1,assembly=b37,length=249250621>
##contig=<ID=2,assembly=b37,length=243199373>
##contig=<ID=3,assembly=b37,length=198022430>
##contig=<ID=4,assembly=b37,length=191154276>
##contig=<ID=5,assembly=b37,length=180915260>
##contig=<ID=6,assembly=b37,length=171115067>
##contig=<ID=7,assembly=b37,length=159138663>
##contig=<ID=8,assembly=b37,length=146364022>
##contig=<ID=9,assembly=b37,length=141213431>
##contig=<ID=10,assembly=b37,length=135534747>
##contig=<ID=11,assembly=b37,length=135006516>
##contig=<ID=12,assembly=b37,length=133851895>
##contig=<ID=13,assembly=b37,length=115169878>
##contig=<ID=14,assembly=b37,length=107349540>
##contig=<ID=15,assembly=b37,length=102531392>
##contig=<ID=16,assembly=b37,length=90354753>
##contig=<ID=17,assembly=b37,length=81195210>
##contig=<ID=18,assembly=b37,length=78077248>
##contig=<ID=19,assembly=b37,length=59128983>
##contig=<ID=20,assembly=b37,length=63025520>
##contig=<ID=21,assembly=b37,length=48129895>
##contig=<ID=22,assembly=b37,length=51304566>
##contig=<ID=X,assembly=b37,length=155270560>
##contig=<ID=Y,assembly=b37,length=59373566>
##contig=<ID=M,assembly=b37,length=16569>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
"""
			new_vcf.write(main_header)
			# samples are sorted numerically
			new_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(self.callsets)+"\n")
			if close:
				new_vcf.close()

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Util script to make a VCF from a GA4GH Request.')
	parser.add_argument('-i', '--input', dest='input', type=str, required=True, default=None, help="GAResponse JSON from file.")
	parser.add_argument('-o', '--output', dest='output', type=str, required=True, default=None, help="VCF File to be created.")

	args = parser.parse_args()

	with open(args.input) as garesponse:
		g = ga4gh2vcf(json.load(garesponse))
		g.writeRecords(args.output)
