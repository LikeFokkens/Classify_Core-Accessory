'''
Given a set of genome alignments, classify contigs/chromosomes are "core" or "accessory".

Background:
In many genomes we can distinguish conserved regions, that are present in all genomes in 
a dataset and often referred to as "Core", and regions that are not conserved in all 
genomes -and thus probablynot absolutely required for survival- and often referred to 
as "Accessory" or "(Conditionally) dispensable" regions.

Here, given a dataset of aligned genomes, we classify a contig/scaffold/chromosome in 
an assembly as "Core" or "Accessory" based on how often it is aligned to regions in the 
other genomes. Note that whether a contig/scaffold/chromosome is designated as core depends 
on the other genomes in the dataset. Typically as a dataset increases in size, more cases 
are found in which (part of) a contig/scaffold/chromosome) is absent and thus as a dataset 
increases the number of core contigs/scaffolds/chromosomes decrease. Whether a 
contig/scaffold/chromosome is designated "Core" or "Accessory" also depends on the two main 
thresholds that can been chosen, namely how much of a contig has to be aligned 
(min_overlap_contig) in what fraction of the genomes in the dataset (min_overlap_dataset) 
for it to be classified as core. Finally you may chose to exclude alignments based on low 
levels of sequence similarity or short length.


Dependencies:
Bedtools (tested with bedtools v2.27.1).

Alignments are assumed to have been generated with nucmer and show_coords (both part of 
MUMmer (tested with version 3.23))

'''

import glob, argparse, os
import nucmercoords_to_bed

description = '''
Given a set of genome alignments, classify contigs/chromosomes are "core" or "accessory".
'''

def init():
	'''
	Parse commandline arguments, check input, setup directory structure
	'''
	
	parser = argparse.ArgumentParser(description=description)
	
	input_options = parser.add_argument_group('Input options')
	input_options.add_argument("-genome_file", dest='genome_file', type = str, \
		help='genomefile (contig {tab} length) of the reference genome.', required=True)
	input_options.add_argument("-ref_name", dest='ref_name', type = str, \
		help='Name of the reference genome.', required=True)
	input_options.add_argument("-in_dir", dest='indir', type = str, \
		help='Name of the folder that contains the .coords files generated by MUMmer.', required=True)
	settings1 = parser.add_argument_group('Alignment filters')
	settings1.add_argument("-min_length", dest='min_length', type=str, default = '1000', \
		help = "Minimum length of alignment to be included")
	settings1.add_argument("-min_identity", dest='min_identity', type=str, default = '90.0', \
		help = "Minimum percent identity of alignment to be included")	

	settings2 = parser.add_argument_group('Classification settings')
	settings2.add_argument("-min_overlap_contig", dest='min_overlap_contig', type=str, default = '0.75', \
		help = "Percent of contig that should be aligned to a query contig")
	settings2.add_argument("-min_overlap_dataset", dest='min_overlap_dataset', type=float, default = 0.90, \
		help = "Percent of queries that have >= min_overlap_contig alignments with a reference contig")
	
	output_options = parser.add_argument_group('Output settings')
	output_options.add_argument("-out_dir", dest='out_dir', type = str, help='Name of the output folder', required=True)
	
	args = parser.parse_args()

	# checks:
	assert os.path.exists(args.in_dir), args.in_dir+" does not exist"

	coords_files = glob.glob(args.in_dir+"/"+args.ref_name+'.vs.*')
	assert len(coords_files) > 0, "Can not find coords files in "+args.in_dir
	assert os.path.exists(args.out_dir), args.out_dir+" does not exist"
	assert os.path.exists(args.genome_file), args.genome_file+" does not exist"

	# setup folder structure and create logfile:
	bed_out_dir = args.out_dir+'/bedfiles/'
	list_out_dir = args.out_dir+'/output/'
	os.system('mkdir -p '+bed_out_dir)
	os.system('mkdir -p '+list_out_dir)

	logfilename = args.out_dir+'/README'
	logfile = open(logfilename, 'w')
	cmnd = ''
	for argv in sys.argv:
		cmnd += argv+' '
	logfile.write('Called with:\n'+cmnd+'\n\n')

	report = '\n\n##################################\n#\n#   SETTINGS\n'
	argsdict = vars(args)
	for var in argsdict.keys():
		report += '#\t'+var+'\t'+str(argsdict[var])+'\n'
	report += '#\n##################################\n\n'
	logfile.write(report)
	
	return args, coords_files, bed_out_dir, list_out_dir, logfile


def coords_to_beds(coords_files, min_identity, min_length, out_dir, logfile):
	'''
	For each coords file in the list, creates a bedfile with non-redundant aligned regions of a query from MUMmer output
	
	Args:
		coords_files (list):  list of tab-separated files generated from delta files (nucmer output) by show-coords.
		min_length (int) : minimum length of the aligned region for it to be included
		min_identify (float) : minimum percent identity of the aligned region for it to be included.
		out_dir (str): name of the output folder. Bedfiles created will be saved here as <coords_file>.bed and <coords_file>.nr.bed
		logfile (File): logfile to which successfully executed commands are written

	Returns:
	list with the names of the bedfiles
	'''

	bedfiles = []

	for coords_file in coords_files:
		bedfile = nucmercoords_to_bed.coords_to_bed(coords_file, min_identity, min_length, out_dir, logfile)
		if bedfile != None: 
			bedfiles.append(bedfile)

	return bedfiles



def bedfiles_to_genome_coverage(bedfiles, genome_file, out_dir, logfile):
	'''
	Calculates for each genomic region how many bedfiles overlap

	Args:
		bedfiles (list): list of bedfiles
		genome_file (str): tab-separated file with the length for each contig/scaffold/chromosome
		out_dir (str): name of the output folder
		logfile (File): logfile to which successfully executed commands are written

	Returns:
	Name of the file with the coverage per region. See bedtools genomecov documentation for more detail
	'''

	# output
	coverage_bed = out_dir+os.path.basename(genome_file)+'.coverage.bed'

	cov_cmnd = 'bedtools genomecov -bga -split -i'
	for bed in bedfiles:
		cov_cmnd += ' '+bed
	cov_cmnd +=' -g '+genome_file + ' > '+coverage_bed

	if os.system(cov_cmnd) == 0:
		logfile.write('\nCalculated coverage with:\n'+cov_cmnd+'\n')

	return coverage_bed


def classify_as_core_or_accessory(genome_coverage_bed, min_number_of_genomes, genome_file, min_overlap_contig, logfile):
	'''
	Lists core and accessory contigs for a given genome

	Args:
		genome_coverage_bed (str) : 
		min_number_of_genomes (int) : 
		genome_file (str) : 
		min_overlap_contig (float) : 
		logfile (File) : 

	Returns:
	The name of the tab-separated textfile that is produced: contig/scaffold/chromosome{tab}core/accessory/unique
	The cumulative size of core regions (more exact than classfying contigs/scaffolds/chromosomes)
	'''

	# filter out regions with enough coverage with awk
	

	cmnd_extract_core = "export LC_NUMERIC=C\nawk -v OFS='\\t' '{if ($4 >= "+min_identity+" && $7 >= "+min_length+") print $12,$1,$2}' "\
	+coords_file+" > "+awk_out_bed
	#calculate the size of these regions, also keep cumulative size
	#test whether this is more than min_overlap


	return None


if __name__ == "__main__":

	args, coords_files, bed_out_dir, list_out_dir, logfile = init()
	bedfiles = coords_to_beds(coords_files, args.min_identity, args.min_length, bed_out_dir, logfile)
	coverage_bed = bedfiles_to_genome_coverage(bedfiles, genome_file, list_out_dir, logfile)

	min_number_of_genomes = round(min_coverage*len(coords_files),0)
	#classify_as_core_or_accessory(coverage_bed, , args.genome_file, args.min_overlap_contig, args.min_coverage, logfile)








