'''
Given a set of genome alignments, classify contigs/chromosomes are "core" or "accessory".

Background:
In many genomes we can distinguish conserved regions, that are present in all genomes in 
a dataset and often referred to as "Core", and regions that are not conserved in all 
genomes -and thus probably not absolutely required for survival- and often referred to 
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
MUMmer (tested with version 3.23)). ALignment filenames are assumed to be formatted as 
follows: REF.vs.QRY.nucmer_maxmatch.coords. The extension .nucmer_maxmatch.coords can be 
changed with -coords_ext.

'''

import glob, argparse, os, sys
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
	input_options.add_argument("-in_dir", dest='in_dir', type = str, \
		help='Name of the folder that contains the .coords files generated by show-coords in MUMmer. '+ \
		'Names of coords files are assumed to be ref_name.vs.qry_name.nucmer_maxmatch.coords by default.'+ \
		'The file extension can be changed with -coords_ext', required=True)
	input_options.add_argument("-coords_ext", dest='coords_ext', type = str, default = '.nucmer_maxmatch.coords', \
		help='Extension of coords files, default is ".nucmer_maxmatch.coords".')
		
	settings1 = parser.add_argument_group('Alignment filters')
	settings1.add_argument("-min_length", dest='min_length', type=str, default = '1000', \
		help = "Minimum length of alignment to be included")
	settings1.add_argument("-min_identity", dest='min_identity', type=str, default = '90.0', \
		help = "Minimum percent identity of alignment to be included")	

	settings2 = parser.add_argument_group('Classification settings')
	settings2.add_argument("-include_self_alignments", dest='include_self', default = False, action='store_true', \
		help = "Set this flag of you also want to include alignments of the reference genome with itself.")
	settings2.add_argument("-min_overlap_contig", dest='min_overlap_contig', type=str, default = '0.75', \
		help = "Percent of contig that should be aligned to a query contig")
	settings2.add_argument("-min_overlap_dataset", dest='min_overlap_dataset', type=str, default = '0.90', \
		help = "Percent of queries that have >= min_overlap_contig alignments with a reference contig")
	
	output_options = parser.add_argument_group('Output settings')
	output_options.add_argument("-out_dir", dest='out_dir', type = str, help='Name of the output folder', required=True)
	
	args = parser.parse_args()

	# check inputs:
	assert os.path.exists(args.in_dir), args.in_dir+" does not exist"
	if args.in_dir[-1] != '/': args.in_dir += '/'
	
	assert os.path.exists(args.genome_file), args.genome_file+" does not exist"
		
	coords_files = glob.glob(args.in_dir+args.ref_name+'.vs.*.coords')
	assert len(coords_files) > 0, "Can not find coords files in "+args.in_dir
	
	#remove any self-alignments of the reference assembly
	if not args.include_self:
		self_coords = args.in_dir+args.ref_name+'.vs.'+args.ref_name+args.coords_ext
		#print(self_coords, len(coords_files))
		if self_coords in coords_files:
			coords_files.remove(self_coords)
			#print(len(coords_files))

	assert os.path.exists(args.out_dir), args.out_dir+" does not exist"
	if args.out_dir[-1] != '/': args.out_dir += '/'


	# setup folder structure and create logfile:
	args.out_dir = args.out_dir+'min_length-'+args.min_length+	'_min_identity-'+args.min_identity+ \
	'_min_overlap_contig-'+args.min_overlap_contig+'_min_overlap_dataset-'+args.min_overlap_dataset+'/'

	bed_out_dir  = args.out_dir+'bedfiles/'
	list_out_dir = args.out_dir+'output/'
	os.system('mkdir -p '+args.out_dir)
	os.system('mkdir -p '+bed_out_dir)
	os.system('mkdir -p '+list_out_dir)

	logfilename = args.out_dir+'README'
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
	cat_cmnd = 'cat'
	for coords_file in coords_files:
		bedfile = nucmercoords_to_bed.coords_to_bed(coords_file, min_identity, min_length, out_dir, logfile)
		if bedfile != None: 
			bedfiles.append(bedfile)
			cat_cmnd+=' '+bedfile
	
	assert len(bedfiles) == len(coords_files), 'failed to generate a bedfile for each .coords file'

	#concatenate all bedfiles and sort
	concat_bed_name = out_dir+'all_alignments.bed'
	cat_cmnd += ' | sort -k1,1 -k2,2n > '+concat_bed_name
	if os.system(cat_cmnd) == 0:
		logfile.write(cat_cmnd+'\n')
		return concat_bed_name

	else:
		print("Failed to concatenate and sort bedfiles")
		print(cat_cmnd) 
	
	return None

	


def bedfiles_to_genome_coverage(bedfile, genome_file, out_dir, logfile):
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

	cov_cmnd = 'bedtools genomecov -bga -split -i '+bedfile+' -g '+genome_file+' > '+coverage_bed

	if os.system(cov_cmnd) == 0:
		logfile.write('\nCalculated coverage with:\n'+cov_cmnd+'\n')
		return coverage_bed
	else:
		print('Failed to generate coverage bed with cmnd:')
		print(cov_cmnd)
		return None
		


def classify_as_core_or_accessory(genome_coverage_bed, min_number_of_genomes, genome_file, min_overlap_contig, outfile, logfile):
	'''
	Lists core and accessory contigs for a given genome

	Args:
		genome_coverage_bed (str) : bedfile with 'coverage', number of alignments with other assemblies
		min_number_of_genomes (int) : minimum number of other assemblies a regions needs to be aligned \
		with for it to be classified as 'core'
		genome_file (str) : file with the size per contig/chromosome/scaffold
		min_overlap_contig (float) : the fraction of contig that needs to be classified as core for the \
		whole contig to be classified as core
		outfile (str) : the name of the file that will contain for each contig whether it is core or accessory
		logfile (File) : logfile, this function adds the cumulative size of core and accessory regions in this assembly

	Returns:
	The name of the tab-separated textfile that is produced: contig/scaffold/chromosome{tab}core/accessory
	The cumulative size of core regions (more exact than classfying contigs/scaffolds/chromosomes)
	'''

	contig2length = {}
	with open(genome_file) as file:
		lines = file.readlines()
		for line in lines:
			contig, length = line.split()
			contig2length[contig] = int(length) 

	contig2cumsize_coreregions = {}
	with open(genome_coverage_bed) as file:
		lines = file.readlines()
		for line in lines:
			contig, start, end, coverage = line.split()
			if not contig in contig2cumsize_coreregions: contig2cumsize_coreregions[contig] = 0

			if int(coverage) >= min_number_of_genomes:
				contig2cumsize_coreregions[contig] += int(end) - int(start)
	

	with open(outfile, 'w') as file:
		file.write('# Contigs/chromosomes are classified as "core" when at least '+str(100*min_overlap_contig)\
			+' percent of the contig is present in at least '+str(min_number_of_genomes)+' genomes.\n')

		cumsize_all_core_regions = 0
		cumsize_contigs          = 0
		for contig in contig2length:
			core_size = contig2cumsize_coreregions[contig]
			if core_size >= min_overlap_contig * contig2length[contig]:
				file.write(contig+'\tcore\t'+str(contig2length[contig])+'\n')
			else:
				file.write(contig+'\taccessory\t'+str(contig2length[contig])+'\n')
			cumsize_all_core_regions += core_size
			cumsize_contigs += contig2length[contig]

	logfile.write('Cumulative size core regions: '+str(cumsize_all_core_regions)+ \
		' out of '+str(cumsize_contigs)+' : '+str(cumsize_all_core_regions/cumsize_contigs)+ \
		', accessory: '+str(cumsize_contigs-cumsize_all_core_regions))


	return cumsize_all_core_regions, cumsize_contigs



if __name__ == "__main__":

	# parse commandline arguments, initialize output folder structure, etc.
	args, coords_files, bed_out_dir, list_out_dir, logfile = init()

	# create for each query genome a bedfile with the aligned regions in the reference genome 
	concat_bedfile = coords_to_beds(coords_files, args.min_identity, args.min_length, bed_out_dir, logfile)

	# calculate the coverage of non-redundant alignments in other genomes
	coverage_bed = bedfiles_to_genome_coverage(concat_bedfile, args.genome_file, list_out_dir, logfile)
	if coverage_bed == None: sys.exit()

	# determine the minimum number of genomes a regions needs to be aligned with in order for it to be classified as core
	min_number_of_genomes = round(float(args.min_overlap_dataset)*len(coords_files),0)
	# name of the file in which we write for each contig whether it is core or accessory
	outfilename = args.out_dir+args.ref_name+'.core-accessory.tab'
	# classify regions and contigs as core or accessory
	cumsize_all_core_regions, cumsize_contigs = classify_as_core_or_accessory(coverage_bed, \
		min_number_of_genomes, args.genome_file, float(args.min_overlap_contig), outfilename, logfile)
	print(args.ref_name, cumsize_all_core_regions, cumsize_contigs-cumsize_all_core_regions)

	logfile.close()







