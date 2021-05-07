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
	input_options.add_argument("-ref_name", dest='ref_name', type = str, \
		help='Name of the reference genome.', required=True)
	input_options.add_argument("-in_dir", dest='in_dir', type = str, \
		help='Name of the folder that contains the .bed files with coverage (=presence in other genomes)')', \
		required=True)

	settings2 = parser.add_argument_group('Classification settings')
	settings2.add_argument("-min_coverage_core", dest='min_coverage_core', type=int, \
		help = "Minimum coverage (number of genomes) for a region to be considered core.", required=True)

	output_options = parser.add_argument_group('Output settings')
	output_options.add_argument("-out_dir", dest='out_dir', type = str, help='Name of the output folder', required=True)

	args = parser.parse_args()

	# check inputs:
	assert os.path.exists(args.in_dir), args.in_dir+" does not exist"
	if args.in_dir[-1] != '/': args.in_dir += '/'

	bed_files = glob.glob(args.in_dir+args.ref_name+'.vs.*.coords')
	assert len(coords_files) > 0, "Can not find coords files in "+args.in_dir


	assert os.path.exists(args.out_dir), args.out_dir+" does not exist"
	if args.out_dir[-1] != '/': args.out_dir += '/'


	# setup folder structure and create logfile:
	list_out_dir = args.out_dir+'output/'
	os.system('mkdir -p '+args.out_dir)
	os.system('mkdir -p '+list_out_dir)

	logfilename = args.out_dir+'README-size'
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

	return args, coords_files, list_out_dir, logfile





def get_core_size(genome_coverage_bedfiles, min_number_of_genomes, outfile, logfile):
	'''
	Lists size of core regions in a given genome

	Args:
		genome_coverage_bedfiles (list) : list of bedfiles with 'coverage', number of alignments with other assemblies
		min_number_of_genomes (int) : minimum number of other assemblies a regions needs to be aligned \
		with for it to be classified as 'core'
		outfile (str) : the name of the file that will contain for strain the sizes of core regions
		logfile (File): logfile

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
