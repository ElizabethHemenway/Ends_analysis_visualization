#!/bin/bash

# ------------------------------------------------------------------------------------
# v2.1 by Colette L. Picard (updated 01/26/2015)

# This script aims to summarize CpG, CHG and CHH methylation within various genomic regions
# (specified by .bed files)
# The methylation data are in BED file format:
# chr	start	end		#unmethylated	#methylated		%methylated
# Chr1	12345	12346	9				1				10.0
# ...

# Inputs and flags: 
# 	-i infile.bed : summarized methylation data in BED format 
# 		= (chr, start, end, # unmethylated, # methylated, fraction methylated)

# 	-r regions.bed : regions of interest to be queried, in BED format
# 		= (chr, start, end, feature name, ., strand) - last 3 fields optional

#	-o outfile.bed : output will be saved to this file (WILL overwrite existing)

#	-m method : (optional) method for summarizing methylation within each region. Options are:
#		1 = (default) weighted mean
#		2 = fraction of sites with at least `cutoff' fraction of reads methylated

#	-c cutoff : (optional) cutoff to use when using method 2. Default = 0.5.

#	-x minReads : sites with fewer than minReads total reads are censored from the analysis.
		
		
# Pipeline:
# (P) Preprocessing: convert all 'chloroplast'>'ChrC', 'mitochondria'>'ChrM', 'chr'>'Chr'
# for consistency. Note not all input files require this.

# (A) Use intersectBed to match each site in methylation.bed with the corresponding region
# 	Command:
# 	bedtools intersect -a methylation.bed -b regions.bed > intersect.bed
	
# (B) Process intersect.bed to obtain a summary of the methylation status within each region,
# 	using one of two methods (see inputs above):
# 	Command:
# 	python summarizeMethylation.py intersect.bed outfile [-method] [-cutoff] [-minReads]
	
	
# Output:
# 	a list of the regions provided in regions.bed, with their summarized methylation levels,
# 	in BED format (chr, start, end, feature name, methylation level, numSites).


# Usage:
# 	sumByFeature.sh -i infile.bed -r regions.bed -o outfile.bed [-m method -c cutoff -x minReads]

# Fields in brackets are optional, defaults are: method = 1 (weighted mean), cutoff = 0.5,
# minReads = 5.

# ------------------------------------------------------------------------------------
# Version history

# v1.0 initial version											June 03, 2013 by Colette Picard
# v2.0:															June 10, 2013 by Colette Picard
#		1) added getopts support to make more user-friendly
#		2) to save space, no longer keep the BED file output - it is deleted after use
#			(therefore, also no longer check whether or not the BED file already exists - too
#			much risk that wrong/bad file could be used and assumed "correct" because of name.
#		3) no longer require user to specify an abbreviation for experiment, region or methylation
#			to facilitate automatic output file naming - just give your own output filenames!
#			** note that running this script with same output filename will overwrite previous **
#		4) working directory and output file destination no longer hardcoded
#		5) all files except for output file now kept as tempfiles - thus should be able to run
#		6) several occurrences of this script at once. Cleanup tested & working properly.
#		** note1: confirmed new version produces same output as old version **
# v2.1:															Dec. 26, 2015 by Colette Picard
#		1) removed requirement for temp/ folder - temporary files are stored in current dir
#		and deleted anyway
#		2) added support for chrc and chrm as alternatives to ChrC, ChrM, etc.
#		3) summarizeMethylation.py no longer assumed to be in sumByFeature subdir, can
#		give path manually with -s or default is same dir as this script
# v2.2:															Jun. 16, 2015 by Colette Picard
#		1) added usage and description info that is printed when the script is called with no args
# v3.0	1) added new flag -O to omit the methylation value within a designated interval for each
#		interval in -r. This assumes that the 4th position in the -r BED file will be e.g. 'Chr1:3502463-3502465'
#		which indicates to ignore any methylation data within that interval for that specific record
#		in the BED file.
# v.4.0 - added option to output all intervals with no overlapping reads of sufficient coverage in the output
#		file (by default these are omitted)
# v.5.0 - utilize Python3 instead of Python2	Mar. 14, 2023 by Elizabeth Hemenway
# ------------------------------------------------------------------------------------

# Description printed when "help" option specified:
read -d '' usage <<"EOF"
sumByFeature v.2.2	06/16/2015
Usage:
sumByFeature.sh [options] -i infile.bed -r regions.bed -o outfile.bed

Summarizes methylation data provided in infile.bed across the regions provided in regions.bed.
Assumes infile.bed contains one entry per single cytosine, and both input files are in BED format.
Outputs a single "methylation score" for each region, along with the number of cytosines (e.g. 
entries in infile.bed) that were used to calculate the score. The score can be calculated in two
different ways:

(1) weighted methylation (DEFAULT)
Calculates the weighted mean of the fraction of methylated reads overlapping each cytosine. This
method gives more weight to more deeply sequenced cytosines within a region.
(2) fraction of cytosines with greater than "cutoff" fraction of reads methylated
Counts the number of cytosines within a region with more than "cutoff" overlapping reads methylated,
and divides this number by the total number of cytosines with data.

In both cases, only cytosines with minimum coverage set by -x (default 5) will be included in the
analysis. Switch between methods using -m, and set the cutoff fraction methylated for method 2 using -c.

User-specified options:
Required arguments:
	-i infile : name of input methylation file in .bed format (chr start end num_unmethyl num_methyl percent_methyl)
	-r regions : file containing regions to summarize methylation over (chr start end name)
	-o outdir : name for output file
Additional arguments:
	-x : minimum reads at a particular cytosine; sites with fewer than this number are censored - default 5
	-m : method to use (1 = weighted mean, 2 = frac cytosines meeting methylation cutoff set by -c) - default 1
	-c : cutoff for minimum fraction of methylated reads at each cytosine for cytosine to be called methylated if using method 2 - default 0.5
	-s : location of helper script summarizeMethylation.py - default same as location of this script
Options:
	-a : include intervals in infile with no valid Cs (based on -x) in the output file
	-O : omit the methylation level of the 1bp (odd) or 2bp (even) at center of each interval
	-h : prints this version and usage information
	
Required installed on user PATH:
	- bedtools

Required helper scripts (if not in same dir as this script, specify location with -s):
	- summarizeMethylation.py by Colette L Picard


------------------------------------------------------------------------------------
EOF

# Clean up on exit
cleanUp () {
	# remove any tempfiles that were created (-f will ignore missing files instead of throwing error)
	rm -f "$meFixed"
	rm -f "$intersect"
	exit 0
}
 
trap cleanUp EXIT

# Initiate environment; regions, infile, and outfile must be supplied by user, all others
# are set to default values

regions=""			# regions filename
infile=""			# input filename
minReads=5			# sites with fewer than minReads total reads are censored for this analysis
method=1			# methylation summarization method (1=weighted mean, 2=cutoff), default = 1
cutoff=0.5			# cutoff to use if method == 2 (default = 0.5)
outfile=""			# output filename
omit=false
outputall=false
sumMe=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )	# location of summarizeMethylation.py (default same as this script)
 
while getopts "i:r:o:c:m:x:s:aOh" opt; do
	case $opt in
		i)	# input filename (required)
			infile="$OPTARG"
			;;
		r)	# regions filename (required)
			regions="$OPTARG"
			;;
		o)	# output filename (required)
			outfile="$OPTARG"
			;;
		c)	# cutoff to use if method == 2 (ignored if method == 1)
			cutoff="$OPTARG"
			;;
		m)	# method to use (1 == weighted mean, 2 == cutoff, see above)
			method="$OPTARG"
			;;
		x)	# sites with fewer than minReads total reads are censored for this analysis
			minReads="$OPTARG"
			;;
		s)	# location of summarizeMethylation.py (default current dir)
			sumMe="$OPTARG"
			;;
		a)	# print usage and version information to stdout and exit
			outputall=true
			;;
		O)	# print usage and version information to stdout and exit
			omit=true
			;;
		h)	# print usage and version information to stdout and exit
			echo "$usage"
			exit 0
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			exit 1
			;;
	esac
done

if [[ $# -eq 0 ]]; then
	printf "%s\n" "$usage"
	exit 0
fi

# check all necessary fields submitted
if [ -z "$infile" ]; then
	echo "Please input methylation file (.bed) using -i flag (e.g. -i infile.bed)"
	exit 1
elif [ -z "$outfile" ]; then
	echo "Please give output file name using -o flag (e.g. -o outfilename.bed)"
	exit 1
elif [ -z "$regions" ]; then
	echo "Please input regions file (.bed) using -r flag (e.g. -r regions.bed)"
	exit 1
fi

# output summary of run
echo "Summarizing methylation within genomic regions in '$regions' using methylation file '$infile'."
if [ "$method" -eq "1" ]; then
	echo "Using weighted mean. Only sites with >$minReads reads used. Outputting results to '$outfile'."
elif [ "$method" -eq "2" ]; then
	echo "Using cutoff of $cutoff. Only sites with >$minReads reads used. Outputting results to '$outfile'."
else
	echo "Invalid method: $method. Supported methods are: 1 = weighted mean, 2 = cutoff."
	exit 1
fi

# (0) Preprocessing:
# first need to convert all occurrences of "chloroplast" to ChrC and "mitochondria" to ChrM
# and make sure all occurrences of chr are capitalized to Chr for consistency (intersectBed
# will not recognize two sites as overlapping if one has chrX and other has ChrX)
echo "Preprocessing methylation input file..."

# Store preprocessed methylation data in tempfile
if [ -f "$infile" ]; then
	meFixed=$(mktemp temp.XXXXXXXXXX) || { echo "Failed to create temp file"; exit 1; }
	perl -p -w -e 's/chrc/ChrC/; s/chrm/ChrM/; s/chloroplast/ChrC/; s/mitochondria/ChrM/; s/chr/Chr/' "$infile" | awk -F$'\t' -v m="$minReads" '$4+$5 > m' > "$meFixed"
else
	echo "Methylation file '$infile' could not be read."
	exit 1
fi

# (1) Use intersectBed to match each site in methylation.bed with the corresponding region;
# store results in another tempfile
intersect=$(mktemp temp.XXXXXXXXXX) || { echo "Failed to create temp file"; exit 1; }
if [ -f "$regions" ]; then
	echo "Step 1: intersectBed running..."
	echo "command is: bedtools intersect -a \"$meFixed\" -b \"$regions\" -wb > \"$intersect\""
	bedtools intersect -a "$meFixed" -b "$regions" -wb > "$intersect"
else
	echo "Region file '$regions' could not be read."
	exit 1
fi

# (2) Run summarizeMethylation.py using supplied parameters
echo "Step 2: summarizing methylation each region..."
outfiletemp=$(mktemp temp.XXXXXXXXXX) || { echo "Failed to create temp file"; exit 1; }
[ "$omit" = "true" ] && omitstr=" -omit" || omitstr=""
echo "command is: "$sumMe"/summarizeMethylation.py \"$intersect\" \"$outfiletemp\" -method \"$method\" -cutoff \"$cutoff\" -minReads \"$minReads\"${omitstr}"
python3 "$sumMe"/summarizeMethylation.py "$intersect" "$outfiletemp" -method "$method" -cutoff "$cutoff" -minReads "$minReads"${omitstr}
[ "$outputall" = "true" ] && { bedtools intersect -b "$meFixed" -a "$regions" -v | awk -F$'\t' '{OFS=FS} {print $1,$2,$3,$4,"","0"}' >> "$outfiletemp"; }

# (3) sort based on chr, start, end of output file
sort -k1,1 -k2n,2 -k3n,3 "$outfiletemp" > "$outfile"

# remove tempfiles
echo "Cleaning up..."
rm "$meFixed"
rm "$intersect"
rm "$outfiletemp"
exit 0

