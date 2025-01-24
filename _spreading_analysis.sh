#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --job-name=TEs     # friendly name for job.
#SBATCH --nodes=1                      # ensure cpus are on one node
#SBATCH --ntasks=1                     # run a single task
#SBATCH --cpus-per-task=2             # number of cpus/threads requested.
#SBATCH --mem=42gb                      # memory requested.
#SBATCH --partition=20                 # partition (queue) to use

# 01/24/2025 - current workflow for TE spreading analysis.

#prep

path="/lab/solexa_gehring/elizabeth/TE_spreading/endo_spreading/"
script_path=$path"scripts/"
outpath=$path"slurmout/"
data=$path"data/"
regions=$path"regions/"

mkdir -p $outpath

m1="r3_1_Col_spiked"
m2="r3_2_Col_spiked"
m3="r3_3_Col_spiked"
wt1="wt_1_Col_spiked"
wt2="wt_2_Col_spiked"
wt3="wt_3_Col_spiked"

#sumby WT methylation over TEs
FILES=$(find $data -type f -name '*min5.bed' -name '*wt*')
for file in $FILES ; do
	#echo $file
	base=$(basename -- "$file" .bed)
	ROI=$(find $regions -name "*.bed")
	for region in $ROI ; do
		regbase=$(basename -- "$region" .bed )
		output="sumby_out/"$base"_sumby_"$regbase.bed
		#echo $path$output
		#$script_path'sumByFeature.sh' -i "$file" -r "$region" -o $path$output -m 1 -x 5 
	done
done

# filter based on WT value, needs to be methylated TE. 
thresh=0.1
mkdir -p $path"sumby_out/"

#$script_path"_filter_sumby.py" $path"sumby_out/" $wt1 $thresh
#$script_path"_filter_sumby.py" $path"sumby_out/" $wt2 $thresh
#$script_path"_filter_sumby.py" $path"sumby_out/" $wt3 $thresh

cd $path"sumby_out/"
#get TE fragments that met threshold ready for ends analysis.

CG="_CG_min5_sumby_TE_fragments_greaterthan_"$thresh".bed"
CHG="_CHG_min5_sumby_TE_fragments_greaterthan_"$thresh".bed"
CHH="_CHH_min5_sumby_TE_fragments_greaterthan_"$thresh".bed"

#bedtools intersect -u -f 1 -a $wt1$CG -b $wt1$CHG > wt_1_TE_fragments_methylated_CG_HG.bed
#bedtools intersect -u -f 1 -a wt_1_TE_fragments_methylated_CG_HG.bed -b $wt1$CHH > wt_1_TE_fragments_methylated.bed

#bedtools intersect -u -f 1 -a $wt2$CG -b $wt2$CHG > wt_2_TE_fragments_methylated_CG_HG.bed
#bedtools intersect -u -f 1 -a wt_2_TE_fragments_methylated_CG_HG.bed -b $wt2$CHH > wt_2_TE_fragments_methylated.bed

#bedtools intersect -u -f 1 -a $wt3$CG -b $wt3$CHG > wt_3_TE_fragments_methylated_CG_HG.bed
#bedtools intersect -u -f 1 -a wt_3_TE_fragments_methylated_CG_HG.bed -b $wt3$CHH > wt_3_TE_fragments_methylated.bed

#bedtools intersect -u -f 1 -a wt_1_TE_fragments_methylated.bed -b wt_2_TE_fragments_methylated.bed > wt_12_TE_fragments_methylated.bed
#bedtools intersect -u -f 1 -a wt_12_TE_fragments_methylated.bed -b wt_3_TE_fragments_methylated.bed > wt_TE_fragments_methylated.bed

#bedtools intersect -wa -f 1 -a $regions"TE_fragments.bed" -b wt_TE_fragments_methylated.bed > TE_fragments_thresh_mC_wt_allreps.bed

#rm *CG_HG* 
#rm *_12_*

#cp TE_fragments_thresh_mC_wt_allreps.bed $regions

TEthresh="TE_fragments_thresh_mC_wt_allreps.bed"

# run ends analysis on thresh-meeting TEs using all data
script_path_2='/lab/solexa_gehring/scripts_and_pipelines/updated_scripts/'

cd $path
#echo $path
## define contexts from data folder
#CG=$(find ./data/ -type f -name '*CG*' -print | awk '{gsub(/ /, ""); if (NR > 1) print prev ","; prev=$0} END {print prev}')
#CGn=$(find ./data/ -type f -name '*CG*' -print | xargs -n 1 basename | awk '{gsub(/ /, ""); if (NR > 1) print prev ","; prev=$0} END {print prev}')

CGdata=./data/wt_1_Col_spiked_CG_min5.bed,./data/wt_2_Col_spiked_CG_min5.bed,./data/wt_3_Col_spiked_CG_min5.bed,./data/r3_1_Col_spiked_CG_min5.bed,./data/r3_2_Col_spiked_CG_min5.bed,./data/r3_3_Col_spiked_CG_min5.bed
CGnames=wt_1_Col_spiked_CG_min5.bed,wt_2_Col_spiked_CG_min5.bed,wt_3_Col_spiked_CG_min5.bed,r3_1_Col_spiked_CG_min5.bed,r3_2_Col_spiked_CG_min5.bed,r3_3_Col_spiked_CG_min5.bed
CHGdata=./data/wt_1_Col_spiked_CHG_min5.bed,./data/wt_2_Col_spiked_CHG_min5.bed,./data/wt_3_Col_spiked_CHG_min5.bed,./data/r3_1_Col_spiked_CHG_min5.bed,./data/r3_2_Col_spiked_CHG_min5.bed,./data/r3_3_Col_spiked_CHG_min5.bed
CHGnames=wt_1_Col_spiked_CHG_min5.bed,wt_2_Col_spiked_CHG_min5.bed,wt_3_Col_spiked_CHG_min5.bed,r3_1_Col_spiked_CHG_min5.bed,r3_2_Col_spiked_CHG_min5.bed,r3_3_Col_spiked_CHG_min5.bed
CHHdata=./data/wt_1_Col_spiked_CHH_min5.bed,./data/wt_2_Col_spiked_CHH_min5.bed,./data/wt_3_Col_spiked_CHH_min5.bed,./data/r3_1_Col_spiked_CHH_min5.bed,./data/r3_2_Col_spiked_CHH_min5.bed,./data/r3_3_Col_spiked_CHH_min5.bed
CHHnames=wt_1_Col_spiked_CHH_min5.bed,wt_2_Col_spiked_CHH_min5.bed,wt_3_Col_spiked_CHH_min5.bed,r3_1_Col_spiked_CHH_min5.bed,r3_2_Col_spiked_CHH_min5.bed,r3_3_Col_spiked_CHH_min5.bed

endsout=$path"ends/"
mkdir -p $endsout
base=TE_fragments_thresh_mC_wt_allreps
#echo $endsout$base

$script_path_2"ends_analysis_eah.sh"  -O 2000 -I 2000 -r $regions'TE_fragments_thresh_mC_wt_allreps.bed' -o $endsout$base'_TE_CHH' -i $CHHdata  -n $CHHnames -w 100  -V 6 -M -R 
$script_path_2"ends_analysis_eah.sh"  -O 2000 -I 2000 -r $regions'TE_fragments_thresh_mC_wt_allreps.bed' -o $endsout$base'_TE_CHG' -i $CHGdata  -n $CHGnames -w 100  -V 6 -M -R 
#$script_path_2'ends_analysis_eah.sh'  -O 2000 -I 2000 -r $regions'TE_fragments_thresh_mC_wt_allreps.bed' -o $endsout$base'_TE_CG' -i $CGdata -n $CGnames -w 100 -V 6 -M -R 

#cluster methylated TEs by ends analysis data
$script_path"cluster_by_ends_analysis_v2.r" $path 3 r3_minus_Col_threshTEs_CG \
$endsout"TE_fragments_thresh_mC_wt_allreps_TE_CG_r3_1_Col_spiked_CG_min5.bed_mat.txt" $endsout"TE_fragments_thresh_mC_wt_allreps_TE_CG_r3_2_Col_spiked_CG_min5.bed_mat.txt" \
$endsout"TE_fragments_thresh_mC_wt_allreps_TE_CG_r3_3_Col_spiked_CG_min5.bed_mat.txt" $endsout"TE_fragments_thresh_mC_wt_allreps_TE_CG_wt_1_Col_spiked_CG_min5.bed_mat.txt" \
$endsout"TE_fragments_thresh_mC_wt_allreps_TE_CG_wt_2_Col_spiked_CG_min5.bed_mat.txt" $endsout"TE_fragments_thresh_mC_wt_allreps_TE_CG_wt_3_Col_spiked_CG_min5.bed_mat.txt"

$script_path"cluster_by_ends_analysis_v2.r" $path 3 r3_minus_Col_threshTEs_CHG \
$endsout"TE_fragments_thresh_mC_wt_allreps_TE_CHG_r3_1_Col_spiked_CHG_min5.bed_mat.txt" $endsout"TE_fragments_thresh_mC_wt_allreps_TE_CHG_r3_2_Col_spiked_CHG_min5.bed_mat.txt" \
$endsout"TE_fragments_thresh_mC_wt_allreps_TE_CHG_r3_3_Col_spiked_CHG_min5.bed_mat.txt" $endsout"TE_fragments_thresh_mC_wt_allreps_TE_CHG_wt_1_Col_spiked_CHG_min5.bed_mat.txt" \
$endsout"TE_fragments_thresh_mC_wt_allreps_TE_CHG_wt_2_Col_spiked_CHG_min5.bed_mat.txt" $endsout"TE_fragments_thresh_mC_wt_allreps_TE_CHG_wt_3_Col_spiked_CHG_min5.bed_mat.txt"

$script_path"cluster_by_ends_analysis_v2.r" $path 3 r3_minus_Col_threshTEs_CHH \
$endsout"TE_fragments_thresh_mC_wt_allreps_TE_CHH_r3_1_Col_spiked_CHH_min5.bed_mat.txt" $endsout"TE_fragments_thresh_mC_wt_allreps_TE_CHH_r3_2_Col_spiked_CHH_min5.bed_mat.txt" \
$endsout"TE_fragments_thresh_mC_wt_allreps_TE_CHH_r3_3_Col_spiked_CHH_min5.bed_mat.txt" $endsout"TE_fragments_thresh_mC_wt_allreps_TE_CHH_wt_1_Col_spiked_CHH_min5.bed_mat.txt" \
$endsout"TE_fragments_thresh_mC_wt_allreps_TE_CHH_wt_2_Col_spiked_CHH_min5.bed_mat.txt" $endsout"TE_fragments_thresh_mC_wt_allreps_TE_CHH_wt_3_Col_spiked_CHH_min5.bed_mat.txt"


#make bed files for browsing out of cluster outputs
$script_path"_cluster2bedfiles.py" $path r3_minus_Col_threshTEs_CG_clusters.csv $regions"TE_fragments.bed" r3_minus_Col_threshTEs_CG
$script_path"_cluster2bedfiles.py" $path r3_minus_Col_threshTEs_CHG_clusters.csv $regions"TE_fragments.bed" r3_minus_Col_threshTEs_CHG
$script_path"_cluster2bedfiles.py" $path r3_minus_Col_threshTEs_CHH_clusters.csv $regions"TE_fragments.bed" r3_minus_Col_threshTEs_CHH

echo "done!"