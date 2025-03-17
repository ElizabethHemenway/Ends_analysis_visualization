#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --job-name=TEs     # friendly name for job.
#SBATCH --nodes=1                      # ensure cpus are on one node
#SBATCH --ntasks=1                     # run a single task
#SBATCH --cpus-per-task=2             # number of cpus/threads requested.
#SBATCH --mem=42gb                      # memory requested.
#SBATCH --partition=20                 # partition (queue) to use

# 01/24/2025 - current workflow for TE spreading analysis
#prep

path="/lab/solexa_gehring/elizabeth/ros1_endo_review_analysis/ends_hmaps/"
script_path="/lab/solexa_gehring/elizabeth/TE_spreading/endo_spreading/scripts/TE_spreading_mC/"
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
		$script_path'sumByFeature.sh' -i "$file" -r "$region" -o $path$output -m 1 -x 5 
	done
done

# filter based on WT value, needs to be methylated TE. 
thresh=0.1

$script_path"_filter_sumby.py" $path"sumby_out/" $wt1 $thresh
$script_path"_filter_sumby.py" $path"sumby_out/" $wt2 $thresh
$script_path"_filter_sumby.py" $path"sumby_out/" $wt3 $thresh

cd $path"sumby_out/"
#get TE fragments that met threshold ready for ends analysis.

CG="_CG_min5_sumby_TE_fragments_greaterthan"$thresh".bed"
CHG="_CHG_min5_sumby_TE_fragments_greaterthan_"$thresh".bed"
CHH="_CHH_min5_sumby_TE_fragments_greaterthan_"$thresh".bed"

bedtools intersect -u -f 1 -a $wt1$CG -b $wt1$CHG > wt_1_TE_fragments_methylated_CG_HG.bed
bedtools intersect -u -f 1 -a wt_1_TE_fragments_methylated_CG_HG.bed -b $wt1$CHH > wt_1_TE_fragments_methylated.bed

bedtools intersect -u -f 1 -a $wt2$CG -b $wt2$CHG > wt_2_TE_fragments_methylated_CG_HG.bed
bedtools intersect -u -f 1 -a wt_2_TE_fragments_methylated_CG_HG.bed -b $wt2$CHH > wt_2_TE_fragments_methylated.bed

bedtools intersect -u -f 1 -a $wt3$CG -b $wt3$CHG > wt_3_TE_fragments_methylated_CG_HG.bed
bedtools intersect -u -f 1 -a wt_3_TE_fragments_methylated_CG_HG.bed -b $wt3$CHH > wt_3_TE_fragments_methylated.bed

bedtools intersect -u -f 1 -a wt_1_TE_fragments_methylated.bed -b wt_2_TE_fragments_methylated.bed > wt_12_TE_fragments_methylated.bed
bedtools intersect -u -f 1 -a wt_12_TE_fragments_methylated.bed -b wt_3_TE_fragments_methylated.bed > wt_TE_fragments_methylated.bed

bedtools intersect -wa -f 1 -a $regions"TE_fragments.bed" -b wt_TE_fragments_methylated.bed > TE_fragments_thresh_mC_wt_allreps.bed

rm *CG_HG* 
rm *_12_*

cp TE_fragments_thresh_mC_wt_allreps.bed $regions
 
TEthresh="TE_fragments_thresh_mC_wt_allreps.bed"



r3CG="_CG_min5_sumby_TE_fragments__1kb_r3_hypergreaterthan_"$thresh".bed"
r3CHG="_CHG_min5_sumby_TE_fragments_1kb_r3_hyper_greaterthan_"$thresh".bed"
r3CHH="_CHH_min5_sumby_TE_fragments_1kb_r3_hyper_greaterthan_"$thresh".bed"

bedtools intersect -u -f 1 -a $wt1$r3CG -b $wt1$r3CHG > wt_1_TE_fragments_1kb_r3_hyper_methylated_CG_HG.bed
bedtools intersect -u -f 1 -a wt_1_TE_fragments_1kb_r3_hyper_methylated_CG_HG.bed -b $wt1$r3CHH > wt_1_TE_fragments_1kb_r3_hyper_methylated.bed

bedtools intersect -u -f 1 -a $wt2$r3CG -b $wt2$r3CHG > wt_2_TE_fragments_1kb_r3_hyper_methylated_CG_HG.bed
bedtools intersect -u -f 1 -a wt_2_TE_fragments_1kb_r3_hyper_methylated_CG_HG.bed -b $wt2$r3CHH > wt_2_TE_fragments_1kb_r3_hyper_methylated.bed

bedtools intersect -u -f 1 -a $wt3$r3CG -b $wt3$r3CHG > wt_3_TE_fragments_1kb_r3_hyper_methylated_CG_HG.bed
bedtools intersect -u -f 1 -a wt_3_TE_fragments_1kb_r3_hyper_methylated_CG_HG.bed -b $wt3$r3CHH > wt_3_TE_fragments_1kb_r3_hyper_methylated.bed

bedtools intersect -u -f 1 -a wt_1_TE_fragments_1kb_r3_hyper_methylated.bed -b wt_2_TE_fragments_1kb_r3_hyper_methylated.bed > wt_12_TE_fragments_1kb_r3_hyper_methylated.bed
bedtools intersect -u -f 1 -a wt_12_TE_fragments_1kb_r3_hyper_methylated.bed -b wt_3_TE_fragments_1kb_r3_hyper_methylated.bed > wt_TE_fragments_1kb_r3_hyper_methylated.bed

bedtools intersect -wa -f 1 -a $regions"TE_fragments_1kb_r3_hyper.bed" -b wt_TE_fragments_1kb_r3_hyper_methylated.bed > TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps.bed

rm *CG_HG* 
rm *_12_*

cp TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps.bed $regions

TEthresh="TE_fragments_thresh_mC_wt_allreps.bed"
r3TEthresh="TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps.bed"

# run ends analysis on thresh-meeting TEs using all data
script_path_2='/lab/solexa_gehring/scripts_and_pipelines/updated_scripts/'

cd $path
#echo $path
## define contexts from data folder, I can't get the find commands to work in this context. fix later. 

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
$script_path_2'ends_analysis_eah.sh'  -O 2000 -I 2000 -r $regions'TE_fragments_thresh_mC_wt_allreps.bed' -o $endsout$base'_TE_CG' -i $CGdata -n $CGnames -w 100 -V 6 -M -R 


base=TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps
#echo $endsout$base

$script_path_2"ends_analysis_eah.sh"  -O 2000 -I 2000 -r $regions'TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps.bed' -o $endsout$base'_TE_CHH' -i $CHHdata  -n $CHHnames -w 100  -V 6 -M -R 
$script_path_2"ends_analysis_eah.sh"  -O 2000 -I 2000 -r $regions'TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps.bed' -o $endsout$base'_TE_CHG' -i $CHGdata  -n $CHGnames -w 100  -V 6 -M -R 
$script_path_2'ends_analysis_eah.sh'  -O 2000 -I 2000 -r $regions'TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps.bed' -o $endsout$base'_TE_CG' -i $CGdata -n $CGnames -w 100 -V 6 -M -R 


#cluster methylated TEs by ends analysis data
$script_path"make_diff_hmaps.r" $path r3_minus_Col_threshTEs_1kb_r3_hyper_CG \
$endsout"TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps_TE_CG_r3_1_Col_spiked_CG_min5.bed_mat.txt" $endsout"TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps_TE_CG_r3_2_Col_spiked_CG_min5.bed_mat.txt" \
$endsout"TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps_TE_CG_r3_3_Col_spiked_CG_min5.bed_mat.txt" $endsout"TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps_TE_CG_wt_1_Col_spiked_CG_min5.bed_mat.txt" \
$endsout"TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps_TE_CG_wt_2_Col_spiked_CG_min5.bed_mat.txt" $endsout"TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps_TE_CG_wt_3_Col_spiked_CG_min5.bed_mat.txt" \
100

$script_path"make_diff_hmaps.r" $path r3_minus_Col_threshTEs_1kb_r3_hyper_CHG \
$endsout"TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps_TE_CHG_r3_1_Col_spiked_CHG_min5.bed_mat.txt" $endsout"TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps_TE_CHG_r3_2_Col_spiked_CHG_min5.bed_mat.txt" \
$endsout"TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps_TE_CHG_r3_3_Col_spiked_CHG_min5.bed_mat.txt" $endsout"TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps_TE_CHG_wt_1_Col_spiked_CHG_min5.bed_mat.txt" \
$endsout"TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps_TE_CHG_wt_2_Col_spiked_CHG_min5.bed_mat.txt" $endsout"TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps_TE_CHG_wt_3_Col_spiked_CHG_min5.bed_mat.txt" \
40

$script_path"make_diff_hmaps.r" $path r3_minus_Col_threshTEs_1kb_r3_hyper_CHH \
$endsout"TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps_TE_CHH_r3_1_Col_spiked_CHH_min5.bed_mat.txt" $endsout"TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps_TE_CHH_r3_2_Col_spiked_CHH_min5.bed_mat.txt" \
$endsout"TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps_TE_CHH_r3_3_Col_spiked_CHH_min5.bed_mat.txt" $endsout"TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps_TE_CHH_wt_1_Col_spiked_CHH_min5.bed_mat.txt" \
$endsout"TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps_TE_CHH_wt_2_Col_spiked_CHH_min5.bed_mat.txt" $endsout"TE_fragments_1kb_r3_hyper_thresh_mC_wt_allreps_TE_CHH_wt_3_Col_spiked_CHH_min5.bed_mat.txt" \
20



## dme TEs
script_path_2='/lab/solexa_gehring/scripts_and_pipelines/updated_scripts/'

path=/lab/solexa_gehring/elizabeth/ros1_endo_review_analysis/dmr_targets_new/dme_TE_ends/
regions=/lab/solexa_gehring/elizabeth/ros1_endo_review_analysis/dmr_targets_new/dme_features/
comp=/lab/solexa_gehring/elizabeth/ros1_endo_review_analysis/dmr_targets_new/DMRs/
data=/lab/solexa_gehring/elizabeth/ros1_endo_review_analysis/ends_TEs/data/
outpath=$path'slurmout/'

mkdir -p $outpath
cd $path

# define contexts from data folder
CHH=$(find $data -type f -name '*CHH*')
CHHn=$(find $data -type f -name '*CHH*' -execdir basename '{}' ';')

CHG=$(find $data -type f -name '*CHG*')
CHGn=$(find $data -type f -name '*CHG*' -execdir basename '{}' ';')

CG=$(find $data -type f -name '*CG*' -o -name '*CpG*')
CGn=$(find $data -type f -name '*CG*' -execdir basename '{}' ';' -o -name '*CpG*' -execdir basename '{}' ';')


CHHdata=$(set -f; echo $CHH)
CHHnames=$(set -f; echo $CHHn)
CHGdata=$(set -f; echo $CHG)
CHGnames=$(set -f; echo $CHGn)
CGdata=$(set -f; echo $CG)
CGnames=$(set -f; echo $CGn)


CGdata=$data'C24xCol_3_Col_spiked_CG_min5.bed',$data'r1xr3_1_C24_pseudo_CG_min5.bed',$data'r1xr3_3_Col_spiked_CG_min5.bed',$data'ros1_sc_1_all_CpG_min5_CGcon_pass_fixed.bed',$data'ros1_sc_2_all_CpG_min5_CGcon_pass_fixed.bed',$data'C24xCol_2_C24_pseudo_CG_min5.bed',$data'ColxC24_3_Col_spiked_CG_min5.bed',$data'wt_sc_2_all_CpG_min5_CGcon_pass_fixed.bed',$data'C24xCol_1_C24_pseudo_CG_min5.bed',$data'dme_ros1_sc_1_all_CpG_min5.bed',$data'ColxC24_1_C24_pseudo_CG_min5.bed',$data'C24xCol_2_Col_spiked_CG_min5.bed',$data'ColxC24_2_C24_pseudo_CG_min5.bed',$data'dme_endo_hs_all_CpG.bed',$data'r3xr1_3_C24_pseudo_CG_min5.bed',$data'r3xr1_2_Col_spiked_CG_min5.bed',$data'r3xr1_3_Col_spiked_CG_min5.bed',$data'ColxC24_2_Col_spiked_CG_min5.bed',$data'ColxC24_3_C24_pseudo_CG_min5.bed',$data'r3xr1_1_C24_pseudo_CG_min5.bed',$data'wt_endo_hs_all_CpG.bed',$data'r1xr3_2_C24_pseudo_CG_min5.bed',$data'r1xr3_2_Col_spiked_CG_min5.bed',$data'r1xr3_1_Col_spiked_CG_min5.bed',$data'dme_sc_1_all_CpG_min5_CGcon_pass_fixed.bed',$data'dme_ros1_sc_2_all_CpG_min5.bed',$data'ColxC24_1_Col_spiked_CG_min5.bed',$data'r1xr3_3_C24_pseudo_CG_min5.bed',$data'wt_sc_1_all_CpG_min5_CGcon_pass_fixed.bed',$data'dme_sc_2_all_CpG_min5_CGcon_pass_fixed.bed',$data'C24xCol_1_Col_spiked_CG_min5.bed',$data'C24xCol_3_C24_pseudo_CG_min5.bed',$data'r3xr1_1_Col_spiked_CG_min5.bed',$data'r3xr1_2_C24_pseudo_CG_min5.bed'
CGnames=C24xCol_3_Col_spiked_CG_min5,r1xr3_1_C24_pseudo_CG_min5,r1xr3_3_Col_spiked_CG_min5,ros1_sc_1_all_CpG_min5_CGcon_pass_fixed,ros1_sc_2_all_CpG_min5_CGcon_pass_fixed,C24xCol_2_C24_pseudo_CG_min5,ColxC24_3_Col_spiked_CG_min5,wt_sc_2_all_CpG_min5_CGcon_pass_fixed,C24xCol_1_C24_pseudo_CG_min5,dme_ros1_sc_1_all_CpG_min5,ColxC24_1_C24_pseudo_CG_min5,C24xCol_2_Col_spiked_CG_min5,ColxC24_2_C24_pseudo_CG_min5,dme_endo_hs_all_CpG,r3xr1_3_C24_pseudo_CG_min5,r3xr1_2_Col_spiked_CG_min5,r3xr1_3_Col_spiked_CG_min5,ColxC24_2_Col_spiked_CG_min5,ColxC24_3_C24_pseudo_CG_min5,r3xr1_1_C24_pseudo_CG_min5,wt_endo_hs_all_CpG,r1xr3_2_C24_pseudo_CG_min5,r1xr3_2_Col_spiked_CG_min5,r1xr3_1_Col_spiked_CG_min5,dme_sc_1_all_CpG_min5_CGcon_pass_fixed,dme_ros1_sc_2_all_CpG_min5,ColxC24_1_Col_spiked_CG_min5,r1xr3_3_C24_pseudo_CG_min5,wt_sc_1_all_CpG_min5_CGcon_pass_fixed,dme_sc_2_all_CpG_min5_CGcon_pass_fixed,C24xCol_1_Col_spiked_CG_min5,C24xCol_3_C24_pseudo_CG_min5,r3xr1_1_Col_spiked_CG_min5,r3xr1_2_C24_pseudo_CG_min5



# define regions (easier to just type in than to figure out iterate right now)
dme_tes=$regions'TEs_1kb_dme.bed'

#run jobs
$script_path_2"ends_analysis_eah.sh"  -O 2000 -I 2000 -r $dme_tes -o 'dme_tes_CHH_in2kbout2kb' -i '$CHHdata'  -n '$CHHnames' -w 100  -V 6 -M -R 
$script_path_2"ends_analysis_eah.sh"  -O 2000 -I 2000 -r $dme_tes -o 'dme_tes_CHG_in2kbout2kb' -i '$CHGdata'  -n '$CHGnames' -w 100  -V 6 -M -R 
$script_path_2"ends_analysis_eah.sh" -O 2000 -I 2000 -r $dme_tes -o dme_tes_CG_in2kbout2kb -i '$CGdata' -n '$CGnames' -w 100 -V 6 -M -R 


#pat outputs for clustering
Colpat1=$path'dme_tes_CG_in2kbout2kb_C24xCol_1_Col_spiked_CG_min5.bed_mat.txt'
Colpat2=$path'dme_tes_CG_in2kbout2kb_C24xCol_2_Col_spiked_CG_min5.bed_mat.txt' 
Colpat3=$path'dme_tes_CG_in2kbout2kb_C24xCol_3_Col_spiked_CG_min5.bed_mat.txt' 
r3pat1=$path'dme_tes_CG_in2kbout2kb_r1xr3_1_Col_spiked_CG_min5.bed_mat.txt' 
r3pat2=$path'dme_tes_CG_in2kbout2kb_r1xr3_2_Col_spiked_CG_min5.bed_mat.txt' 
r3pat3=$path'dme_tes_CG_in2kbout2kb_r1xr3_3_Col_spiked_CG_min5.bed_mat.txt'

Col1=$path'dme_tes_CG_in2kbout2kb_wt_1_Col_spiked_CG_min5.bed_mat.txt'
Col2=$path'dme_tes_CG_in2kbout2kb_wt_2_Col_spiked_CG_min5.bed_mat.txt' 
Col3=$path'dme_tes_CG_in2kbout2kb_wt_3_Col_spiked_CG_min5.bed_mat.txt' 
r31=$path'dme_tes_CG_in2kbout2kb_r3_1_Col_spiked_CG_min5.bed_mat.txt'
r32=$path'dme_tes_CG_in2kbout2kb_r3_2_Col_spiked_CG_min5.bed_mat.txt' 
r33=$path'dme_tes_CG_in2kbout2kb_r3_3_Col_spiked_CG_min5.bed_mat.txt' 


cd $regions'ends/'

$script_path'make_diff_hmaps.r' $regions dme_TEs_CG $r31 $r32 $r33 $Col1 $Col2 $Col3 50


Col1=$path'dme_tes_CHG_in2kbout2kb_wt_1_Col_spiked_CHG_min5.bed_mat.txt'
Col2=$path'dme_tes_CHG_in2kbout2kb_wt_2_Col_spiked_CHG_min5.bed_mat.txt' 
Col3=$path'dme_tes_CHG_in2kbout2kb_wt_3_Col_spiked_CHG_min5.bed_mat.txt' 
r31=$path'dme_tes_CHG_in2kbout2kb_r3_1_Col_spiked_CHG_min5.bed_mat.txt'
r32=$path'dme_tes_CHG_in2kbout2kb_r3_2_Col_spiked_CHG_min5.bed_mat.txt' 
r33=$path'dme_tes_CHG_in2kbout2kb_r3_3_Col_spiked_CHG_min5.bed_mat.txt' 
$regions'make_diff_hmaps.r' $regions'ends/' dme_TEs_CHG $r31 $r32 $r33 $Col1 $Col2 $Col3 20

Col1=$path'dme_tes_CHH_in2kbout2kb_wt_1_Col_spiked_CHH_min5.bed_mat.txt'
Col2=$path'dme_tes_CHH_in2kbout2kb_wt_2_Col_spiked_CHH_min5.bed_mat.txt' 
Col3=$path'dme_tes_CHH_in2kbout2kb_wt_3_Col_spiked_CHH_min5.bed_mat.txt' 
r31=$path'dme_tes_CHH_in2kbout2kb_r3_1_Col_spiked_CHH_min5.bed_mat.txt'
r32=$path'dme_tes_CHH_in2kbout2kb_r3_2_Col_spiked_CHH_min5.bed_mat.txt' 
r33=$path'dme_tes_CHH_in2kbout2kb_r3_3_Col_spiked_CHH_min5.bed_mat.txt' 
$regions'make_diff_hmaps.r' $regions'ends/' dme_TEs_CHH $r31 $r32 $r33 $Col1 $Col2 $Col3 10
