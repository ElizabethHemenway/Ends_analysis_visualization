#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys

# to filter sumby TE dataframes for min mC
# ./_filter_sumby.py /lab/solexa_gehring/elizabeth/TE_spreading/leaf_spreading/sumby/ wt_1_leaf 0.1
for arg in sys.argv:
    print(arg)
path=sys.argv[1] # working directory
name=sys.argv[2] # name of file up to mC class
min_mC=float(sys.argv[3]) # filter by this as minimum fraction mC
min_mC_str=str(min_mC)

col=['chr','start','end','feature','frac_mC','numC'] #'strand', 'DMR_chr','DMR_start','DMR_end']
tab="\t"

CHH=pd.read_csv(path+name+'_CHH_min5_sumby_TE_fragments.bed', sep=tab, header=None, names=col)
pass_CHH=CHH[CHH['frac_mC']>=min_mC]
pass_CHH.to_csv(path+name+'_CHH_min5_sumby_TE_fragments_greaterthan_'+min_mC_str+'.bed', sep=tab, header=None, index=None)

CHG=pd.read_csv(path+name+'_CHG_min5_sumby_TE_fragments.bed', sep=tab, header=None, names=col)
pass_CHG=CHG[CHG['frac_mC']>=min_mC]
pass_CHG.to_csv(path+name+'_CHG_min5_sumby_TE_fragments_greaterthan_'+min_mC_str+'.bed', sep=tab, header=None, index=None)

CG=pd.read_csv(path+name+'_CG_min5_sumby_TE_fragments.bed', sep=tab, header=None, names=col)
pass_CG=CG[CG['frac_mC']>=min_mC]
pass_CG.to_csv(path+name+'_CG_min5_sumby_TE_fragments_greaterthan_'+min_mC_str+'.bed', sep=tab, header=None, index=None)

print("done!")