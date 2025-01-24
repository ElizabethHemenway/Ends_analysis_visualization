#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys

# after clustering, make a bed file of each cluster for browsing in IGV and downstream analysis

for arg in sys.argv:
    print(arg)
path=sys.argv[1] # where outputs should go
clust=sys.argv[2] # name of cluster .csv with path
regions=sys.argv[3] # name of regions.bed file with path
output=sys.argv[4] # sample prefix (ex. r3pat_minus_Colpat_r3TEs_CG)

bed_columns = ['chr','start','end']
big_bed_col=['chr','start','end','feature','score','strand'] #'strand', 'DMR_chr','DMR_start','DMR_end']
tab="\t"


# after clustering in R, read back in cluster assignments to make bed files. 
clusters=pd.read_csv(clust, header=0, names=['feature','cluster'])
TEs=pd.read_csv(regions, sep=tab, header=None, names=big_bed_col)

def make_bed (df, saveas=False):
    bed=pd.DataFrame()
    bed['chr']=df['chr']
    bed['start']=df['start']
    bed['end']=df['end']
    bed['feature']=df["feature"]
    #bed['strand']=df["strand"]
    
    
    if saveas is not False:
        bed.to_csv(saveas, sep="\t", index=None, header=None)
        print('Bed file saved here: '+saveas)
        return bed
    else:
        return bed
        
clustTE=TEs.merge(clusters, on='feature', how="left")
print(clustTE)

clusterlist=clusters['cluster'].unique().astype('int')
for c in clusterlist:
	clustx=clustTE[clustTE['cluster']==c]
	c_str=c.astype('str')
	make_bed(clustx, output+"_cluster"+c_str+"_browse.bed")


print("done!")