#!/usr/bin/env python

#####################################################################################
#@author: shuowang
# This file is a component of MetaBGC (Metagenomic identifier of Biosynthetic Gene Clusters)
# (contact Francine Camacho at camachofrancine@gmail.com).
#####################################################################################

import pandas as pd
from sklearn.cluster import DBSCAN
import json
import re

def mbgccluster(table,max_dist,min_samples,cpu):
    df = pd.read_csv(table, delimiter="\t")
    mat = df.iloc[:,1:df.shape[1]].values
    read_names = df.iloc[:,0].values
    cl = DBSCAN(eps=max_dist, min_samples=min_samples, metric="correlation",
                algorithm="brute", n_jobs=cpu).fit_predict(mat)
    d = dict(zip(read_names, cl.tolist()))
    out_file = re.sub("(.*)\\..*", r"\1_DBSCAN.json", table)
    with open(out_file, "w") as h:
        json.dump(d, h, indent=4)
    return out_file
