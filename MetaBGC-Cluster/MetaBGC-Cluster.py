#!/usr/bin/env python3


import pandas as pd
from sklearn.cluster import DBSCAN
import json
import re
import argparse


def reads_DBSCAN(args):
    df = pd.read_csv(args.table, delimiter="\t")
    mat = df.iloc[:,1:df.shape[1]].values
    read_names = df.iloc[:,0].values
    cl = DBSCAN(eps=args.max_dist, min_samples=args.min_samples, metric="correlation", 
                algorithm="brute", n_jobs=args.threads).fit_predict(mat)
    d = dict(zip(read_names, cl.tolist()))
    with open(re.sub("(.*)\\..*", r"\1_DBSCAN.json", args.table), "w") as h:
        json.dump(d, h, indent=4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="cluster reads with DBSCAN")
    parser.add_argument("--table", type=str, required=True, 
                        help="path of tab-delimited abundance table")
    parser.add_argument("--max_dist", type=float, default=0.1, 
                        help="maximum Pearson distance between two reads "\
                             "to be in the same cluster. Default is 0.1")
    parser.add_argument("--min_samples", type=float, default=1, 
                        help="minimum number of samples required for a cluster. "\
                             "Default is 1. If min_samples > 1, noise are "\
                             "labelled as -1")
    parser.add_argument("--threads", type=int, default=1, 
                        help="number of threads for parallel jobs. Default is 1")
    args = parser.parse_args()
    reads_DBSCAN(args)
