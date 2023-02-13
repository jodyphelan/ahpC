#! /usr/bin/env python

# Load useful libraries
import json
from collections import defaultdict
import argparse
import os
from tqdm import tqdm
import sys
import csv
import pathogenprofiler as pp
import tbprofiler
import math
import numpy as np


lineages = [
    "lineage1","lineage2","lineage3","lineage4",
    "lineage5","lineage6","Other"
]

subsetB_variants = set([
    ('fabG1','c.-15C>T'),
    ('inhA','c.-154G>A'),
    ('fabG1','c.-8T>C'),
    ('fabG1','c.-8T>A'),
    ('fabG1','c.-16A>G')
])

subsetC_variants = set([
    ('katG','p.Ser315Thr'),
    ('katG','p.Ser315Asn')
])

variant_exclusion = [
    ('katG','p.Arg463Leu'),
    ('katG','p.Val469Leu'),
    ('ahpC','c.-88G>A')
]

import math
class Tab:
    """Class to store 2x2 table and calculate different metrics"""
    def __init__(self):
        self.tp = 0
        self.fp = 0
        self.tn = 0
        self.fn = 0
    def sens(self):
        """Return Sensitivity"""
        return self.tp/(self.tp+self.fn)
    def spec(self):
        """return Specificity"""
        return self.tn/(self.tn+self.fp)
    def prev(self):
        """Return prevalence"""
        return (self.tp+self.fn)/(self.tn+self.fp+self.tp+self.fn)
    def ppv(self):
        """Return the PPV and 95% confidence intervals"""
        p = self.prev()
        Sp = self.spec()
        Se = self.sens()
        n1 = self.tp + self.fn
        n0 = self.tn + self.fp
        
        if (self.tp+self.fp)==0:
            return (np.nan,np.nan,np.nan)
        
        ppv = (Se * p)/((Se * p) + ((1-Sp) * (1-p)))

        if self.tp>0 and self.fp==0:
            return (ppv,np.nan,np.nan)
        if self.tp==0 and self.fp>0:
            return (0,np.nan,np.nan)

        logit_ppv = math.log((Se * p)/((1 - Sp)*(1 - p)))
        logit_ppv_var = ((1 - Se)/Se) * (1/n1) + (Sp/(1-Sp)) * (1/n0)
        logit_ppv_lb = (math.e**(logit_ppv - 1.96*math.sqrt(logit_ppv_var)))/(1+(math.e**(logit_ppv - 1.96*math.sqrt(logit_ppv_var))))
        logit_ppv_up = (math.e**(logit_ppv + 1.96*math.sqrt(logit_ppv_var)))/(1+(math.e**(logit_ppv + 1.96*math.sqrt(logit_ppv_var))))

        
        return (ppv,logit_ppv_lb,logit_ppv_up)
         

def main(args):

    dst = {}
    for row in csv.DictReader(open(args.meta)):
        if row['dst']!="NA":
            dst[row["wgs_id"]] = int(row['dst'])
    
    sample2vars = defaultdict(set)
    var2samples = defaultdict(set)
    lin2samples = defaultdict(set)
    var2genome_pos = {}

    def num_gene_variants(s,gene):
        """Return number of variants in a gene for a sample"""
        return len([d for d in sample2vars[s] if d[0]==gene])

    def has_katG_var(s):
        """Check if sample has a katG variant"""
        if "katG" in [d[0] for d in sample2vars[s]]:
            return True
        else:
            return False

    def var_in_lineage(var,samps):
        """Check if variant is in lineage"""
        res = {}
        for lin in lin2samples:
            samps_in_lin = set([s for s in lin2samples[lin] if s in samps])
            res[lin] = {}
            res[lin]["with"] = len(samps_in_lin.intersection(var2samples[var]))
            res[lin]["without"] = len(samps_in_lin) - res[lin]['with']
            res[lin]["%"] = res[lin]['with']/len(samps_in_lin)*100 if len(samps_in_lin) else np.nan
        return res

    def get_stats(var,samps):
        """Create table instance and poplate cells with numbers"""
        tab = Tab()
        for s in samps:
            if var in sample2vars[s] and dst[s]==1:
                tab.tp+=1
            elif var in sample2vars[s] and dst[s]==0:
                tab.fp+=1
            elif var not in sample2vars[s] and dst[s]==1:
                tab.fn+=1
            elif var not in sample2vars[s] and dst[s]==0:
                tab.tn+=1
            else:
                raise ValueError

        return tab

    for row in csv.DictReader(open(args.variants)):
        
        # This section skips loading variants on several conditions
        if ";" in row['sublin']: continue # mixed sample
        if row['sample_id'] not in dst: continue # samples without dst
        if float(row['freq'])<args.af: continue # minority variants (70% default)
        if row['type']=="synonymous_variant": continue # synonymous variants
        if row['gene'] not in ('katG','fabG1','inhA','ahpC'): continue # only relevant genes
        key = (row['gene'],row['change'])
        if key in variant_exclusion: continue # excludes variants defined as neutral
        
        sample2vars[row['sample_id']].add(key)
        var2samples[key].add(row['sample_id'])
        lin = row['sublin'].split(".")[0]
        if lin not in lineages:
            lin = "Other"
        lin2samples[lin].add(row['sample_id'])
        var2genome_pos[key] = int(row['genome_pos'])

    def get_result_row(var,samps):
        tab = get_stats(var,samps)
        r_samples = [s for s in samps if dst[s]==1 and s in var2samples[var]]
        s_samples = [s for s in samps if dst[s]==0 and s in var2samples[var]]
        if var==('ahpC','c.-48G>A'):
            print(vars(tab))
            print(r_samples)
        res = {
            "gene": var[0],
            "change":var[1],
            "genome_pos":var2genome_pos[var],
            "r_samples": len(r_samples),
            "s_samples": len(s_samples),
            "r_samples_with_katG": len([s for s in r_samples if has_katG_var(s)==True]),
            "r_samples_without_katG": len([s for s in r_samples if has_katG_var(s)==False]),
            "s_samples_with_katG": len([s for s in s_samples if has_katG_var(s)==True]),
            "s_samples_without_katG": len([s for s in s_samples if has_katG_var(s)==False]),
            "ppv": "%.2f" % (tab.ppv()[0]*100),
            "ppv_lb": "%.2f" % (tab.ppv()[1]*100),
            "ppv_up": "%.2f" % (tab.ppv()[2]*100)
        }
        lin_info = var_in_lineage(var,samps)

        for lin in lineages:
            res[lin+"_%"] = "%.2f" % lin_info[lin]['%']
            res[lin+"_with"] = "%.2f" % lin_info[lin]['with']
            res[lin+"_without"] = "%.2f" % lin_info[lin]['without']
        return res

    rows = []
    # variants_to_test = [('ahpC','p.Tyr34Cys')]
    # variants_to_test = [var for var in var2samples if var[0]=="ahpC"]
    variants_to_test = sorted([var for var in var2samples if var2genome_pos[var] in range(2726101,2726193)],key=lambda x:var2genome_pos[x])
    for var in tqdm(variants_to_test):
        subsetA = [s for s in list(sample2vars) if num_gene_variants(s,"ahpC")<=1]
        subsetB = [s for s in subsetA if len(sample2vars[s].intersection(subsetB_variants))==0]
        subsetC = [s for s in subsetB if len(sample2vars[s].intersection(subsetC_variants))==0]

        row = get_result_row(var,subsetA)
        row.update({"subset":"A"})
        rows.append(row)
        row = get_result_row(var,subsetB)
        row.update({"subset":"B"})
        rows.append(row)
        row = get_result_row(var,subsetC)
        row.update({"subset":"C"})
        rows.append(row)

    with open(args.out,"w") as O:
        writer = csv.DictWriter(O,fieldnames = list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

# Set up the parser
parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--meta',type=str,help='File with samples',required = True)
parser.add_argument('--variants',type=str,help='Directory containing results',required = True)
parser.add_argument('--out',type=str,default="out.csv",help='Directory containing results')
parser.add_argument('--af',default=0.75,type=float,help='Directory containing results')
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
