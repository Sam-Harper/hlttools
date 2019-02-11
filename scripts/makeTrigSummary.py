#! /usr/bin/env python

import os
import argparse
import json
import md5
import sys
from FWCore.PythonUtilities.LumiList import LumiList
from timeit import default_timer as timer


class HLTPath :
    def __init__(self,name):
        self.name = str(name)
        self.lumi_prescaled = 0.
        self.lumi_active = 0.
        self.first_run = 999999
        self.last_run = 0
    def fill(self,lumi_prescaled,lumi_active,run):
        self.lumi_prescaled+=lumi_prescaled
        self.lumi_active+=lumi_active
        if run < self.first_run: self.first_run = run
        if run > self.last_run: self.last_run = run


class L1PreScaleCache :
    """
    HLT paths have many seeds and many paths share the same seed. It is time consuming to workout
    the lowest unprescaled L1 seed and the same information is accessed again so this caches it
    For each hash of a set of seeds, it stores the lowest prescale factor (with zero counting as 
    max int as that is effectively what it is) for each prescale column so it can be quickly accessed
    """
    def __init__(self,nr_ps_col):
        self.cache = {}
        self.nr_ps_col = nr_ps_col
        
    def lowest_l1_prescale(self,l1_seeds,l1_ps_tbl,ps_col):
        seeds_md5  = md5.new(l1_seeds).digest()
        if seeds_md5 not in self.cache:
            self.cache[seeds_md5] = [0]*self.nr_ps_col
        seed_cache = self.cache[seeds_md5]
        if seed_cache[ps_col]==0:
            lowest_l1_ps = int(sys.maxint)
            seeds_split = l1_seeds.split(" OR ") 
            #iterating through the list backwards is faster as seeds tend to be written in 
            #assending thresholds which means the last ones are more likely to be unprescaled
            #and thus the function can end early
            seeds_split.reverse() 
            ps_index = int(ps_col)+2
            for seed in seeds_split:
                seed = seed.rstrip().lstrip()
                try:
                    l1_ps = int(get_l1_prescales(l1_ps_tbl,seed)[int(ps_index)])
                except TypeError:
                    print "seed error not found ",seed
                    print "ps table",l1_ps_tbl
                    l1_ps = 0
                if l1_ps !=0 and l1_ps < lowest_l1_ps: lowest_l1_ps = l1_ps
                if lowest_l1_ps == 1: break #1 is the lowest, no need to keep looking
            seed_cache[ps_col] = int(lowest_l1_ps)
        return seed_cache[ps_col]
            

def uncompact_list(compact_list):
    
    new_list = []
    for range_ in compact_list:
        for x in range(range_[0],range_[1]+1):
            new_list.append(x)
    new_list.sort()
    return new_list
    
def get_ps_col(ps_cols,lumi):
    for ps_col in ps_cols:
        if ps_col <0: continue
        for lumi_range in ps_cols[ps_col]:
            if lumi>=lumi_range[0] and lumi<=lumi_range[1]:
                return ps_col
    #for some lumisections, wbm doesnt record them/have the prescale
    #this is always at the start/end of a run
    #so we guess by taking the column of the lumi before/after it
    #we brute force it, all possible lumi sections are written and
    #we take the first/last ones
    all_lumis = []
    for ps_col in ps_cols:
        if ps_col>=0: all_lumis.extend(uncompact_list(ps_cols[ps_col]))
    all_lumis.sort()
    if all_lumis:
        if lumi < all_lumis[0]: return get_ps_col(ps_cols,all_lumis[0])
        elif lumi > all_lumis[-1]: return get_ps_col(ps_cols,all_lumis[-1])

    return -2

def make_hlt_ps_dict(ps_tbl):
    ps_dict = {}
    for line in ps_tbl:
        pathname =get_pathname_from_ps_tbl(line[1])
        if pathname in ps_dict:
            raise RuntimeError("pathname {} already exists in the ps dict".format(pathname))
        ps_dict[pathname] = line
    return ps_dict

def get_hlt_prescales(ps_tbl,pathname):
    for line in ps_tbl:
        if get_pathname_from_ps_tbl(line[1]) == pathname:
            return line
    return None

def get_nr_ps_col(ps_tbl):
    """given a prescale table, returns how many columns it has
    does this by the length of the header
    detects if its L1 or HLT by the name of the second entry of the header which effects 
    the translation of header length to number of columns
    """
    try:
        if ps_tbl[0][1]=="L1 Algo Name":
            offset = -2
        elif ps_tbl[0][1]=="HLT Path Name":
            offset = -3
        else:
            raise RuntimeError("prescale table header not recognised as L1 or HLT".format(ps_tbl[0]))

        nr_col = len(ps_tbl[0])+offset
        return nr_col
    except IndexError:
        return 0
        
def get_l1_prescales(l1_ps_tbl,l1_seed):
    for line in l1_ps_tbl:
        if line[1] == l1_seed:
            return line
    return None
            
def get_inst_lumi(run,lumi,pu_data):
    pu_run_data = pu_data[run]
    for entry in pu_run_data:
        if lumi==entry[0]: return float(entry[1])
    return 0.

def get_lowest_l1_prescale(l1_seeds,l1_ps_tbl,ps_col,lowest_l1_ps_cache={}):
    """function doesnt handle ANDed seeds"""
    if l1_seeds.find(" AND ")!=-1 and l1_seeds.find("L1_FirstCollisionInTrain AND (NOT L1_FirstCollisionInOrbit)")==-1:
        raise RuntimeError("function can not handle AND in the seeds: {}".format(l1_seeds))

    seeds_md5  = md5.new(l1_seeds).digest()
    if seeds_md5 not in lowest_l1_ps_cache:
        lowest_l1_ps = 9999999999
        seeds_split = l1_seeds.split(" OR ") 
        ps_index = int(ps_col)+2
        for seed in seeds_split:
            seed = seed.rstrip().lstrip()
            try:
                l1_ps = int(get_l1_prescales(l1_ps_tbl,seed)[int(ps_index)])
            except TypeError:
                print "seed error not found ",seed
                print "ps table",l1_ps_tbl
                l1_ps = 0
            if l1_ps !=0 and l1_ps < lowest_l1_ps: lowest_l1_ps = l1_ps
        lowest_l1_ps_cache[seeds_md5] = int(lowest_l1_ps)
        
    return lowest_l1_ps_cache[seeds_md5]


def process_path(hltpath,run,lumi,ps_col,l1_ps_tbl,hlt_ps_dict,pu_data,lowest_l1_ps_cache):
   
    hlt_path_prescales = hlt_ps_dict[hltpath.name]
   
    if int(ps_col) < 0: print "ps column is <0",run,lumi,ps_col  
    ps_index = int(ps_col)+2
    #hlt_path_prescales has the L1 seeds as the last entry
    l1_prescale = lowest_l1_ps_cache.lowest_l1_prescale(l1_seeds=hlt_path_prescales[-1],l1_ps_tbl=l1_ps_tbl,ps_col=int(ps_col))
    hlt_prescale = int(hlt_path_prescales[ps_index])
    inst_lumi = get_inst_lumi(run=run,lumi=lumi,pu_data=pu_data)
    if hlt_prescale!=0:
        hltpath.fill(inst_lumi/hlt_prescale/l1_prescale,inst_lumi,int(run))
       

def rm_hlt_version(name):
    version_start = name.rfind("_v")
    if version_start == -1: 
        return name
    else:
        return name[:version_start+2]    


def get_pathname_from_ps_tbl(entry):
    hlt_path = entry.split()[0]
    return rm_hlt_version(hlt_path)

def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start

def is_physics_menu(hlt_menu):
    return hlt_menu.find("/cdaq/physics/Run2016/25ns") ==0 or hlt_menu.find("/cdaq/physics/Run2017/2e34/") ==0

def get_hlt_menu_version(hlt_menu):
    if not is_physics_menu(hlt_menu):
#        print "error, hlt menu name",hlt_menu," not a 2016/2017 collisions menu"
        return ""
    version = hlt_menu.split("/")[5]
    pos_minor = find_nth(version,".",2)
    if pos_minor>=0: version = version[:pos_minor]
    return version




#they are anded together
def combine_grls(grl1,grl2):
    lumis1 = LumiList(compactList=grl1)
    lumis2 = LumiList(compactList=grl2)
    
    new_lumis = lumis1 & lumis2
#    print new_lumis.compactList
    return new_lumis.compactList
        



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='prints out runs with HLT')
    parser.add_argument('--hlt_data',required=True,help='hlt data json')
    parser.add_argument('--l1_ps_data',required=True,help='l1 prescale data')
    parser.add_argument('--hlt_ps_data',required=True,help='hlt prescale data')
    parser.add_argument('--grl',required=True,help='good run list')
    parser.add_argument('--pu_data',required=True,help='pileup json')
    parser.add_argument('--datasets',required=True,help='dataset def')
    parser.add_argument('--min_runnr',default=0,type=int,help='min run number')
    parser.add_argument('--max_runnr',default=999999,type=int,help='max run number')
    parser.add_argument('--era',default="2016",help='era')
    parser.add_argument('--out',default="out.json",help='output json file')
    args = parser.parse_args()
    

    good_lumis = LumiList(filename=args.grl)

    path_to_dataset = {}
    with open(args.datasets) as f:
        path_to_dataset = json.load(f)

    hlt_data = {}
    with open(args.hlt_data) as f:
        hlt_data = json.load(f)

    l1_ps_data = {}
    with open(args.l1_ps_data) as f:
        l1_ps_data = json.load(f)

    hlt_ps_data = {}
    with open(args.hlt_ps_data) as f:
        hlt_ps_data = json.load(f)

    pu_data = {}
    with open(args.pu_data) as f:
        pu_data = json.load(f)

    runs = hlt_data.keys()
    runs.sort()

    hlt_menu_runs = {}

    hlt_paths = {}

    lowest_l1_prescales = {}
    print "loaded all data, starting processing runs"
    start_time = timer()
    nr_runs = len(runs)
    for run_indx,run in enumerate(runs):
        print "processing run {} {} / {} time: {:.1f}s".format(run,run_indx,nr_runs,timer()-start_time)
        if int(run) < args.min_runnr or int(run) > args.max_runnr: continue
            
        if not good_lumis.contains(int(run)): continue
        run_hlt_data = hlt_data[run]
        hlt_menu = run_hlt_data["hlt_menu"]
        trig_mode = run_hlt_data["trig_mode"]

        l1_ps_tbl = l1_ps_data[trig_mode]
        hlt_ps_tbl = hlt_ps_data[hlt_menu]
        hlt_ps_dict = make_hlt_ps_dict(hlt_ps_tbl)
        
        if trig_mode not in lowest_l1_prescales:
            lowest_l1_prescales[trig_mode] = L1PreScaleCache(get_nr_ps_col(l1_ps_tbl))
        runs_lowest_l1_ps = lowest_l1_prescales[trig_mode]
        

        good_lumis_compact = good_lumis.getCompactList()[run]
        good_lumis_unpacked = uncompact_list(good_lumis_compact)
        ps_cols = run_hlt_data['ps_cols']


        for lumi in good_lumis_unpacked:
            lowest_l1_ps_cache = {}
            ps_col = get_ps_col(ps_cols,lumi)

            for line in hlt_ps_tbl:
                hlt_pathname = get_pathname_from_ps_tbl(line[1]) 
                if hlt_pathname.find("HLT_")==0 and (hlt_pathname.find("Ele")!=-1 or hlt_pathname.find("Pho")!=-1 or hlt_pathname.find("pho")!=-1 or hlt_pathname.find("SC")!=-1):
                   # if hlt_pathname.find("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_")!=0: continue
                    if hlt_pathname not in hlt_paths:
                        hlt_paths[hlt_pathname]=HLTPath(hlt_pathname)
                    hlt_path = hlt_paths[hlt_pathname]

                    process_path(hlt_path,run,lumi,ps_col,l1_ps_tbl,hlt_ps_dict,pu_data,runs_lowest_l1_ps)

    hlt_path_names = hlt_paths.keys()
    hlt_path_names.sort()
    datasets = {}
    for hlt_path  in hlt_path_names:
        dataset = path_to_dataset.get(hlt_path,"NotFound")
        if dataset not in datasets:
            datasets[dataset] = []
        datasets[dataset].append(hlt_path)
    for dataset in datasets:
        datasets[dataset].sort()
    
    dataset_names = datasets.keys()
    dataset_names.sort()
    
    print "| path | active lumi (fb-1) | effective lumi (fb-1) | first run | last run | menu version | dataset | details | "
    for dataset in dataset_names:
        for hlt_path_name in datasets[dataset]:
#    for hlt_path_name in hlt_path_names:
            print hlt_path_name
            hlt_path = hlt_paths[hlt_path_name]
            try:
                hlt_version = get_hlt_menu_version(hlt_data[str(hlt_path.first_run)]["hlt_menu"])
            except KeyError:
                hlt_version = "n/a"
            ps_lumi = hlt_path.lumi_prescaled/1.0E9
            ps_lumi_format = '{:.2f}' if ps_lumi >= 0.1 else '{:.4f}'
            link_name = hlt_path_name[0:min(len(hlt_path_name),31)]
            print ("| {} | {:.2f} | "+ps_lumi_format+" | {} | {} | {} | {} | [[https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTPathDetails#{}][details]] |").format(hlt_path.name,hlt_path.lumi_active/1.0E9,hlt_path.lumi_prescaled/1.0E9,hlt_path.first_run,hlt_path.last_run,hlt_version,dataset,link_name)
    
    hlt_json  = {}
    for path in hlt_paths:
        hlt_path = hlt_paths[path]
        hlt_json[path]={}
        hlt_json[path][args.era]={}
        hlt_json[path][args.era]["lumi_active"] = hlt_path.lumi_active/1.0E9
        hlt_json[path][args.era]["lumi_effective"] = hlt_path.lumi_prescaled/1.0E9
        hlt_json[path][args.era]["first_run"] = hlt_path.first_run
        hlt_json[path][args.era]["last_run"] = hlt_path.last_run
        try:
            hlt_json[path][args.era]["menu_version"] = get_hlt_menu_version(hlt_data[str(hlt_path.first_run)]["hlt_menu"])
        except KeyError:
            hlt_json[path][args.era]["menu_version"] = "none"
        try:
            hlt_json[path][args.era]["dataset"] = path_to_dataset[path]
        except KeyError:
            hlt_json[path][args.era]["dataset"] = "NotFound"
    with open(args.out,'w') as f:
        json.dump(hlt_json,f)
