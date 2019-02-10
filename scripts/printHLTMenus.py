#! /usr/bin/env python

import os
import argparse
import json

def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start

def is_physics_menu(hlt_menu):
    return hlt_menu.find("/cdaq/physics/Run2016/25ns") ==0 or hlt_menu.find("/cdaq/physics/Run2017/2e34/") ==0 or hlt_menu.find("/cdaq/physics/Run2018/2e34/") ==0

def get_hlt_menu_version(hlt_menu):

    if not is_physics_menu(hlt_menu):
#        print "error, hlt menu name",hlt_menu," not a 2016/2017 collisions menu"
        return ""
    version = hlt_menu.split("/")[5]
    pos_minor = find_nth(version,".",2)
    if pos_minor>=0: version = version[:pos_minor]
    if version=="V3.0": version = "v3.0" #hack to fix a bug
    return version






if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='prints out runs with HLT')
    parser.add_argument('--hlt_data',help='hlt data json')
    parser.add_argument('--grl',help='good run list')
    parser.add_argument('--min_runnr',default=0,type=int,help='min run number')
    parser.add_argument('--max_runnr',default=999999,type=int,help='max run number')
    parser.add_argument('--short',action='store_true',help='just prints the major/minor version')
    args = parser.parse_args()


    good_lumis = {}
    with open(args.grl) as f:
        good_lumis = json.load(f)

    hlt_data = {}
    with open(args.hlt_data) as f:
        hlt_data = json.load(f)

    runs = hlt_data.keys()
    runs.sort()

    hlt_menu_runs = {}

    for run in runs:
        hlt_menu = hlt_data[run]["hlt_menu"]
        if args.short:
            hlt_menu = get_hlt_menu_version(hlt_menu)
        if hlt_menu not in hlt_menu_runs:
            hlt_menu_runs[hlt_menu] = []
        hlt_menu_runs[hlt_menu].append(int(run))

        version = get_hlt_menu_version(hlt_menu)
        if version=="v3.1" or version == "v3.2":
            print run,hlt_menu

    hlt_menus = hlt_menu_runs.keys()
    hlt_menus.sort()

    for hlt_menu in hlt_menus:
        runs = hlt_menu_runs[hlt_menu]
        runs.sort()
        if is_physics_menu(hlt_menu):
            print "{} {}-{}".format(hlt_menu,runs[0],runs[-1])

    major_version_dict = {}
    physics_menus = (x for x in hlt_menus if is_physics_menu(x))

    for hlt_menu in physics_menus:
        version = get_hlt_menu_version(hlt_menu)
        if version not in major_version_dict:
            major_version_dict[version] = [hlt_menu,min(hlt_menu_runs[hlt_menu]),max(hlt_menu_runs[hlt_menu])]
        else:
            max_run = max(hlt_menu_runs[hlt_menu])
            min_run = min(hlt_menu_runs[hlt_menu])
            if max_run > major_version_dict[version][2]: 
                major_version_dict[version][0] = hlt_menu
                major_version_dict[version][2] = max_run
            if min_run < major_version_dict[version][1]:
                major_version_dict[version][1] = min_run
    
    versions = major_version_dict.keys()
    versions.sort()
    print "| version | run range | last hlt menu |"
    for version in versions:
        print "| {} | {}-{} | {} |".format(version,major_version_dict[version][1],major_version_dict[version][2],major_version_dict[version][0])
        
