#!/usr/bin/env python
import json
import argparse


def get_active_lumi_str(path_data):
    out_str = "active lumi:"
    eras = path_data.keys()
    eras.sort()
    for era in eras:
        out_str+=" {:.2f} fb<sup>-1</sup> ({})".format(path_data[era]["lumi_active"],era)
    return out_str

def get_effective_lumi_str(path_data):
    out_str = "effective lumi:"
    eras = path_data.keys()
    eras.sort()
    for era in eras:
        eff_lumi = path_data[era]["lumi_effective"]
        if eff_lumi >= 0.1:
            out_str+=" {:.2f} fb<sup>-1</sup> ({})".format(eff_lumi,era)
        else:
            out_str+=" {:.4f} fb<sup>-1</sup> ({})".format(eff_lumi,era)
    return out_str

def get_dataset_str(path_data):
    out_str = "dataset:"
    eras = path_data.keys()
    eras.sort()
    for era in eras:
        out_str+=" {} ({})".format(path_data[era]["dataset"],era)
    return out_str

def get_run_range_str(path_data):
    out_str = "run range in menu:"
    eras = path_data.keys()
    eras.sort()
    for era in eras:
        out_str +=" {} - {} ({})".format(path_data[era]['first_run'],path_data[era]['last_run'],era)
    return out_str

def get_sel_str(sel):
    out_str = ""
    ver_str = ""
    pre_sel = ""
    menu_versions = sel.keys()
    menu_versions.sort()
    for menu_version in menu_versions:
        try:
            if sel[menu_version]==pre_sel:
                ver_str+=", "+menu_version
            else:
                if ver_str != "":
                    out_str+= "menu versions : {} <br>\n".format(ver_str)
                    out_str+= pre_sel+"\n"
                pre_sel = sel[menu_version]
                ver_str = str(menu_version)
        except KeyError:
            pass
    out_str += "menu versions : {} <br>\n".format(ver_str.replace("p","."))
    out_str += pre_sel+"\n"
    return out_str
    
def combine_dicts(all_dicts):
    """combines the various dictionaries"
    they are all keyed by path name without version
    we need to combine selection jsons and the era data dictionaries
    the era data dictionaries will be one for each era
    the selection dictionary will also be unique for now
    """
    comb_dict = {}
    for dict_ in all_dicts:
        for path in dict_:
            if path not in comb_dict:
                comb_dict[path] = {}
            for key in dict_[path]:
                if key!="selection":
                    if 'meta_data' not in comb_dict[path]:
                        comb_dict[path]['meta_data'] = {}
                    comb_dict[path]['meta_data'][key] = dict(dict_[path][key])
                else:
                    comb_dict[path][key] = dict(dict_[path][key])

    return comb_dict



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='reads in various HLT data and generates a summary twiki')
    parser.add_argument('hlt_jsons',nargs="+",help='hlt data json')
    args = parser.parse_args()
    
    all_dicts = []
    for hlt_json in args.hlt_jsons:
        with open(hlt_json) as f:
            all_dicts.append(json.load(f))

    hlt_data = combine_dicts(all_dicts)

   # print "dict "
   # print hlt_data
    
    print '''
Welcome to the E/g HLT 2016-2018 Path Information Page

This twiki is automatically generated from HLT menus & wbm info, please do not edit this twiki directly as your edits will be lost next time the auto generation script is run. Please contact E/gamma for any edits you wish to make. <br>

warning, this is still underconstruction<br> 

The information on this twiki is thought to be accurate except for the cavates below but again as its auto generated, there may be mistakes in edge cases for paths, please let us know if you find any. Additionally, the lumi numbers are approximate and should only be treated as guidence rather than analysis quality numbers. Just remember this is simply some python code which parses the HLT and tries to print what it things the HLT would do when given this config and as we have to reimpliment the logic in python, there could be bugs  <br>

Known issues: 
   * variable definations are hard coded and do not full evolve in time
      * main issue is tracking when rho correction went from variable to filter
   * paths moving dataset is only tracker per year basis and so it if moves in the year, this will not be recorded
      * only happens in 2017 for HLT_DoubleMu20_7_Mass0to30_Photon23_v which migrates from DoubleMuon -> MuonEG
   * code changes in modules may be not taken into account (although unlikely)
   * there is some disagreement on what constitutes barrel/endcap start end (eg 1.479, 1.5 , 2.5, 2.65) so to make life easier, we fuzz those boundaries so exact eta values may be slightly off,to be fixed
   * does not handle DZ, path leg combination filters
   * does not handle anything but e/gamma
   * should work for standard paths but werid complex paths may have edge cases
   * right now does not distingush all cases where there is a discontected filter, that is a filter in the path but the filters after it do not depend on it. So it means you may be requiring an object with Et>50  GeV and an object with H/E<0.15 which is different to an object with Et>50 GeV and H/E <0.15. Think we've got most of them but some remain (and those that remain are likely bugs)
   * its a little slow to render, however not sure how to make it faster without splitting it into different pages

---++ Interpreting The Twiki
This section describes how interpret the results for each trigger. 

---+++ Layout of an CMS Trigger

A CMS trigger path consists of a bunch of producers (modules which produce superclusters, id variables like sigmaIEtaIEta, isolations etc) and filters which cut on those quantities. As soon as a filter returns false (assuming it is not ignored or negated), the execution of the path stops and is marked as failing. Paths can and do share producer and filter modules. There are a relatively small number of producers shared between all E/gamma paths that require the variable it produces while there are many different filters. 

HLT Filters (which are a subset of standard EDFilters) do the following. First they take an input collection of objects, apply a selection on those objects and then put a collection of those objects passing the selection in the event. This collection of passing objects is then often used as an input for subsequent filters and is how we probably AND cuts together. For example instead of requiring 1 electron with Et>40 GeV and 1 electron with sigmaIEtaIEta<0.011 this ensures that we actually require 1 electron with Et>40 GeV and sigmaIEtaIEta<0.011, ie its the same electron passing the two cuts. Secondly the filter returns a bool which indicates whether the filter passes or fails. What this means is filter specific but in E/gamma it means the required number of objects passing the selection were found. If this filter returns false, unless the filter is ignored or negated, the execution of the path stops. An important point to realise is that the filter will write out a list of all the objects that passes its selection, even if the number of required objects is not found. So for a double electron filter requiring H/E<0.1, if there is only one electron with H/E<0.1, it will still write out that electron as passing the filter and that information will be saved in the event however as there is only 1 electron, the filter will return false and the path executation will stop with the path not accepting the event. 

The filters global decision can also be ignored or negated. Filter negation is not often used in e/gamma however the ability to ignore the pass is used to make the OR of selection requirements in the diphoton paths. An example is as follows, if you want H/E<0.1 AND ( sigmaIEtaIEta<0.11 OR a R9>0.9 ) the filter is constructed as : SigmaIEtaIEtaFilter(running on H/E input), ignore the result, allow the path to continue, R9Filter(running on H/E input), ignore the result, allow the path to continue and then a special filter which looks at the SigmaIEtaIEta and R9 filter outputs and if there a photon passing either one of them, the filter will accept. This ORing is represented by colour coding the parts of the chain that are ORed together. 

In the HLT filters exist as chains of filters, with each subsequent filter in the chain reading in the objects the previous filter has passed. Simple single object paths will just have a single filter chain while multiobject paths will have multiple chains.  This is representing in this twiki by having a seperate filter table for each chain. 




'''
    print "%TOC%"

    path_names = hlt_data.keys()
    path_names.sort()

    for path_name in path_names:
        path_data = hlt_data[path_name].get('meta_data',{})
        path_sel = hlt_data[path_name].get('selection',{})
        print "---++ "+path_name
        prefix = "   * "
#        print prefix+"recommended signal trigger:"+path_data.get('recommend','')
        print prefix+"purpose:"+path_data.get('purpose','')
        print prefix+get_active_lumi_str(path_data)
        print prefix+get_effective_lumi_str(path_data)
        print prefix+get_run_range_str(path_data)
        print prefix+get_dataset_str(path_data)
        print "\nselection:<br>"
        print get_sel_str(path_sel),

    
