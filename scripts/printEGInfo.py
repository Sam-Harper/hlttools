#!/usr/bin/env python

def is_egamma_path(process,path_name):
    """A path is an E/gamma path if it contains hltEgammaCandidates or hltEgammaCandidatesUnseeded in the modules it runs"""   
    path = getattr(process,path_name)
    return "hltEgammaCandidates" in path.moduleNames() or "hltEgammaCandidatesUnseeded" in path.moduleNames()

def get_nice_var_name(name):

    return {
        "hltEgammaHoverE":"H (for H/E, cone<0.14)",
        "hltEgammaGsfTrackVars:DetaSeed":"DEtaInSeed",
        "hltEgammaGsfTrackVars:Dphi":"DPhiIn",
        "hltEgammaGsfTrackVars:MissingHits":"#miss hits",
        "hltEgammaGsfTrackVars:Chi2":"trk chi2",
        "hltEgammaGsfTrackVars:OneOESuperMinusOneOP":"1/E-1/p",
        "hltEgammaEcalPFClusterIso":"ECAL Clus Iso",
        "hltEgammaHcalPFClusterIso":"HCAL Clus Iso",
        "hltEgammaEleGsfTrackIso":"Ele Trk Iso",
        "hltEgammaClusterShape:sigmaIEtaIEta5x5":"sigmaIEtaIEta (full5x5)"

    }.get(name, str(name))  

class EGammaStdCut:
    def __init__(self,op_type="",divide_by_var="E",const_term=-1.0,linear_term=-1.0,quad_term=-1.0,term_combine_op="||"):
        self.op_type = op_type
        self.divide_by_var = divide_by_var
        self.const_term = const_term
        self.linear_term = linear_term
        self.quad_term = quad_term
        self.term_combine_op = term_combine_op
    def label(self):
        label_str = "{}".format(self.op_type)
        first_term = True
        ignore_value = -1.0 if self.term_combine_op=="||" else 0.
        if self.const_term!=ignore_value:
            label_str+=" {}".format(self.const_term)
            first_term = False
        if self.linear_term!=ignore_value:
            if not first_term: label_str+=" {} ".format(self.term_combine_op)
            label_str+=" {} * {}".format(self.linear_term,self.divide_by_var)
        if self.quad_term!=ignore_value:
            if not first_term: label_str+=" {} ".format(self.term_combine_op)
            label_str+=" {} * {}^{{2}}".format(self.quad_term,self.divide_by_var)
            
        #we need to deal with the special case where its all zeros
        if first_term and ignore_value==0.:
            label_str+=" 0."

        return label_str



class EGammaFilter:
    def __init__(self,filt_mod=None):
        if(filt_mod) fill(filt_mod)
        
    def fill(filt_mod):
        self.filt_type = filt_mod.type_()
        if self.filt_type=="HLTEgammaGenericFilter":
            self.cut_params = filt_mod.parameters_()
            self.var = file_mod.varTag.value()
            #now remove the parameters which are not part of the selection from cut_params
            for x in ['saveTags','l1EGCand','candTag']: del self.cut_params[x]
        
        

class EGammaSelChain:
    def __init__(self):
        self.operator = ""
        self.divide_by_et = True
        self.eta_range = []
        self.cuts = ""

def printHLT1Photon(filt,indent=6):
    min_pt = filt.getParameter("MinPt").value()
    min_e = filt.getParameter("MinE").value()
    max_eta  = filt.getParameter("MaxEta").value()
    label_str=""
    if min_pt!=-1: label_str+="{}E_{{T}} >= {}".format(" "*indent,min_pt)
    if min_e!=-1: label_str+=" E >= {}".format(min_e)
    if max_eta!=-1: label_str+=" |eta| <= {}".format(max_eta)
    print label_str

def printHLTEgammaTriggerFilterObjectWrapper(filt,indent=6):
    print "{}L1 match NOT required".format(" "*indent)

def printHLTEgammaL1TMatchFilterRegional(filt,indent=6):
    print "{}L1 match required".format(" "*indent)

def printHLTEgammaEtFilter(filt,indent=6):
    print "{}E_{{T}} >= {} (EB) >= {} (EE)".format(" "*indent,filt.getParameter("etcutEB").value(),filt.getParameter("etcutEE").value())
    
def printHLTElectronPixelMatchFilter(filt,indent=6):
    label_str= "{}PixelMatch".format(" "*indent)
    if filt.getParameter("useS").value(): label_str+=" using S2 additional matching"
    label_str+= " required" if not filt.getParameter("pixelVeto").value() else " vetoed"
    print label_str

def printHLTEgammaGenericQuadraticFilter(filt,indent=6):
    op_str = "<=" if filt.getParameter("lessThan").value() else ">="
    et_str = "E_{T}" if filt.getParameter("useEt").value() else "E"
    var_name = get_nice_var_name(filt.getParameter("varTag").value())
    cut_barrel = EGammaStdCut(op_str,et_str,filt.getParameter("thrRegularEB").value(),filt.getParameter("thrOverEEB").value(),filt.getParameter("thrOverE2EB").value(),"+")
    cut_endcap = EGammaStdCut(op_str,et_str,filt.getParameter("thrRegularEE").value(),filt.getParameter("thrOverEEE").value(),filt.getParameter("thrOverE2EE").value(),"+")
    print "{}{} {} (EB) {} (EE)".format(" "*indent,var_name,cut_barrel.label(),cut_endcap.label())

def printHLTEgammaGenericFilter(filt,indent=6):
    op_str = "<=" if filt.getParameter("lessThan").value() else ">="
    et_str = "E_{T}" if filt.getParameter("useEt").value() else "E"
    var_name = get_nice_var_name(filt.getParameter("varTag").value())
    cut_barrel = EGammaStdCut(op_str,et_str,filt.getParameter("thrRegularEB").value(),filt.getParameter("thrOverEEB").value(),filt.getParameter("thrOverE2EB").value())
    cut_endcap = EGammaStdCut(op_str,et_str,filt.getParameter("thrRegularEE").value(),filt.getParameter("thrOverEEE").value(),filt.getParameter("thrOverE2EE").value())
    print "{}{} {} (EB) {} (EE)".format(" "*indent,var_name,cut_barrel.label(),cut_endcap.label())

    

def process_filter_sel(filt,indent=6):
    filt_type = filt.type_()
    if filt_type in ["HLTBool","HLTPrescaler","HLTTriggerTypeFilter","HLTL1TSeed"]: return 
    try:
        globals()["print"+filt.type_()](filt,indent=6)
    except KeyError:
        print "{}filter type {} does not have a defined print function print{}, please add it".format(" "*indent,filt.type_(),filt.type_())

def print_filter_sel(process,path_name):
    print "  {}:".format(path_name)
    path = getattr(process,path_name)
    for filter_name in str(path).split("+"):
        if filter_name.find("ignore(")==0:
            print "{}skipping {} as can not yet handle ignored filters".format(" "*6,filter_name)
            continue
        filt = getattr(process,filter_name)
        if type(filt).__name__=="EDFilter":
            process_filter_sel(filt,indent=6)


        

def filter_module_names(process,path_name):
    filter_names = []
    path = getattr(process,path_name)
    for filter_name in path.moduleNames():
        filt = getattr(process,filter_name)
        if type(filt).__name__=="EDFilter":
            filter_names.append(filter_name)
    return filter_names
    


import argparse
import importlib

def main():

    parser = argparse.ArgumentParser(description='dumps the save tag filters of a menu')

    parser.add_argument('hlt_menu_name',help="the python file containing the hlt menu")
    args = parser.parse_args()

    with open(args.hlt_menu_name) as f:
        exec f.read()

#   print process
#    mod = importlib.import_module(args.hlt_menu_name)
#    process = getattr(mod,"process")

    for path_name in process.pathNames().split():
        if path_name.find("HLT_")==0 and is_egamma_path(process,path_name):
            print_filter_sel(process,path_name)

# print_filter_sel(process,"HLT_Ele27_WPTight_Gsf_v4")
if __name__ == "__main__":
    main()

