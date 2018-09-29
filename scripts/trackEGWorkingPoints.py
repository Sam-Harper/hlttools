#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

def is_egamma_path(process,path_name):
    """A path is an E/gamma path if it contains hltEgammaCandidates or hltEgammaCandidatesUnseeded in the modules it runs"""   
    path = getattr(process,path_name)
    return "hltEgammaCandidates" in path.moduleNames() or "hltEgammaCandidatesUnseeded" in path.moduleNames()

def rm_hlt_version(name):
    version_start = name.rfind("_v")
    if version_start == -1: 
        return name
    else:
        return name[:version_start+2]    

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
    def __init__(self,op_type="",divide_by_var="E",const_term=-1.0,linear_term=-1.0,quad_term=-1.0,term_combine_op="||",min_eta = 0, max_eta = 2.65):
        self.op_type = op_type
        self.divide_by_var = divide_by_var
        self.const_term = const_term
        self.linear_term = linear_term
        self.quad_term = quad_term
        self.term_combine_op = term_combine_op
        self.min_eta = min_eta
        self.max_eta = max_eta
    
    def valid_for_eta(self,eta):
        return eta >= self.min_eta and eta< self.max_eta

    def __str__(self):
        label_str = "{}".format(self.op_type)
        first_term = True
        ignore_value = -1.0 if self.term_combine_op=="||" else 0.
        if self.const_term!=ignore_value:
            label_str+=" {}".format(self.const_term)
            first_term = False
        if self.linear_term!=ignore_value:
            if not first_term: label_str+=" {} ".format(self.term_combine_op)
            label_str+=" {} * {}".format(self.linear_term,self.divide_by_var)
            first_term = False
        if self.quad_term!=ignore_value:
            if not first_term: label_str+=" {} ".format(self.term_combine_op)
            label_str+=" {} * {}^{{2}}".format(self.quad_term,self.divide_by_var)
            first_term = False

        #we need to deal with the special case where its all zeros
        if first_term and ignore_value==0.:
            label_str+=" 0."

        return label_str

    def label(self):
        return str(self)

class EGammaCut:
    def __init__(self,filt=None):
        if filt.type_()=="HLTElectronPixelMatchFilter":
            self.var = "pixel match"
            self.cuts = ["pass tmp"]
        else:
            self.var = get_nice_var_name(filt.getParameter("varTag").value()) 
            self.cuts = []

            op_str = "<=" if filt.getParameter("lessThan").value() else ">="
            et_str = "E_{T}" if filt.getParameter("useEt").value() else "E" 
            term_op = "||" if  filt.type_() == "HLTEgammaGenericFilter" else "+"
            self.cuts.append(EGammaStdCut(op_str,et_str,filt.getParameter("thrRegularEB").value(),filt.getParameter("thrOverEEB").value(),filt.getParameter("thrOverE2EB").value(),term_op,min_eta=0,max_eta=1.479))
            self.cuts.append(EGammaStdCut(op_str,et_str,filt.getParameter("thrRegularEE").value(),filt.getParameter("thrOverEEE").value(),filt.getParameter("thrOverE2EE").value(),term_op,min_eta=0,max_eta=2.65))

    def get_cut(self,eta):
        for cut in self.cuts:
            try:
                if cut.valid_for_eta(eta): return cut
            except AttributeError:
                return cut
        return None

class EGammaCutColl:
    
    def __init__(self):
        self.eta_bins = [0,1.479,2.65]
        self.cuts = []
        self.ncands = 0
        self.l1_seeded = False
        self.min_et_eb = 0
        self.min_et_ee = 0

    def fill(self,process,filter_names,l1_seeded=True):
        for filter_name in filter_names:
            filt = getattr(process,filter_name)
            if filt.type_() == "HLTEgammaEtFilter":
                self.min_et_eb = max(self.min_et_eb,filt.getParameter("etcutEB").value())
                self.min_et_ee = max(self.min_et_ee,filt.getParameter("etcutEE").value())
                self.ncands = max(self.ncands,filt.getParameter("ncandcut").value())
            elif filt.type_() == "HLT1Photon":
                self.min_et_eb = max(self.min_et_eb,filt.getParameter("MinPt").value())
                self.min_et_ee = max(self.min_et_ee,filt.getParameter("MinPt").value())
                self.ncands = max(self.ncands,filt.getParameter("MinN").value())
                self.eta_bins[-1] = min(self.eta_bins[-1],filt.getParameter("MaxEta").value())
            else:
                self.cuts.append(EGammaCut(filt))
                self.ncands = max(self.ncands,filt.getParameter("ncandcut").value())

    def __str__(self):
        out_str = "E_{{T}} > {} GeV (EB), {} GeV (EE), #cands = {}, L1 seeded = {}\n".format(self.min_et_eb,self.min_et_ee,self.ncands,self.l1_seeded)
        out_str += "| *var* |"
        for binnr in range(0,len(self.eta_bins)-1):
            out_str += " *{} < &#124;&eta;&#124; < {}* |".format(self.eta_bins[binnr],self.eta_bins[binnr+1])
        out_str += "\n"
        for cut in self.cuts:
            out_str += "| {} | ".format(cut.var)
            for binnr in range(0,len(self.eta_bins)-1):
                out_str += " {} |".format(cut.get_cut(self.eta_bins[binnr]))
            out_str+="\n"
        return out_str
    
class EGammaFilter:
    def __init__(self,filt_mod=None):
        if(filt_mod): self.fill(filt_mod)
        
    def fill(self,filt_mod):
        self.filt_type = filt_mod.type_()
        self.cut_params = {}
        self.var = None
        if self.filt_type in ["HLTEgammaGenericFilter","HLTEgammaGenericQuadraticFilter"]:
            self.cut_params = filt_mod.parameters_()
            self.var = filt_mod.varTag.value().replace("Unseeded","")
            #now remove the parameters which are not part of the selection from cut_params
            for x in ['saveTags','l1EGCand','candTag','ncandcut']: del self.cut_params[x]
            #now clean large values (999 is the same as 99999 etc)
            for param in self.cut_params:
                if type(self.cut_params[param]).__name__=="double" and self.cut_params[param].value()>990:
                    self.cut_params[param].setValue(999.)

        elif self.filt_type=="HLTElectronPixelMatchFilter":
            cut_params = filt_mod.parameters_()
            self.cut_params = {"useS" : cut_params["useS"],"pixelVeto" : cut_params["pixelVeto"] }
            self.var = "pixelMatch"

        else:
            print "debug",self.filt_type,"not recognised"
       
    def __eq__(self, other):
        return self.filt_type == other.filt_type and self.cut_params == other.cut_params
        if self.filt_type == other.filt_type:
            same = True
            for cut_param_name in self.cut_params:
                if self.cut_params[cut_param_name] != other.cut_params[cut_param_name]:
                    print "miss match",self.var,cut_param_name," self ",self.cut_params[cut_param_name],"other",other.cut_params[cut_param_name]
                    same = False
            return same
        return False
                
            
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def __str__(self):
        if self.filt_type=="HLTEgammaGenericFilter":
            op_str = "<=" if self.cut_params["lessThan"].value() else ">="
            et_str = "E_{T}" if self.cut_params["useEt"].value() else "E"
            var_name = get_nice_var_name(self.var)
            cut_barrel = EGammaStdCut(op_str,et_str,self.cut_params["thrRegularEB"].value(),self.cut_params["thrOverEEB"].value(),self.cut_params["thrOverE2EB"].value())
            cut_endcap = EGammaStdCut(op_str,et_str,self.cut_params["thrRegularEE"].value(),self.cut_params["thrOverEEE"].value(),self.cut_params["thrOverE2EE"].value())
            return  "{} {} (EB) {} (EE)".format(var_name,cut_barrel.label(),cut_endcap.label())
        else:
            return ""
class EGammaFilterChain:
    def __init__(self):
        self.filters = []
        self.paths_using = []
    def add_filter(self,filt_mod):
        if filt_mod.type_() in ["HLTEgammaGenericFilter",'HLTEgammaGenericQuadraticFilter','HLTElectronPixelMatchFilter']:
            self.filters.append(EGammaFilter(filt_mod))
        elif filt_mod.type_() not in ['HLT1Photon','HLTEgammaEtFilter']:
            print "debug",filt_mod.type_(),"not recognised"

    def add_path(self,path):
        self.paths_using.append(path)

    def __eq__(self,other):
        return self.filters == other.filters

    def __ne__(self,other):
        return not self.__eq__(other)

    def __getitem__(self,index):
        return self.filters[index]

    def __str__(self):
        print_str = "start chain "+str(self.paths_using)
        for filt in self.filters:
            print_str+="\n  "+str(filt)
        return print_str

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

def rm_filter_modifiers(filt_name):
    filt_start = filt_name.find("(")
    filt_end = filt_name.rfind(")")
    if filt_start!=-1: return filt_name[filt_start+1:filt_end]
    else: return str(filt_name)

def is_valid_egid_filt_type(filt): 
    if type(filt).__name__=="EDFilter":
        if filt.type_() in ['HLTTriggerTypeFilter','HLTBool','HLTPrescaler','HLTTriggerTypeFilter','HLTL1TSeed',"CaloJetSelector","CandViewCountFilter","CandViewSelector","EtMinCaloJetSelector","EtaRangeCaloJetSelector","HLT1CaloJet","HLT1CaloMET","HLT1PFJet","HLT1PFMET","HLT1PFTau","HLT2PFJetPFJet","HLT2PhotonMET","HLT2PhotonPFMET","HLT2PhotonPFTau","HLT2PhotonPhotonDZ","HLT2PhotonTau","HLTCaloJetTag","HLTCaloJetVBFFilter","HLTEgammaAllCombMassFilter","HLTEgammaCombMassFilter","HLTEgammaDoubleLegCombFilter","HLTElectronMuonInvMassFilter","HLTHtMhtFilter","HLTMhtFilter","HLTMuonIsoFilter","HLTMuonL1TFilter","HLTMuonL2FromL1TPreFilter","HLTMuonL3PreFilter","HLTPFJetCollectionsFilter","HLTPFJetTag","HLTPFTauPairDzMatchFilter","HLTPMMassFilter","JetVertexChecker","LargestEtCaloJetSelector","PFTauSelector","PrimaryVertexObjectFilter","VertexSelector",'HLTEgammaL1TMatchFilterRegional','HLTEgammaTriggerFilterObjectWrapper']: return False
        else: return True
    else: return False

def get_prev_filt_name(filt):
    if filt.type_() in ['HLTEgammaGenericFilter','HLTEgammaGenericQuadraticFilter','HLTElectronPixelMatchFilter']:
        return filt.candTag.value()
    elif filt.type_() in ['HLTEgammaL1TMatchFilterRegional','HLT1Photon','HLTEgammaTriggerFilterObjectWrapper']:
        return None
    elif filt.type_() in ['HLTEgammaEtFilter']:
        return None#filt.inputTag.value()
    else:
        print filt.type_()+" not recognised"
        return None
        raise RuntimeError(filt.type_()+" not recognised")

    

def split_into_chains(process,filt_names):
    chains = []
    for filt_name in filt_names:
        filt_name = rm_filter_modifiers(filt_name)
        filt = getattr(process,filt_name)
        if not is_valid_egid_filt_type(filt): continue
        prev_filt_name = get_prev_filt_name(filt)
 #       print filt_name,prev_filt_name
        if prev_filt_name==None:
            chains.append([filt_name])
        found=False
        for chain in chains:
            if prev_filt_name in chain:
                found=True
                chain.append(filt_name)
    return chains

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


def print_trig_sel(process,path_name):

    path = getattr(process,path_name)    
    
    filter_chains = []
    chains_filter_names = split_into_chains(process,str(path).split("+"))

    print "---++ "+rm_hlt_version(path_name)
    for chain in chains_filter_names:
        cutcoll = EGammaCutColl()
        cutcoll.fill(l1_seeded=True,process=process,filter_names=chain)
  
        print cutcoll

def make_wpdefs(process,wpoints,path_name):

    path = getattr(process,path_name)    
    
    filter_chains = []
    chains_filter_names = split_into_chains(process,str(path).split("+"))
#    print "path : ",path_name,chains_filter_names
    for chain in chains_filter_names:
        filter_chain = EGammaFilterChain()
        for filter_name in chain:
            filt = getattr(process,filter_name)
            filter_chain.add_filter(filt)
        filter_chains.append(filter_chain)
    #print filter_chain
    #print wpoints
#    wpoints = [filter_chain]
    if filter_chain not in wpoints:
        filter_chain.add_path(path_name)
        wpoints.append(filter_chain)
    else:
        for wp in wpoints:
            if filter_chain == wp:
                wp.add_path(path_name)
                break

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
    wpoints = []
    for path_name in process.pathNames().split():
        if path_name.find("HLT_")==0 and is_egamma_path(process,path_name):
#            if path_name.find("WPLoose")==-1: continue
         #   print_filter_sel(process,path_name)
            print_trig_sel(process,path_name)
#            make_wpdefs(process,wpoints,path_name)
#    print len(wpoints)
#    for wp in wpoints:
#        print wp
# print_filter_sel(process,"HLT_Ele27_WPTight_Gsf_v4")
if __name__ == "__main__":
    main()

