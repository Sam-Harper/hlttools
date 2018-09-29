#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms
import math
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
        "hltEgammaHToverET":"H<sub>T</sub> (for H/E, cone<0.14)",
        "hltEgammaGsfTrackVars:DetaSeed":"DEtaInSeed",
        "hltEgammaGsfTrackVars:Dphi":"DPhiIn",
        "hltEgammaGsfTrackVars:MissingHits":"#miss hits",
        "hltEgammaGsfTrackVars:Chi2":"trk chi2",
        "hltEgammaGsfTrackVars:OneOESuperMinusOneOP":"1/E-1/p",
        "hltEgammaEcalPFClusterIso":"ECAL Clus Iso (rho corr,cone<0.3)",
        "hltEgammaEcalPFClusterIsoR02" : "ECAL Clus Iso (rho corr,cone<0.2)",
        "hltEgammaHcalPFClusterIso" : "HCAL Clus Iso (rho corr,cone<0.3)",
        "hltEgammaHcalPFClusterIsoR02" : "HCAL Clus Iso (rho corr,cone<0.2)",
        "hltEgammaEleGsfTrackIso":"Ele Trk Iso (cone<0.3)",
        "hltEgammaEleGsfTrackIsoR02" : "Ele Trk Iso (cone<0.2)",
        "hltEgammaClusterShape:sigmaIEtaIEta5x5":"sigmaIEtaIEta (full5x5)",
        "hltEgammaClusterShape":"sigmaIEtaIEta (fractions)",
        "hltEgammaR9ID":"R9 (fractions)",
        "hltEgammaR9ID:r95x5":"R9 (full5x5)",
        "hltEgammaHollowTrackIso":"Pho Trk Iso (hollow cone, cone<0.29)",
        "hltEgammaHoverERhoCorr" : "H (for H/E, cone<0.14, rho corr)",
        
        

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
            cut = "pass" if filt.pixelVeto.value() else "veto" 
            if filt.useS.value():
                s2_cal = lambda x : "{:.1f}".format((math.atanh(x)*10)**2) if x<1.0 else "inf" if x==1.0 else "-inf"  
                cut+=" with old s2 < {} BPIX, < {} BPIX-FPIX, < {} FPIX".format( s2_cal(filt.tanhSO10BarrelThres.value()),s2_cal(filt.tanhSO10InterThres.value()),s2_cal(filt.tanhSO10ForwardThres.value()))
            self.cuts = [cut]
            
        else:
        #    print filt,filt.type_()
            if "varTag" in filt.parameterNames_():
                self.var = get_nice_var_name(filt.getParameter("varTag").value().replace("Unseeded","")) 
            elif "isoTag" in filt.parameterNames_():
                self.var = get_nice_var_name(filt.getParameter("isoTag").value().replace("Unseeded","")) 
            else:
                raise RuntimeError(str(filt)+" "+filt.type_()+" has no known way to get its varible")

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
                if "l1EGCand" in filt.parameterNames_():
                    self.l1_seeded = filt.getParameter("l1EGCand").value().find("Unseeded")==-1
                else:
                    self.l1_seeded = filt.getParameter("L1IsoCand").value().find("Unseeded")==-1
            elif filt.type_() == "HLT1Photon":
                self.min_et_eb = max(self.min_et_eb,filt.getParameter("MinPt").value())
                self.min_et_ee = max(self.min_et_ee,filt.getParameter("MinPt").value())
                self.ncands = max(self.ncands,filt.getParameter("MinN").value())
                self.eta_bins[-1] = min(self.eta_bins[-1],filt.getParameter("MaxEta").value())
                self.l1_seeded = filt.getParameter("inputTag").value().find("Unseeded")==-1
            else:
                self.cuts.append(EGammaCut(filt))
                self.ncands = max(self.ncands,filt.getParameter("ncandcut").value())

    def __str__(self):
        out_str = "E_{{T}} > {} GeV (EB), > {} GeV (EE), #cands = {}, L1 seeded = {}\n".format(self.min_et_eb,self.min_et_ee,self.ncands,self.l1_seeded)
        out_str += "|  *var*  |"
        for binnr in range(0,len(self.eta_bins)-1):
            out_str += "  *{} < &#124;&eta;&#124; < {}*  |".format(self.eta_bins[binnr],self.eta_bins[binnr+1])
        out_str += "\n"
        for cut in self.cuts:
            out_str += "|  {}  | ".format(cut.var)
            for binnr in range(0,len(self.eta_bins)-1):
                out_str += "  {}  |".format(cut.get_cut(self.eta_bins[binnr]))
            out_str+="\n"
        return out_str
    
def rm_filter_modifiers(filt_name):
    filt_start = filt_name.find("(")
    filt_end = filt_name.rfind(")")
    if filt_start!=-1: return filt_name[filt_start+1:filt_end]
    else: return str(filt_name)

def is_valid_egid_filt_type(filt): 
    if type(filt).__name__=="EDFilter":
        #so we have a black list rather than a white list so we dont miss new E/gamma ID modules
        if filt.type_() in ['HLTTriggerTypeFilter','HLTBool','HLTPrescaler','HLTTriggerTypeFilter','HLTL1TSeed',"CaloJetSelector","CandViewCountFilter","CandViewSelector","EtMinCaloJetSelector","EtaRangeCaloJetSelector","HLT1CaloJet","HLT1CaloMET","HLT1PFJet","HLT1PFMET","HLT1PFTau","HLT2PFJetPFJet","HLT2PhotonMET","HLT2PhotonPFMET","HLT2PhotonPFTau","HLT2PhotonPhotonDZ","HLT2PhotonTau","HLTCaloJetTag","HLTCaloJetVBFFilter","HLTEgammaAllCombMassFilter","HLTEgammaCombMassFilter","HLTEgammaDoubleLegCombFilter","HLTElectronMuonInvMassFilter","HLTHtMhtFilter","HLTMhtFilter","HLTMuonIsoFilter","HLTMuonL1TFilter","HLTMuonL2FromL1TPreFilter","HLTMuonL3PreFilter","HLTPFJetCollectionsFilter","HLTPFJetTag","HLTPFTauPairDzMatchFilter","HLTPMMassFilter","JetVertexChecker","LargestEtCaloJetSelector","PFTauSelector","PrimaryVertexObjectFilter","VertexSelector",'HLTEgammaL1TMatchFilterRegional','HLTEgammaTriggerFilterObjectWrapper',"HLT2PhotonMuonDZ","HLT2MuonPhotonDZ","MuonSelector","HLTMuonDimuonL3Filter","HLTDisplacedmumuFilter"]: return False
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




def get_path_sel(process,path_name):

    sel_str = ""
    path = getattr(process,path_name)    
    
    filter_chains = []
    chains_filter_names = split_into_chains(process,str(path).split("+"))
    for chain in chains_filter_names:
        cutcoll = EGammaCutColl()
        cutcoll.fill(l1_seeded=True,process=process,filter_names=chain)
           
        sel_str += str(cutcoll) +"\n"
    return sel_str



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
import json

def main():

    parser = argparse.ArgumentParser(description='dumps the save tag filters of a menu')
  #  parser.add_argument('hlt_menu_name',help="the python file containing the hlt menu")
    parser.add_argument('--out',required=True,help='output file to write the json')
    args = parser.parse_args()


#   print process
#    mod = importlib.import_module(args.hlt_menu_name)
#    process = getattr(mod,"process")

    menu_versions = ["2016_v1p1","2016_v1p2","2016_v2p1","2016_v2p2","2016_v3p0","2016_v3p1","2016_v4p1","2016_v4p2"]
 #   menu_versions = ["2018_test"]
    hlt_sel = {}

    for version in menu_versions:
        hlt_menu = "hltMenus/hltMenu_{}.py".format(version)
  #      hlt_menu = "testMenu2018.py"
        with open(hlt_menu) as f:
            exec f.read()
        for path_name in process.pathNames().split():
            if path_name.find("HLT_")==0 and is_egamma_path(process,path_name):
                path_sel = get_path_sel(process,path_name)
                path_name_woversion =  rm_hlt_version(path_name)
                if path_name_woversion not in hlt_sel:
                    hlt_sel[path_name_woversion] = {}
                hlt_sel[path_name_woversion][version] = path_sel
        del process

    print "warning, this is still underconstruction<br> Known issues: ORed paths not listed as ORed, L1 seeding not identified"
    print "%TOC%"
    
    path_names = hlt_sel.keys()
    path_names.sort()
    for path_name in path_names:
        print "---++ "+path_name
        ver_str = ""
        pre_sel = ""
        for menu_version in menu_versions:
            try:
                if hlt_sel[path_name][menu_version]==pre_sel:
                    ver_str+=", "+menu_version
                else:
                    if ver_str != "":
                        print "menu versions : {} <br>".format(ver_str)
                        print pre_sel
                    pre_sel = hlt_sel[path_name][menu_version]
                    ver_str = str(menu_version)
            except KeyError:
                pass
        print "menu versions : {} <br>".format(ver_str)
        print pre_sel
    with open(args.out,'w') as f:
        json.dump(hlt_sel,f)
            
if __name__ == "__main__":
    main()

