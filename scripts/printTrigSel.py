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
        "hltEgammaSolidTrackIso": "Pho Trk Iso (solid cone, cone<0.29)",
        "hltEgammaHoverERhoCorr" : "H (for H/E, cone<0.14, rho corr)",
        "hltEgammaPixelMatchVars:s2" : "PM S2",
        
        

    }.get(name, str(name))  


class EGammaStdCut:
    def __init__(self,op_type="",divide_by_var="E",const_term=-1.0,linear_term=-1.0,quad_term=-1.0,rho_term=0,term_combine_op="||",min_eta = 0, max_eta = 2.65):
        self.op_type = op_type
        self.divide_by_var = divide_by_var
        self.const_term = const_term
        self.linear_term = linear_term
        self.quad_term = quad_term
        self.rho_term = rho_term
        self.term_combine_op = term_combine_op
        self.min_eta = min_eta
        self.max_eta = max_eta
    
    def valid_for_eta(self,eta):
        return eta >= self.min_eta and eta< self.max_eta

    def __str__(self):
        label_str = "{}".format(self.op_type)
        first_term = True
        ignore_value = -1.0 if self.term_combine_op=="||" else 0.
        #do we need to rho_corr the individual terms (ie when ORing) or can we just at it at the end
        rho_corr_terms = self.term_combine_op=="||" and self.rho_term != 0.
        if self.const_term!=ignore_value:
            label_str+=" {}".format(self.const_term)
            if rho_corr_terms:
                label_str+=" {}*rho".format(self.rho_term)
            first_term = False
        if self.linear_term!=ignore_value:
            if not first_term: label_str+=" {} ".format(self.term_combine_op)
            if rho_corr_terms:
                label_str+=" {} * {}".format(self.linear_term,self.divide_by_var)
            else:
                label_str+=" ({} + {}*rho) * {}".format(self.linear_term,self.rho_term,self.divide_by_var)
            first_term = False
        if self.quad_term!=ignore_value:
            if not first_term: label_str+=" {} ".format(self.term_combine_op)
            if rho_corr_terms:
                label_str+=" {} * {}^{{2}}".format(self.quad_term,self.divide_by_var)
            else:
                label_str+=" ({} + {}*rho) * {}^{{2}}".format(self.quad_term,self.rho_term,self.divide_by_var)
            first_term = False

        if  not rho_corr_terms and self.rho_term!=0.:
            label_str+=" + {}*rho"
            first_term = False

        #we need to deal with the special case where its all zeros
        if first_term and ignore_value==0.:
            label_str+=" 0."

        return label_str

    def label(self):
        return str(self)

class EGammaCut:
    def __init__(self,filt=None,subchain=0,ignored=False):
        self.subchain = subchain
        self.ignored = ignored
        if filt.type_()=="HLTElectronPixelMatchFilter":
            self.var = "pixel match"
            cut = "pass" if filt.pixelVeto.value()==False else "veto" 
            if filt.useS.value():
                s2_cal = lambda x : "{:.1f}".format((math.atanh(x)*10)**2) if x<1.0 else "inf" if x==1.0 else "-inf"  
                cut+=" with old s2 < {} BPIX, < {} BPIX-FPIX, < {} FPIX".format( s2_cal(filt.tanhSO10BarrelThres.value()),s2_cal(filt.tanhSO10InterThres.value()),s2_cal(filt.tanhSO10ForwardThres.value()))
            self.cuts = [cut]
        elif filt.type_()=="HLTDisplacedEgammaFilter":
            self.var = "displaced ID"
            self.cuts = "pass"
            
        elif filt.type_()=="HLTEgammaGenericQuadraticEtaFilter":
            self.cuts = []  
            self.var = get_nice_var_name(filt.getParameter("varTag").value().replace("Unseeded","")) 
            #right we're going to simplify this and limit ourselfs to certain cases
            #mainly as I think other cases wont occur and so not to waste time coding for them
            if len(filt.energyLowEdges.value())!=1 or filt.energyLowEdges.value()[0]!=0.:
                raise ValueError("can only handle the case of a single energy threshold at zero, for "+str(filt)+" got "+str(filt.energyLowEdges.value()))
            if len(filt.absEtaLowEdges.value())!=4 or filt.etaBoundaryEB12.value()!=filt.absEtaLowEdges.value()[1] or filt.etaBoundaryEE12.value()!=filt.absEtaLowEdges.value()[3] or filt.absEtaLowEdges.value()[0]!=0. or filt.absEtaLowEdges.value()[2]!=1.479:
                raise ValueError("can only handle the case where the rho corr bins are exactly the same as a cut bins\n"+filt.dumpPython())
            

            op_str = "<=" if filt.lessThan.value() else ">="
            et_str = "E_{T}" if filt.useEt.value() else "E" 
            term_op = "+"
            rho_term = 0.
            eta_low_edges = filt.absEtaLowEdges.value()
            self.cuts.append(EGammaStdCut(op_str,et_str,filt.getParameter("thrRegularEB1").value(),filt.getParameter("thrOverEEB1").value(),filt.getParameter("thrOverE2EB1").value(),rho_term,term_op,min_eta=eta_low_edges[0],max_eta=eta_low_edges[1]))
            self.cuts.append(EGammaStdCut(op_str,et_str,filt.getParameter("thrRegularEB2").value(),filt.getParameter("thrOverEEB2").value(),filt.getParameter("thrOverE2EB2").value(),rho_term,term_op,min_eta=eta_low_edges[1],max_eta=eta_low_edges[2]))
            self.cuts.append(EGammaStdCut(op_str,et_str,filt.getParameter("thrRegularEE1").value(),filt.getParameter("thrOverEEE1").value(),filt.getParameter("thrOverE2EE1").value(),rho_term,term_op,min_eta=eta_low_edges[2],max_eta=eta_low_edges[3]))
            self.cuts.append(EGammaStdCut(op_str,et_str,filt.getParameter("thrRegularEE2").value(),filt.getParameter("thrOverEEE2").value(),filt.getParameter("thrOverE2EE2").value(),rho_term,term_op,min_eta=eta_low_edges[3],max_eta=2.65))
            
        
        
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
            rho_term = 0.
            self.cuts.append(EGammaStdCut(op_str,et_str,filt.getParameter("thrRegularEB").value(),filt.getParameter("thrOverEEB").value(),filt.getParameter("thrOverE2EB").value(),rho_term,term_op,min_eta=0,max_eta=1.479))
            self.cuts.append(EGammaStdCut(op_str,et_str,filt.getParameter("thrRegularEE").value(),filt.getParameter("thrOverEEE").value(),filt.getParameter("thrOverE2EE").value(),rho_term,term_op,min_eta=0,max_eta=2.65))
           
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
            filt = getattr(process,filter_name['name'])
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
                self.cuts.append(EGammaCut(filt,subchain=filter_name['subchain'],ignored=filter_name['modifier']=="ignore"
))
                self.ncands = max(self.ncands,filt.getParameter("ncandcut").value())

    def __str__(self):
        out_str = "E_{{T}} > {} GeV (EB), > {} GeV (EE), #cands = {}, L1 seeded = {}\n".format(self.min_et_eb,self.min_et_ee,self.ncands,self.l1_seeded)
        out_str += "|  *var*  |"
        for binnr in range(0,len(self.eta_bins)-1):
            out_str += "  *{} < &#124;&eta;&#124; < {}*  |".format(self.eta_bins[binnr],self.eta_bins[binnr+1])
        out_str += "\n"
        for cut in self.cuts:
            colours = ['%TEAL%','%PURPLE%','%MAROON%']
            colour_nr = cut.subchain % len(colours)
            colour_str = "" if not cut.ignored else colours[colour_nr]
            out_str += "|  {}{}  | ".format(colour_str,cut.var)
            for binnr in range(0,len(self.eta_bins)-1):
                out_str += "  {}{}  |".format(colour_str,cut.get_cut(self.eta_bins[binnr]))
            out_str+="\n"
        return out_str
    
def rm_filter_modifiers(filt_name):
    filt_start = filt_name.find("(")
    filt_end = filt_name.rfind(")")
    if filt_start!=-1: return filt_name[filt_start+1:filt_end],filt_name[0:filt_start]
    else: return str(filt_name),None

def is_valid_egid_filt_type(filt): 
    if type(filt).__name__=="EDFilter":
        #so we have a black list rather than a white list so we dont miss new E/gamma ID modules
        if filt.type_() in ['HLTTriggerTypeFilter','HLTBool','HLTPrescaler','HLTTriggerTypeFilter','HLTL1TSeed',"CaloJetSelector","CandViewCountFilter","CandViewSelector","EtMinCaloJetSelector","EtaRangeCaloJetSelector","HLT1CaloJet","HLT1CaloMET","HLT1PFJet","HLT1PFMET","HLT1PFTau","HLT2PFJetPFJet","HLT2PhotonMET","HLT2PhotonPFMET","HLT2PhotonPFTau","HLT2PhotonPhotonDZ","HLT2PhotonTau","HLTCaloJetTag","HLTCaloJetVBFFilter","HLTEgammaAllCombMassFilter","HLTEgammaCombMassFilter","HLTEgammaDoubleLegCombFilter","HLTElectronMuonInvMassFilter","HLTHtMhtFilter","HLTMhtFilter","HLTMuonIsoFilter","HLTMuonL1TFilter","HLTMuonL2FromL1TPreFilter","HLTMuonL3PreFilter","HLTPFJetCollectionsFilter","HLTPFJetTag","HLTPFTauPairDzMatchFilter","HLTPMMassFilter","JetVertexChecker","LargestEtCaloJetSelector","PFTauSelector","PrimaryVertexObjectFilter","VertexSelector",'HLTEgammaL1TMatchFilterRegional','HLTEgammaTriggerFilterObjectWrapper',"HLT2PhotonMuonDZ","HLT2MuonPhotonDZ","MuonSelector","HLTMuonDimuonL3Filter","HLTDisplacedmumuFilter","HLTMuonTrkL1TFilter","HLT2MuonMuonDZ","HLTPFJetVBFFilter"]: return False
        else: return True
    else: return False

def get_prev_filt_name(filt):
    if filt.type_() in ['HLTEgammaGenericFilter','HLTEgammaGenericQuadraticFilter','HLTElectronPixelMatchFilter','HLTEgammaGenericQuadraticEtaFilter']:
        return filt.candTag.value()
    elif filt.type_() in ['HLTDisplacedEgammaFilter']:
        return filt.inputTag.value();
    elif filt.type_() in ['HLTEgammaL1TMatchFilterRegional','HLT1Photon','HLTEgammaTriggerFilterObjectWrapper']:
        return None
    elif filt.type_() in ['HLTEgammaEtFilter']:
        return None#filt.inputTag.value()
    else:
        print filt.type_()+" not recognised"
        return None
        raise RuntimeError(filt.type_()+" not recognised")

    

def split_into_chains(process,filt_names):
    """
    function takes the list of filters (which are assumed to be in path order) 
    and groups them into distinct chains of filters
    aka  [A -> B -> C -> D],[G -> H -> I] etc
    currently a module if in a chain ifs input module is also in that chain
    so the case  A -> B -> C   and A -> B -> D create as a chain A -> B -> C -> D 
    It does have the concept of "subchain" where if the input filter is immediately before 
    it, it keeps the sub chain number, otherwise it increases the sub chain number
    so in the above A, B C, would be sub chain 0 while D would be sub chain 1
    """
    chains = []
    for filt_name in filt_names:
        filt_name,modifier = rm_filter_modifiers(filt_name)
        if modifier != None and modifier != "ignore":
            raise ValueError('filter "{}" has a modifier "{}" that this code can not handle'.format(filt_name,modifier))
        filt = getattr(process,filt_name)
        if not is_valid_egid_filt_type(filt): continue
        prev_filt_name = get_prev_filt_name(filt)

        filt_data = {'name' : filt_name,'subchain' : 0, 'prev_filt' : prev_filt_name, 'modifier' : modifier}
        if prev_filt_name==None:
            chains.append([filt_data])
        else:
            found=False
            for chain in chains:
                if prev_filt_name == chain[-1]['name']:
                    found=True
                    filt_data['subchain'] = chain[-1]['subchain']
                    chain.append(filt_data)
                elif any( data['name'] == prev_filt_name for data in chain):
                    found=True
                    filt_data['subchain'] = chain[-1]['subchain']+1
                    chain.append(filt_data)
                if found: break
                
            if not found:
                chains.append([filt_data])

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
    menu_versions = ["2017_v4p2"]
 #   menu_versions = ["2018_test"]
    hlt_sel = {}

    for menu_version in menu_versions:
        hlt_menu = "hltMenus/hltMenu_{}.py".format(menu_version)
  #      hlt_menu = "testMenu2018.py"
        with open(hlt_menu) as f:
            exec f.read()
        for path_name in process.pathNames().split():
            if path_name.find("HLT_")==0 and is_egamma_path(process,path_name):
             #   if path_name.find("HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_v")!=0: continue
                path_sel = get_path_sel(process,path_name)
                path_name_no_ver =  rm_hlt_version(path_name)
                if path_name_no_ver not in hlt_sel:
                    hlt_sel[path_name_no_ver] = {}
                    hlt_sel[path_name_no_ver]['selection'] = {}
                hlt_sel[path_name_no_ver]['selection'][menu_version] = path_sel
        del process

    print '''
Welcome to the E/g HLT 2016-2018 Path Information Page

This twiki is automatically generated from HLT menus & wbm info, please do not edit this twiki directly as your edits will be lost next time the auto generation script is run. Please contact E/gamma for any edits you wish to make. <br>

warning, this is still underconstruction<br> 

The information on this twiki is thought to be accurate but again as its auto generated, there may be mistakes in edge cases for paths, please let us know if you find any. Additionally, the lumi numbers are approximate and should only be treated as guidence rather than analysis quality numbers. <br>

Known issues: 
   * this is simply a python file which parses the HLT and tries to print what it things the HLT would do when given this config
   * variable definations are hard coded and do not evolve in time
   * code changes in modules may be not taken into account (although unlikely)
   * does not handle DZ, path leg combination filters
   * does not handle inferor leptons or anything icky and hadronic
   * should work for standard paths but werid complex paths may have edge cases
   * right now does not distingush all cases where there is a discontected filter, that is a filter in the path but the filters after it do not depend on it. So it means you may be requiring an object with Et>50  GeV and an object with H/E<0.15 which is different to an object with Et>50 GeV and H/E <0.15. Think we've got most of them but some remain (and those that remain are likely bugs)
'''
    print "%TOC%"
    
    path_names = hlt_sel.keys()
    path_names.sort()
    for path_name in path_names:
        print "---++ "+path_name
        ver_str = ""
        pre_sel = ""
        for menu_version in menu_versions:
            try:
                if hlt_sel[path_name]['selection'][menu_version]==pre_sel:
                    ver_str+=", "+menu_version
                else:
                    if ver_str != "":
                        print "menu versions : {} <br>".format(ver_str)
                        print pre_sel
                    pre_sel = hlt_sel[path_name]['selection'][menu_version]
                    ver_str = str(menu_version)
            except KeyError:
                pass
        print "menu versions : {} <br>".format(ver_str.replace("p","."))
        print pre_sel
    with open(args.out,'w') as f:
        json.dump(hlt_sel,f)
            
if __name__ == "__main__":
    main()

