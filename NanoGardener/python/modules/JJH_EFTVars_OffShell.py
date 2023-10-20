import ROOT
import math 
import ctypes
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.modules.common.collectionMerger import collectionMerger
from LatinoAnalysis.NanoGardener.framework.BranchMapping import mappedOutputTree, mappedEvent

import os.path

class JJH_EFTVars_OffShell(Module):
    def __init__(self, branch_map=''):
        
        self.cmssw_base = os.getenv('CMSSW_BASE')
        self.cmssw_arch = os.getenv('SCRAM_ARCH')

        ROOT.gSystem.AddIncludePath("-I"+self.cmssw_base+"/interface/")
        ROOT.gSystem.AddIncludePath("-I"+self.cmssw_base+"/src/")
        #ROOT.gSystem.AddIncludePath("-I"+self.cmssw_base+"/lib/")
        ROOT.gSystem.Load("libJHUGenMELAMELA.so") 
        ROOT.gSystem.Load(self.cmssw_base+"/src/JHUGenMELA/MELA/data/"+self.cmssw_arch+"/libmcfm_707.so")

        try:
            ROOT.gROOT.LoadMacro(self.cmssw_base+'/src/LatinoAnalysis/Gardener/python/variables/melaHiggsEFT.C+g')
        except RuntimeError:
            ROOT.gROOT.LoadMacro(self.cmssw_base+'/src/LatinoAnalysis/Gardener/python/variables/melaHiggsEFT.C++g')
      
        self.mela = ROOT.Mela(13, 125,  ROOT.TVar.SILENT) 

        self._branch_map = branch_map

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = mappedOutputTree(wrappedOutputTree, mapname=self._branch_map)
        self.newbranches = [
          'hm','me_qqh_hsm','me_ggh_hsm', 'D2jVBF'
          ]
        
        for nameBranches in self.newbranches :
          self.out.branch(nameBranches , "F");

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        event = mappedEvent(event, mapname=self._branch_map)

        Lepton = Collection(event, "Lepton")
        nLepton = len(Lepton)

        Jet   = Collection(event, "CleanJet")
        nJet = len(Jet)

        OrigJet = Collection(event, "Jet")

        hm = -999
        me_qqh_hsm = -999 
        me_ggh_hsm = -999 
        D2jVBF = -999     # newly added 
       

        if nJet > 1 and nLepton > 1 :
        #if nJet == 2 and nLepton == 2 :
         
         L1 = ROOT.TLorentzVector()
         L2 = ROOT.TLorentzVector()
         L1.SetPtEtaPhiM(Lepton[0].pt, Lepton[0].eta, Lepton[0].phi, 0)
         L2.SetPtEtaPhiM(Lepton[1].pt, Lepton[1].eta, Lepton[1].phi, 0)

         LL = ROOT.TLorentzVector()
         LL = L1 + L2

         MET_phi   = event.PuppiMET_phi
         MET_pt    = event.PuppiMET_pt

         NuNu = ROOT.TLorentzVector()
         nunu_px = MET_pt*math.cos(MET_phi)
         nunu_py = MET_pt*math.sin(MET_phi)
         nunu_pz = LL.Pz()
         nunu_m  = LL.M()     # change this value a/c to gen studies done for reco matrix element
         nunu_e  = math.sqrt(nunu_px*nunu_px + nunu_py*nunu_py + nunu_pz*nunu_pz + nunu_m*nunu_m)
         NuNu.SetPxPyPzE(nunu_px, nunu_py, nunu_pz, nunu_e)

         Higgs = ROOT.TLorentzVector()
         Higgs = LL + NuNu
         hm  = Higgs.M()

         print ('hm', hm)

         indx_j1 = 0
         indx_j2 = 1 

         J1 = ROOT.TLorentzVector()
         J2 = ROOT.TLorentzVector() 
         indx_oj1 = Jet[indx_j1].jetIdx
         indx_oj2 = Jet[indx_j2].jetIdx
         J1.SetPtEtaPhiM(Jet[indx_j1].pt, Jet[indx_j1].eta, Jet[indx_j1].phi, OrigJet[indx_oj1].mass)
         J2.SetPtEtaPhiM(Jet[indx_j2].pt, Jet[indx_j2].eta, Jet[indx_j2].phi, OrigJet[indx_oj2].mass)

         daughter_coll = ROOT.SimpleParticleCollection_t() 
         associated_coll = ROOT.SimpleParticleCollection_t()

         daughter = ROOT.SimpleParticle_t(25, Higgs)
         associated1 = ROOT.SimpleParticle_t(0, J1)
         associated2 = ROOT.SimpleParticle_t(0, J2)

         daughter_coll.push_back(daughter)                                                         
         associated_coll.push_back(associated1)
         associated_coll.push_back(associated2)

         self.mela.setCandidateDecayMode(ROOT.TVar.CandidateDecay_Stable)   
         self.mela.setInputEvent(daughter_coll, associated_coll, 0, 0)
         self.mela.setCurrentCandidateFromIndex(0)

         ME_VBF = ROOT.melaHiggsEFT(self.mela, ROOT.TVar.JHUGen, ROOT.TVar.JJVBF, 0, 1) 
        
         #print(type(ME_VBF)) 
         #print ('ME_VBF[0]',ME_VBF[0])
         #print ('ME_VBF[1]',ME_VBF[1])
         #print ('check ME_VBF',ME_VBF)

         me_qqh_hsm   = ME_VBF.find("me_hsm").second
         #print ('me_qqh_hsm', me_qqh_hsm)
        
         ME_QCD = ROOT.melaHiggsEFT(self.mela, ROOT.TVar.JHUGen, ROOT.TVar.JJQCD, 1, 1)
         #me_qcd_hsm = ROOT.melaHiggsEFT_OffShell(self.mela, ROOT.TVar.JHUGen, ROOT.TVar.JJQCD, 1, 1)
         me_ggh_hsm   = ME_QCD.find("me_hsm").second
         #print ('me_ggh_hsm', me_ggh_hsm)

         #D2jVBF = 1/(1+(me_ggh_hsm)/(me_qqh_hsm))
         # Calculate D2jVBF with a default value of 0.0 when the denominator is zero

         #if me_qqh_hsm + me_ggh_hsm != 0:
         #   D2jVBF = me_qqh_hsm / (me_qqh_hsm + me_ggh_hsm)
         #else:
         #   D2jVBF = 0.0
         if me_qqh_hsm != 0:
            D2jVBF = 1/(1+((me_ggh_hsm)/(me_qqh_hsm)))
         else:
            D2jVBF = -200.0
             
             
         #print ('D2jVBF', D2jVBF)

         self.mela.resetInputEvent()
         
         

        self.out.fillBranch( 'hm',  hm )
        self.out.fillBranch( 'me_qqh_hsm',  me_qqh_hsm )
        self.out.fillBranch( 'me_ggh_hsm',  me_ggh_hsm )
        self.out.fillBranch('D2jVBF',  D2jVBF)
        

        return True











