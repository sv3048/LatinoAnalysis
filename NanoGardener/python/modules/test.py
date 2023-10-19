import ROOT
import os
import re
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True



# Open the ROOT file
file = ROOT.TFile("/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/Summer20UL18_106x_nAODv9_Full2018v9/AddLHE_MEs/nanoLatino_GluGluHToWWToLNuQQ_M1000__part0.root", "READ")

# Access the TTree containing the event data
tree = file.Get("Events")
if not tree:
    print("Error: Unable to find the TTree in the ROOT file.")
    file.Close()
    exit(1)

# Create a TBranch for PdgId
LHEPart_pdgId = ROOT.std.vector('int')()
tree.SetBranchAddress("LHEPart_pdgId", LHEPart_pdgId)

# Initialize event count
eventCount = 0

# Loop over the events
numEntries = tree.GetEntries()
for i in range(numEntries):
    tree.GetEntry(i)

    # Check if PdgId 21 (jGluonJEt) is in the list
    if 21 in LHEPart_pdgId:
        eventCount += 1

# Close the ROOT file
file.Close()

# Print the number of events with PdgId 21
print("Number of events with PdgId 21 (jGluonJEt)",eventCount)
                                                      