#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <TFile.h>
#include <TTree.h>

#include "DataFormatsITSMFT/ROFRecord.h"
#include "SimulationDataFormat/MCTrack.h"

#endif

void Read_clsROF_3() {

  // distribution of MC events over the Read-Out-Frames

  // get the number of MC events
  using o2::MCTrack;
  TFile fileK("o2sim_Kine.root");
  TTree* kineTree = (TTree*)fileK.Get("o2sim");
  std::vector<o2::MCTrack> mcTrkVec, *mcTrkVecP = &mcTrkVec;
  kineTree->SetBranchAddress("MCTrack",&mcTrkVecP);
  int nEvents = kineTree->GetEntries();
  printf("Number of events %d \n", nEvents);
  fileK.Close();

  using ROFRec = o2::itsmft::ROFRecord;
  using MC2ROF = o2::itsmft::MC2ROFRecord;

  TFile fileC("mftclusters.root");
  TTree *clsTree = (TTree*)fileC.Get("o2sim");
  std::vector<o2::itsmft::MC2ROFRecord> mc2rofVec, *mc2rofVecP = &mc2rofVec;
  if (clsTree->GetBranch("MFTClustersMC2ROF")) {
    clsTree->SetBranchAddress("MFTClustersMC2ROF", &mc2rofVecP);
  } else {
    printf("No Monte-Carlo information in this file\n");
    return;
  }

  std::vector<ROFRec> rofRecVec, *rofRecVecP = &rofRecVec;
  clsTree->SetBranchAddress("MFTClustersROF", &rofRecVecP);

  int nEntries = clsTree->GetEntries();
  printf("Number of entries in clusters tree %d \n", nEntries);
  
  clsTree->GetEntry(0);

  int nROFRec = (int)rofRecVec.size();
  printf("Found %d ROF records\n", nROFRec);
  std::vector<int> mcEvMin(nROFRec, nEvents);
  std::vector<int> mcEvMax(nROFRec, -1);

  for (int imc = mc2rofVec.size(); imc--;) {
    const auto& mc2rof = mc2rofVec[imc];
    //mc2rof.print();
    printf("ROFRecord ID %d min %d max %d from event %d \n",
	   mc2rof.rofRecordID,
	   mc2rof.minROF, mc2rof.maxROF, mc2rof.eventRecordID);
    for (int irofd = (mc2rof.maxROF - mc2rof.minROF + 1); irofd--;) {
      int irof = mc2rof.rofRecordID + irofd;
      if (mcEvMin[irof] > imc) {
        mcEvMin[irof] = imc;
      }
      if (mcEvMax[irof] < imc) {
        mcEvMax[irof] = imc;
      }
    }
  }

  for (int irof = 0; irof < nROFRec; irof++) {
    if (mcEvMin[irof] > mcEvMax[irof]) {
      printf("In ROF %3d found no MC events!\n", irof); 
    } else {
      printf("In ROF %3d found MC events %2d to %2d \n", irof, mcEvMin[irof], mcEvMax[irof]);
    }
  }

  fileC.Close();
  
}
