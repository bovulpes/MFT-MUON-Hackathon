#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <TFile.h>
#include <TTree.h>

#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/Digit.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTrack.h"

#endif

void Read_digROF_3() {

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

  // DataFormats/Detectors/ITSMFT/common/include/DataFormatsITSMFT/ROFRecord.h
  using ROFRec = o2::itsmft::ROFRecord;
  using MC2ROF = o2::itsmft::MC2ROFRecord;
  // DataFormats/Detectors/ITSMFT/common/include/DataFormatsITSMFT/Digit.h
  using Digit = o2::itsmft::Digit;

  TFile fileD("mftdigits.root");
  TTree *digTree = (TTree*)fileD.Get("o2sim");
  std::vector<MC2ROF> mc2rofVec, *mc2rofVecP = &mc2rofVec;
  std::vector<o2::itsmft::Digit> digVec, *digVecP = &digVec;
  digTree->SetBranchAddress("MFTDigit", &digVecP);
  if (digTree->GetBranch("MFTDigitMC2ROF")) {
    digTree->SetBranchAddress("MFTDigitMC2ROF", &mc2rofVecP);
  } else {
    printf("No Monte-Carlo information in this file\n");
    return;
  }
  
  std::vector<ROFRec> rofRecVec, *rofRecVecP = &rofRecVec;
  digTree->SetBranchAddress("MFTDigitROF", &rofRecVecP);

  int nEntries = digTree->GetEntries();
  printf("Number of entries in digits tree %d \n", nEntries);
  
  digTree->GetEntry(0);

  int nROFRec = (int)rofRecVec.size();
  printf("Found %d ROF records\n", nROFRec);
  std::vector<int> mcEvMin(nROFRec, nEvents);
  std::vector<int> mcEvMax(nROFRec, -1);

  for (int imc = mc2rofVec.size(); imc--;) {
    const auto& mc2rof = mc2rofVec[imc];
    //mc2rof.print();
    printf("ROFRecord ID %3d   min %3d   max %3d   from event %2d \n",
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

  fileD.Close();
  
}
