#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <TFile.h>
#include <TTree.h>

#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/Digit.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"

#endif

void Read_digROF_2() {

  // read digits by Read-Out-Frame
  
  // DataFormats/Detectors/ITSMFT/common/include/DataFormatsITSMFT/ROFRecord.h
  using ROFRec = o2::itsmft::ROFRecord;
  using MC2ROF = o2::itsmft::MC2ROFRecord;
  // DataFormats/Detectors/ITSMFT/common/include/DataFormatsITSMFT/Digit.h
  using Digit = o2::itsmft::Digit;

  TFile* fileD = TFile::Open("mftdigits.root");
  TTree *digTree = (TTree*)fileD->Get("o2sim");
  std::vector<o2::itsmft::Digit> digVec, *digVecP = &digVec;
  digTree->SetBranchAddress("MFTDigit", &digVecP);
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* labels = nullptr;
  o2::dataformats::IOMCTruthContainerView* labelROOTbuffer = nullptr;
  o2::dataformats::ConstMCTruthContainer<o2::MCCompLabel> constlabels;
  // for backward compatibility we check what is stored in the file
  auto labelClass = digTree->GetBranch("MFTDigitMCTruth")->GetClassName();
  bool oldlabelformat = false;
  if (digTree->GetBranch("MFTDigitMCTruth")) {
    if (TString(labelClass).Contains("IOMCTruth")) {
      // new format
      digTree->SetBranchAddress("MFTDigitMCTruth", &labelROOTbuffer);
    } else {
      // old format
     digTree->SetBranchAddress("MFTDigitMCTruth", &labels);
     oldlabelformat = true;
    }
  } else {
    printf("No Monte-Carlo information in this file\n");
    return;
  }
  
  std::vector<ROFRec> rofRecVec, *rofRecVecP = &rofRecVec;
  digTree->SetBranchAddress("MFTDigitROF", &rofRecVecP);

  int nEntries = digTree->GetEntries();
  printf("Number of entries in digits tree %d \n", nEntries);
  
  digTree->GetEntry(0);
  if (!oldlabelformat) {
    labelROOTbuffer->copyandflatten(constlabels);
  }
  int nROFRec = (int)rofRecVec.size();
  printf("Found %d ROF records\n", nROFRec);

  int srcID, trkID, evnID;
  bool fake;
  for (int irof = 0; irof < nROFRec; irof++) {
    const auto& rofRec = rofRecVec[irof];
    //rofRec.print();
    printf("ROFRecord ID %3d   digit entries %5d   first digit entry %5d \n",
	   rofRec.getROFrame(), rofRec.getNEntries(), rofRec.getFirstEntry());
    int ndig = rofRec.getNEntries();
    int fdig = rofRec.getFirstEntry();
    for (int idig = fdig; idig < (fdig + ndig); idig++) {
      const auto digit = digVec[idig];
      const auto& label = oldlabelformat ? (labels->getLabels(idig)[0]) : (constlabels.getLabels(idig)[0]);

      if (!label.isNoise()) {
        srcID = label.getSourceID();
        trkID = label.getTrackID();
        evnID = label.getEventID();
	      printf("Digit %5d   in chip %03d   evn %2d   trk %4d \n", idig, digit.getChipIndex(), evnID, trkID);
      }
    }
  }
  
}
