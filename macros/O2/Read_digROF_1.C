#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <TFile.h>
#include <TTree.h>

#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/Digit.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"

#endif

void Read_digROF_1() {

  // read digits from the vector + read MC labels
  
  // DataFormats/Detectors/ITSMFT/common/include/DataFormatsITSMFT/Digit.h
  using Digit = o2::itsmft::Digit;

  TFile* fileD = TFile::Open("mftdigits.root");
  TTree *digTree = (TTree*)fileD->Get("o2sim");
  std::vector<o2::itsmft::Digit> digVec, *digVecP = &digVec;
  digTree->SetBranchAddress("MFTDigit", &digVecP);
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* labels = nullptr;
  if (digTree->GetBranch("MFTDigitMCTruth")) {
    digTree->SetBranchAddress("MFTDigitMCTruth", &labels);
  } else {
    printf("No Monte-Carlo information in this file\n");
    return;
  }
  
  int nEntries = digTree->GetEntries();
  printf("Number of entries in digits tree %d \n", nEntries);

  digTree->GetEntry(0);

  int nDigits = digVec.size();
  printf("Number of digits %d \n", nDigits);
  
  int srcID, trkID, evnID;
  bool fake;
  for (int idig = 0; idig < nDigits; idig++) {
    auto digit = digVec[idig];
    auto& label = (labels->getLabels(idig))[0];
    if (!label.isNoise()) {
      label.get(trkID, evnID, srcID, fake);
      printf("Digit %5d   in chip %03d   evn %2d   trk %4d \n", idig, digit.getChipIndex(), evnID, trkID);
    }
  }
  
}
