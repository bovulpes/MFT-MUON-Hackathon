#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <TFile.h>
#include <TTree.h>

#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"

#endif

void Read_clsROF_2() {

  // read clusters by Read-Out-Frame

  // class CompClusterExt : public CompCluster
  // DataFormats/Detectors/ITSMFT/common/include/DataFormatsITSMFT/CompCluster.h
  using o2::itsmft::CompClusterExt;
  using ROFRec = o2::itsmft::ROFRecord;

  TFile fileC("mftclusters.root");
  TTree *clsTree = (TTree*)fileC.Get("o2sim");
  std::vector<CompClusterExt> clsVec, *clsVecP = &clsVec;
  clsTree->SetBranchAddress("MFTClusterComp", &clsVecP);
  std::vector<o2::itsmft::MC2ROFRecord> mc2rofVec, *mc2rofVecP = &mc2rofVec;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* labels = nullptr;
  if (clsTree->GetBranch("MFTClusterMCTruth")) {
    clsTree->SetBranchAddress("MFTClusterMCTruth", &labels);
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

  int srcID, trkID, evnID;
  bool fake;
  for (int irof = 0; irof < nROFRec; irof++) {
    const auto& rofRec = rofRecVec[irof];
    //rofRec.print();
    printf("ROFRecord ID %3d cluster entries %4d first cluster entry %4d \n",
	   irof, rofRec.getNEntries(), rofRec.getFirstEntry());
    int ncls = rofRec.getNEntries();
    int fcls = rofRec.getFirstEntry();
    for (int icls = fcls; icls < (fcls + ncls); icls++) {
      auto cluster = clsVec[icls];
      auto& label = (labels->getLabels(icls))[0];
      if (!label.isNoise()) {
	label.get(trkID, evnID, srcID, fake);
	printf("Cluster %5d   chip ID %03d   evn %2d   mctrk %4d \n", icls, cluster.getChipID(), evnID, trkID);
      }
    }
  }

  fileC.Close();
  
}
