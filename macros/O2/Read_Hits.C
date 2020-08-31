#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <TFile.h>
#include <TTree.h>

#include "ITSMFTSimulation/Hit.h"

#endif

void Read_Hits() {

  using o2::itsmft::Hit;

  // class Hit : public o2::BasicXYZEHit<Float_t, Float_t>
  // class BasicXYZEHit : public BasicXYZVHit<T, E, V>
  // class BasicXYZVHit : public BaseHit
  // DataFormats/simulation/include/SimulationDataFormat/BaseHits.h
  using HitVec = std::vector<Hit>;

  TFile fileH("o2sim_HitsMFT.root");
  TTree* hitTree = (TTree*)fileH.Get("o2sim");
  // to allow the "automatic" loop over the hits
  // for (auto hit : hitArray) { ... }
  std::vector<Hit> hitVec, *hitVecP = &hitVec;
  hitTree->SetBranchAddress("MFTHit", &hitVecP);
  
  int nEvents = hitTree->GetEntries();
  printf("Number of events %d \n", nEvents);

  for (int iev = 0; iev < nEvents; ++iev) {
    hitTree->GetEntry(iev);
    int nHits = hitVec.size();
    printf("Event %d has %5d hits\n", iev, nHits);
    int nh = 0;
    for (auto hit : hitVec) {
      printf("Hit %5d   trkID %4d   detID %03d   x,y,z  %7.3f  %7.3f  %7.3f \n", nh, hit.GetTrackID(), hit.GetDetectorID(), hit.GetX(), hit.GetY(), hit.GetZ());
      ++nh;
    }
  }
  
  fileH.Close();
  
}
