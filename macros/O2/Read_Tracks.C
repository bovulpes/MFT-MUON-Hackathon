#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <TFile.h>
#include <TTree.h>

#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "MFTBase/GeometryTGeo.h"
#include "MathUtils/Cartesian3D.h"
#include "MathUtils/Utils.h"

#endif

void Read_Tracks() {

  // class CompClusterExt : public CompCluster
  // DataFormats/Detectors/ITSMFT/common/include/DataFormatsITSMFT/CompCluster.h
  using o2::itsmft::CompClusterExt;

  // Geometry and matrix transformations
  std::string inputGeom = "o2sim_geometry.root";
  o2::base::GeometryManager::loadGeometry(inputGeom);
  auto gman = o2::mft::GeometryTGeo::Instance();
  gman->fillMatrixCache(o2::utils::bit2Mask(o2::TransformType::L2G));
  
  // Cluster pattern dictionary
  std::string dictfile = "MFTdictionary.bin";
  o2::itsmft::TopologyDictionary dict;
  std::ifstream file(dictfile.c_str());
  if (file.good()) {
    printf("Running with dictionary: %s \n", dictfile.c_str());
    dict.readBinaryFile(dictfile);
  } else {
    printf("Can not run without dictionary !\n");
    return;
  }

  // Clusters
  
  TFile fileC("mftclusters.root");
  TTree *clsTree = (TTree*)fileC.Get("o2sim");
  std::vector<CompClusterExt> clsVec, *clsVecP = &clsVec;
  clsTree->SetBranchAddress("MFTClusterComp", &clsVecP);
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* clsLabels = nullptr;
  if (clsTree->GetBranch("MFTClusterMCTruth")) {
    clsTree->SetBranchAddress("MFTClusterMCTruth", &clsLabels);
  } else {
    printf("No Monte-Carlo information in this file\n");
    return;
  }

  int nEntries = clsTree->GetEntries();
  printf("Number of entries in clusters tree %d \n", nEntries);
  
  clsTree->GetEntry(0);

  int nClusters = clsVec.size();
  printf("Number of clusters %d \n", nClusters);
  
  // Tracks
  // class TrackMFT : public o2::track::TrackParCovFwd
  // class TrackParCovFwd : public TrackParFwd
  // DataFormats/Detectors/ITSMFT/MFT/include/DataFormatsMFT/TrackMFT.h
  // DataFormats/Reconstruction/include/ReconstructionDataFormats/TrackFwd.h
  TFile fileT("mfttracks.root");
  TTree *trackTree = (TTree*)fileT.Get("o2sim");
  std::vector<o2::mft::TrackMFT> trackVec, *trackVecP = &trackVec;
  trackTree->SetBranchAddress("MFTTrack", &trackVecP);
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* trkLabels = nullptr;
  if (trackTree->GetBranch("MFTTrackMCTruth")) {
    trackTree->SetBranchAddress("MFTTrackMCTruth", &trkLabels);
  } else {
    printf("No Monte-Carlo information in this file\n");
    return;
  }
    
  std::vector<int> trackExtClsVec, *trackExtClsVecP = &trackExtClsVec;
  trackTree->SetBranchAddress("MFTTrackClusIdx", &trackExtClsVecP);

  trackTree->GetEntry(0);

  int srcID, trkID, evnID;
  bool fake;
  int iTrack = 0;
  for (auto &track : trackVec) {
    auto trkX = track.getX();
    auto trkY = track.getY();
    auto trkZ = track.getZ();
    auto trkLabel = trkLabels->getLabels(iTrack);
    //trkLabel[0].print();
    auto eventID = trkLabel[0].getEventID();
    auto outParam = track.getOutParam();
    auto trkOutX = outParam.getX();
    auto trkOutY = outParam.getY();
    auto trkOutZ = outParam.getZ();
    printf("Track %3d   isCA %1d   x,y,z-in  %7.3f  %7.3f  %7.3f  x,y,z-out  %7.3f  %7.3f  %7.3f   ev %2d  labels ", iTrack, track.isCA(), trkX, trkY, trkZ, trkOutX, trkOutY, trkOutZ, eventID);
    for (auto ilab = 0; ilab < trkLabel.size(); ++ilab) {
      auto trkID = trkLabel[ilab].getTrackID();
      printf("%4d   ", trkID);
    }
    printf("\n");
    auto ncls = track.getNumberOfPoints();
    auto offset = track.getExternalClusterIndexOffset();
    printf("Number of clusters %2d \n", ncls);
    for (int icls = 0; icls < ncls; ++icls) {
      auto clsEntry = trackExtClsVec[offset + icls];
      auto cluster = clsVec[clsEntry];
      auto& clsLabel = (clsLabels->getLabels(clsEntry))[0];
      if (!clsLabel.isNoise()) {
	clsLabel.get(trkID, evnID, srcID, fake);
	auto chipID = cluster.getChipID(); 
	auto pattID = cluster.getPatternID();
	Point3D<float> locC;
	int npix = 0;
	if (pattID == o2::itsmft::CompCluster::InvalidPatternID || dict.isGroup(pattID)) {
	  // temporary fix ...
	  printf("temporary fix for group pattern ...\n");
	  locC = dict.getClusterCoordinates(cluster);
	  
	  //o2::itsmft::ClusterPattern patt(pattIt);
	  //locC = dict.getClusterCoordinates(cluster, patt);
	} else {
	  locC = dict.getClusterCoordinates(cluster);
	  npix = dict.getNpixels(pattID);
	}
	// Transformation to the local --> global
	auto gloC = gman->getMatrixL2G(chipID) * locC;
	printf("Cluster %5d   chip ID %03d   evn %2d   mctrk %4d   x,y,z  %7.3f  %7.3f  %7.3f \n", icls, cluster.getChipID(), evnID, trkID, gloC.X(), gloC.Y(), gloC.Z());
      }
    }
    ++iTrack;
  }
  
  fileT.Close();
  fileC.Close();

}

