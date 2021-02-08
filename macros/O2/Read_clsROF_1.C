#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <TFile.h>
#include <TTree.h>

#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "MFTBase/GeometryTGeo.h"
#include "MathUtils/Cartesian3D.h"
#include "MathUtils/Utils.h"

#endif

void Read_clsROF_1() {

  // read clusters from the vector + read MC labels

  // class CompClusterExt : public CompCluster
  // DataFormats/Detectors/ITSMFT/common/include/DataFormatsITSMFT/CompCluster.h
  using o2::itsmft::CompClusterExt;

  // Geometry and matrix transformations
  std::string inputGeom = "o2sim_geometry.root";
  o2::base::GeometryManager::loadGeometry(inputGeom);
  auto gman = o2::mft::GeometryTGeo::Instance();
  gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
  
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

  TFile fileC("mftclusters.root");
  TTree *clsTree = (TTree*)fileC.Get("o2sim");
  std::vector<CompClusterExt> clsVec, *clsVecP = &clsVec;
  clsTree->SetBranchAddress("MFTClusterComp", &clsVecP);
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* labels = nullptr;
  o2::dataformats::IOMCTruthContainerView* labelROOTbuffer = nullptr;
  o2::dataformats::ConstMCTruthContainer<o2::MCCompLabel> constlabels;
  // for backward compatibility we check what is stored in the file
  auto labelClass = clsTree->GetBranch("MFTClusterMCTruth")->GetClassName();
  bool oldlabelformat = false;
  if (clsTree->GetBranch("MFTClusterMCTruth")) {
    if (TString(labelClass).Contains("IOMCTruth")) {
      // new format
      clsTree->SetBranchAddress("MFTClusterMCTruth", &labelROOTbuffer);
    } else {
      // old format
     clsTree->SetBranchAddress("MFTClusterMCTruth", &labels);
     oldlabelformat = true;
    }
  } else {
    printf("No Monte-Carlo information in this file\n");
    return;
  }

  int nEntries = clsTree->GetEntries();
  printf("Number of entries in clusters tree %d \n", nEntries);
  
  clsTree->GetEntry(0);
  if (!oldlabelformat) {
    labelROOTbuffer->copyandflatten(constlabels);
  }
  int nClusters = clsVec.size();
  printf("Number of clusters %d \n", nClusters);
  
  int srcID, trkID, evnID;
  //bool fake;
  for (int icls = 0; icls < nClusters; icls++) {
    auto cluster = clsVec[icls];
    const auto& label = oldlabelformat ? (labels->getLabels(icls)[0]) : (constlabels.getLabels(icls)[0]);
    //auto& label = (labels->getLabels(icls))[0];
    if (!label.isNoise()) {
      //label.get(trkID, evnID, srcID, fake);
      srcID = label.getSourceID();
      trkID = label.getTrackID();
      evnID = label.getEventID();
      auto chipID = cluster.getChipID(); 
      auto pattID = cluster.getPatternID();
      o2::math_utils::Point3D<float> locC;
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

  fileC.Close();
  
}
