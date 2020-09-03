#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <TFile.h>
#include <TTree.h>

#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsITSMFT/ClusterPattern.h"

#endif

void Read_CompClsDict() {

  // read the entries from the cluster dictionary
  
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

  int nEntries = dict.getSize();
  printf("Found %d entries in the dictionary.\n", nEntries);

  unsigned long hash;
  float errX, errY, xCOG, yCOG;
  int nPixels;
  bool isGroup;
  double freq;
  o2::itsmft::ClusterPattern clsPatt;
  for (int ipatt = 0; ipatt < nEntries; ++ipatt) {
    hash = dict.getHash(ipatt);
    errX = dict.getErrX(ipatt);
    errY = dict.getErrZ(ipatt);
    xCOG = dict.getXCOG(ipatt);
    yCOG = dict.getZCOG(ipatt);
    nPixels = dict.getNpixels(ipatt);
    freq = dict.getFrequency(ipatt);
    isGroup = dict.isGroup(ipatt);
    clsPatt = dict.getPattern(ipatt);
    if (!isGroup) continue;
    printf("Pattern %5d is:\n", ipatt);
    printf("Hash: %20ld  ErrX: %6.4f ErrY: %6.4f XCOG: %7.4f YCOG: %7.4f Freq: %5.3f isGroup %d \n", hash, errX, errY, xCOG, yCOG, freq, isGroup);
    std::cout << clsPatt;
    //if (ipatt > 10) break;
  }
  
}
