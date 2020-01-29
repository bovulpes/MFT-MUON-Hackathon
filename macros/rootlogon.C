{
  cout << "Loading libraries ..." << endl;

  gSystem->Load("libVMC");
  gSystem->Load("libMinuit");
  gSystem->Load("libTree");
  gSystem->Load("libProofPlayer");
  gSystem->Load("libXMLParser");
  gSystem->Load("libPhysics");

  gSystem->Load("libGui");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libCDB");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libSTEER");

  //gSystem->Load("libEVGEN");
  
  cout << "Setting include path ..." << endl;
  TString includePath = "-I${ALICE_ROOT}/STEER ";
  includePath        += "-I${ALICE_ROOT}/STEER/STEER ";
  includePath        += "-I${ALICE_ROOT}/STEER/STEERBase ";
  includePath        += "-I${ALICE_ROOT}/STEER/CDB ";
  includePath        += "-I${ALICE_ROOT}/STEER/ESD ";
  includePath        += "-I${ALICE_ROOT}/RAW ";
  includePath        += "-I${ALICE_ROOT}/FASTSIM ";
  includePath        += "-I${ALICE_ROOT}/EVGEN ";
  includePath        += "-I${ALICE_ROOT}/SHUTTLE/TestShuttle ";
  includePath        += "-I${ALICE_ROOT}/ITS ";
  includePath        += "-I${ALICE_ROOT}/MUON ";
  includePath        += "-I${ALICE_ROOT}/MUON/mapping ";
  includePath        += "-I${ALICE_ROOT}/RAW ";
  includePath        += "-I${ALICE_ROOT}/include ";
  includePath        += "-I${ALICE_ROOT}/MFT ";
  
  gSystem->SetIncludePath(includePath.Data());
}
