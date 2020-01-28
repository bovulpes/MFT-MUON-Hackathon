void runSimulation(Int_t nevents=1,
                   Int_t runNumber=169099) {
  
  // AliLog::SetGlobalDebugLevel(1);
  
  AliSimulation *simulator = new AliSimulation();
  gSystem->Load("libFITbase.so");
  gSystem->Load("libFITsim.so");
  
  TDatime dt;
  UInt_t seed = dt.Get();
  
  simulator->SetSeed(seed);
  simulator->SetRunNumber(runNumber);
  simulator->SetMakeDigits("MUON MFT");
  //simulator->SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON ZDC PMD T0 VZERO FMD MFT");
  simulator->SetMakeSDigits("MUON MFT");
  simulator->SetMakeDigitsFromHits("ITS");
  simulator->SetRunQA(":");
  simulator->SetRunHLT("");
  
  gRandom->SetSeed(seed);
  simulator->SetDefaultStorage("local:///opt/alice/AliRootOCDB/OCDB");
    
  simulator->SetSpecificStorage("GRP/GRP/Data",
                                Form("local://%s",gSystem->pwd()));
  
  TStopwatch timer;
  timer.Start();
  simulator->Run(nevents);
  timer.Stop();
  timer.Print();
  AliCodeTimer::Instance()->Print();

}
