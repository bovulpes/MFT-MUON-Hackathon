void runReconstruction() {
  
  TDatime dt;
  UInt_t seed = dt.Get();
  
  gRandom->SetSeed(seed);
  
  AliReconstruction *reco = new AliReconstruction("galice.root");
  
  // switch off cleanESD
  reco->SetCleanESD(kFALSE);

  reco->SetDefaultStorage("local:///opt/alice/AliRootOCDB/OCDB");
  
  // GRP from local OCDB
  reco->SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd()));

  reco->SetOption("MUON MFT","SAVEDIGITS");
  reco->SetRunQA(":");
  
  reco->SetWriteESDfriend(kFALSE);
  reco->SetStopOnError(kFALSE);
  TStopwatch timer;
  timer.Start();
  reco->Run();
  timer.Stop();
  timer.Print();
  
  delete reco;
  AliCodeTimer::Instance()->Print();

}
