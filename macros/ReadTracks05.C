#include <fstream>
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TGeoManager.h>
#include <TVector3.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TList.h>
#include <TMatrix.h>

#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliMFTCATrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrack.h"
#include "AliMuonForwardTrack.h"
#include "AliMUONESDInterface.h"
#include "AliMUONTrackExtrap.h"
#include "AliCDBManager.h"
#include "AliMUONCDB.h"
#include "AliMUONConstants.h"
#include "AliESDVertex.h"
#include "AliStack.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliHeader.h"
#include "AliMFTCACell.h"

void ReadTracks05(Int_t listNr = 7, Int_t runNumber = 169099) {

  TString baseDir = TString("/home/vulpescu/alice/aliroot/work/sim/04/1");

  Bool_t useVertexFromESD = kFALSE; // if not use saved MC vertex
 
  Bool_t test01 = kFALSE;

  Bool_t prn = kFALSE;

  Bool_t useStack = kTRUE;

  Bool_t onlyRef = kFALSE;

  if (test01) useStack = kTRUE;

  Double_t mmChi2Cut = 2.0;
  Double_t mmDcaTheXcut = 2.5;
  Double_t mmDcaTheYcut = 2.5;
  Double_t mmVtxTheXcut = 0.6;
  Double_t mmVtxTheYcut = 0.6;
  /*
  // no cuts
  mmDcaTheXcut = 180.0;
  mmDcaTheYcut = 180.0;
  mmVtxTheXcut = 180.0;
  mmVtxTheYcut = 180.0;
  */
  //
  Short_t muonQ, mftQ, trackQ, chi2minMuonQ, chi2minMftQ;
  Double_t fXVertex=0;
  Double_t fYVertex=0;
  Double_t fZVertex=0;
  Double_t errXVtx=0;
  Double_t errYVtx=0;
  Int_t nMisQMuonMft = 0;
  Int_t nEvents, nparticles, nprimarypart, ntracks, nTracks;
  Int_t nMUONTracks, nMFTTracks, nMFTCells, iEventS, iEventE;
  Double_t xVertex, yVertex, zVertex;
  Double_t xVertexMC, yVertexMC, zVertexMC;
  Double_t xESDTrack, yESDTrack, zESDTrack;
  Double_t xDCAESDTrack, yDCAESDTrack, zDCAESDTrack;
  Double_t esdTrackPxyzVtx[3], esdTrackPVtx, esdTrackPtVtx;
  Double_t esdTrackPxyzVtx1[3], esdTrackPVtx1, esdTrackPtVtx1;
  Double_t esdTrackTheVtx, esdTrackPhiVtx;
  Double_t esdTrackTheVtx1, esdTrackPhiVtx1;
  Double_t esdTrackPxyzDca[3], esdTrackPDca, esdTrackPtDca;
  Double_t esdTrackPxyzDca1[3], esdTrackPDca1, esdTrackPtDca1;
  Double_t esdTrackPxyzUnc[3], esdTrackPUnc;
  Double_t mmTrackPxyzDca[3], mmTrackPDca, mmTrackPtDca;
  Double_t mmTrackPxyzDca1[3], mmTrackPDca1, mmTrackPtDca1;
  Double_t stackPxyzVtx[3], stackPVtx, stackPtVtx;
  Double_t stackPxyzVtx1[3], stackPVtx1, stackPtVtx1;
  Double_t stackXyzVtx[3], stackXyzVtx1[3];
  Double_t mmTrackXyzDca[3], mmTrackXyzDca1[3];
  Double_t esdTrackXyzVtx[3], esdTrackXyzDca[3], esdTrackXyzUnc[3];
  Double_t mftTrackEnd[3], esdTrackXyzVtxMft[3], esdTrackXyzDcaMft[3];
  Double_t mftTrackTheX, mftTrackTheY, theXyDifX, theXyDifY, theXyDif;
  Double_t theXyDifMin;
  Bool_t matchMM_MC, matchMM_Q, matchMM_Rxy, foundMCLabelMFT;
  Bool_t esdToMuonTrack, chi2minUpdate, matchMM_DcaTheXY, matchMM_VtxTheXY;
  Double_t esdTrackDcaSqErrX, esdTrackDcaSqErrY;
  Double_t esdTrackVtxSqErrX, esdTrackVtxSqErrY;
  Double_t esdTrackDcaTheX, esdTrackDcaTheY, esdTrackVtxTheX, esdTrackVtxTheY;
  Double_t rxyDcaMFT, rxDcaMFT, ryDcaMFT, errxyDcaMFT;
  Double_t rxyVtxMFT, rxVtxMFT, ryVtxMFT, errxyVtxMFT;
  TMatrixD covDcaMFT(5,5), covVtxMFT(5,5);
  Double_t chi2min;
  Int_t muonMCLabel, chi2minMuonMCLabel, chi2minMftMCLabel, pdgCodeMUON, pdgCodeMFT;
  Int_t theXyDifMinMuonMCLabel, theXyDifMinMftMCLabel;
  Double_t muonTrackPxyz0[3], muonTrackPxyz1[3], muonTrackPxyz2[3];
  Double_t muonTrackPxyz3[3], muonTrackPxyz4[3];
  Double_t muonTrackPt0, muonTrackPt1, muonTrackPt2, muonTrackPt3, muonTrackPt4;
  Double_t the, eta, phi;
  TVector3 mftTrackV, esdTrackV;
  Bool_t skipit;

  Double_t frxy = 3.0;

  //

  TH1F *hPtM = new TH1F("hPtM","Pt of ESDMuon tracks",100,0.,20.);
  TH1F *hPtMM = new TH1F("hPtMM","Pt of MUON+MFT tracks",100,0.,20.);
  TH1F *hPt1s1 = new TH1F("hPt1s1","Pt for same MC label and same charge",100,0.,5.);
  TH1F *hPt1s0 = new TH1F("hPt1s0","Pt for same MC label",100,0.,5.);

  // pT from track param at first MUON cluster
  TH1F *hPt0 = new TH1F("hPt0","Pt of all MUON tracks",100,0.,5.);
  TH1F *hPt1 = new TH1F("hPt1","Pt of MUON + MFT tracks ",100,0.,5.);
  TH1F *hPt2 = new TH1F("hPt2","Pt of MUON + MFT tracks MC matched",100,0.,5.);
  TH1F *hPt3 = new TH1F("hPt3","Pt of MUON + MFT tracks MC matched for chi2min",100,0.,5.);
  TH1F *hPt4 = new TH1F("hPt4","Pt of MUON + MFT tracks MC and Q matched",100,0.,5.);
  TH1F *hPt5 = new TH1F("hPt5","Pt of MUON + MFT tracks MC for chi2min and Q matched",100,0.,5.);
  TH1F *hPt6 = new TH1F("hPt6","Pt of MUON + MFT tracks MC for theXyDifMin matched",100,0.,5.);
  TH1F *hPt7 = new TH1F("hPt7","Pt of MUON + MFT tracks MC for theXyDifMin and chi2min matched",100,0.,5.);
  TH1F *hPt8 = new TH1F("hPt8","Pt of MUON + MFT tracks not MC matched for chi2min",100,0.,5.);
  TH1F *hPt9 = new TH1F("hPt9","Pt of MUON + MFT tracks not MC matched for theXyDifMin",100,0.,5.);

  // pT at vertex (Branson)
  TH1F *hVtxPt0 = new TH1F("hVtxPt0","Pt of all MUON tracks",100,0.,5.);
  TH1F *hVtxPt1 = new TH1F("hVtxPt1","Pt of MUON + MFT tracks ",100,0.,5.);
  TH1F *hVtxPt2 = new TH1F("hVtxPt2","Pt of MUON + MFT tracks MC matched",100,0.,5.);
  TH1F *hVtxPt3 = new TH1F("hVtxPt3","Pt of MUON + MFT tracks MC matched for chi2min",100,0.,5.);
  TH1F *hVtxPt4 = new TH1F("hVtxPt4","Pt of MUON + MFT tracks MC and Q matched",100,0.,5.);
  TH1F *hVtxPt5 = new TH1F("hVtxPt5","Pt of MUON + MFT tracks MC for chi2min and Q matched",100,0.,5.);
  TH1F *hVtxPt6 = new TH1F("hVtxPt6","Pt of MUON + MFT tracks MC for theXyDifMin matched",100,0.,5.);
  TH1F *hVtxPt7 = new TH1F("hVtxPt7","Pt of MUON + MFT tracks MC for theXyDifMin and chi2min matched",100,0.,5.);
  TH1F *hVtxPt8 = new TH1F("hVtxPt8","Pt of MUON + MFT tracks not MC matched for chi2min",100,0.,5.);
  TH1F *hVtxPt9 = new TH1F("hVtxPt9","Pt of MUON + MFT tracks not MC matched for theXyDifMin",100,0.,5.);

  // pT at DCA
  TH1F *hDcaPt0 = new TH1F("hDcaPt0","Pt of all MUON tracks",100,0.,5.);
  TH1F *hDcaPt1 = new TH1F("hDcaPt1","Pt of MUON + MFT tracks ",100,0.,5.);
  TH1F *hDcaPt2 = new TH1F("hDcaPt2","Pt of MUON + MFT tracks MC matched",100,0.,5.);
  TH1F *hDcaPt3 = new TH1F("hDcaPt3","Pt of MUON + MFT tracks MC matched for chi2min",100,0.,5.);
  TH1F *hDcaPt4 = new TH1F("hDcaPt4","Pt of MUON + MFT tracks MC and Q matched",100,0.,5.);
  TH1F *hDcaPt5 = new TH1F("hDcaPt5","Pt of MUON + MFT tracks MC for chi2min and Q matched",100,0.,5.);
  TH1F *hDcaPt6 = new TH1F("hDcaPt6","Pt of MUON + MFT tracks MC for theXyDifMin matched",100,0.,5.);
  TH1F *hDcaPt7 = new TH1F("hDcaPt7","Pt of MUON + MFT tracks MC for theXyDifMin and chi2min matched",100,0.,5.);
  TH1F *hDcaPt8 = new TH1F("hDcaPt8","Pt of MUON + MFT tracks not MC matched for chi2min",100,0.,5.);
  TH1F *hDcaPt9 = new TH1F("hDcaPt9","Pt of MUON + MFT tracks not MC matched for theXyDifMin",100,0.,5.);

  TH1F *hEta0 = new TH1F("hEta0","Eta of all MUON tracks",100,-4.,-2.); 
  TH1F *hEta1 = new TH1F("hEta1","Eta of all MUON + MFT tracks",100,-4.,-2.); 

  TH1F *hChi2m0a = new TH1F("hChi2m0a","Chi2 for diff MC label",200,0,10);
  TH1F *hChi2m0b = new TH1F("hChi2m0b","Chi2 for diff MC label",200,0,10);
  TH1F *hChi2m1a = new TH1F("hChi2m1a","Chi2 for same MC label",200,0,10);
  TH1F *hChi2m1b = new TH1F("hChi2m1b","Chi2 for same MC label",200,0,10);
  TH1F *hChi2muon = new TH1F("hChi2muon","Chi2 MUON track",200,0,10);
  TH2F *hDxyVtx0 = new TH2F("hDxyVtx0","Vertex dx,dy from ESD track",100,-0.1,+0.1,100,-0.1,+0.1);
  TH2F *hDxyVtx1 = new TH2F("hDxyVtx1","Vertex dx,dy from MUON+MFT track",100,-0.1,+0.1,100,-0.1,+0.1);

  TH1F *hDP0 = new TH1F("hDP0","P diff ESD-MC",100,-5.,+5.);
  TH1F *hDP1 = new TH1F("hDP1","P diff (MUON+MFT)-MC",100,-5.,+5.);
  TH2F *hDxy0 = new TH2F("hDxy0","XY at DCA diff ESD-MC",100,-0.2,+0.2,100,-0.2,+0.2);
  TH2F *hDxy1 = new TH2F("hDxy1","XY at DCA diff (MUON+MFT)-MC",100,-0.2,+0.2,100,-0.2,+0.2);

  TH2F *hDxyDcaMft0 = new TH2F("hDxyDcaMft0","MFT dx, dy from ESD track",100,-32.,+32.,100,-32.,+32.);
  TH2F *hDxyVtxMft0 = new TH2F("hDxyVtxMft0","MFT dx, dy from ESD track",100,-32.,+32.,100,-32.,+32.);
  TH2F *hDxyDcaMft1 = new TH2F("hDxyDcaMft1","MFT dx, dy from ESD track",100,-32.,+32.,100,-32.,+32.);
  TH2F *hDxyVtxMft1 = new TH2F("hDxyVtxMft1","MFT dx, dy from ESD track",100,-32.,+32.,100,-32.,+32.);
  TH1F *hErrxyMFT = new TH1F("hErrxyMFT","Errxy (R) diff",100,-10.,+10.);
  TH2F *hErrxyDifMFT0 = new TH2F("hErrxyDifMFT0","Errxy vs diff",100,0.,10.,100,0.,10.);
  TH2F *hErrxyDifMFT1 = new TH2F("hErrxyDifMFT1","Errxy vs diff",100,0.,10.,100,0.,10.);

  TH2F *hTheXyDcaDif0 = new TH2F("hTheXyDcaDif0","MUON MFT diff thex=NPB, they=BP (DCA)",100,-5,+5,100,-5,+5);
  TH2F *hTheXyDcaDif1 = new TH2F("hTheXyDcaDif1","MUON MFT diff thex=NPB, they=BP (DCA)",100,-5,+5,100,-5,+5);
  TH2F *hTheXyVtxDif0 = new TH2F("hTheXyVtxDif0","MUON MFT diff they=NBP, they=BP (VTX)",100,-5,+5,100,-5,+5);
  TH2F *hTheXyVtxDif1 = new TH2F("hTheXyVtxDif1","MUON MFT diff they=NBP, they=BP (VTX)",100,-5,+5,100,-5,+5);

  TH2F *hTheXyVtxDcaDif0 = new TH2F("hTheXyVtxDcaDif0","Correlation thex",100,-5,+5,100,-5,+5);
  TH2F *hTheXyVtxDcaDif1 = new TH2F("hTheXyVtxDcaDif1","Correlation thex",100,-5,+5,100,-5,+5);

  TH3F *hEff0 = new TH3F("hEff0","Eff muons from JPsi, 0",10,0.,10.,10,170.,177.,10,0.,360.);
  TH3F *hEff1 = new TH3F("hEff1","Eff muons from JPsi, 1",10,0.,10.,10,170.,177.,10,0.,360.);

  TH1F *hMuPt0 = new TH1F("hMuPt0","Pt of all MUON muon tracks",200,0.,5.);
  TH1F *hMuEta0 = new TH1F("hMuEta0","Eta of all MUON muon tracks",100,-4.,-2.); 

  TList *histList = new TList();

  histList->Add(hPt1s0);
  histList->Add(hPt1s1);
  histList->Add(hChi2m0a);
  histList->Add(hChi2m0b);
  histList->Add(hChi2m1a);
  histList->Add(hChi2m1b);
  histList->Add(hChi2muon);
  histList->Add(hDxyVtx0);
  histList->Add(hDxyVtx1);
  histList->Add(hDP0);
  histList->Add(hDP1);
  histList->Add(hDxy0);
  histList->Add(hDxy1);
  histList->Add(hDxyDcaMft0);
  histList->Add(hDxyVtxMft0);
  histList->Add(hDxyDcaMft1);
  histList->Add(hDxyVtxMft1);
  histList->Add(hErrxyMFT);
  histList->Add(hErrxyDifMFT0);
  histList->Add(hErrxyDifMFT1);
  histList->Add(hPtM);
  histList->Add(hPtMM);
  histList->Add(hPt0);
  histList->Add(hPt1);
  histList->Add(hPt2);
  histList->Add(hPt3);
  histList->Add(hPt4);
  histList->Add(hPt5);
  histList->Add(hPt6);
  histList->Add(hPt7);
  histList->Add(hPt8);
  histList->Add(hPt9);
  histList->Add(hDcaPt0);
  histList->Add(hDcaPt1);
  histList->Add(hDcaPt2);
  histList->Add(hDcaPt3);
  histList->Add(hDcaPt4);
  histList->Add(hDcaPt5);
  histList->Add(hDcaPt6);
  histList->Add(hDcaPt7);
  histList->Add(hDcaPt8);
  histList->Add(hDcaPt9);
  histList->Add(hVtxPt0);
  histList->Add(hVtxPt1);
  histList->Add(hVtxPt2);
  histList->Add(hVtxPt3);
  histList->Add(hVtxPt4);
  histList->Add(hVtxPt5);
  histList->Add(hVtxPt6);
  histList->Add(hVtxPt7);
  histList->Add(hVtxPt8);
  histList->Add(hVtxPt9);
  histList->Add(hEta0);
  histList->Add(hEta1);
  histList->Add(hTheXyDcaDif0);
  histList->Add(hTheXyDcaDif1);
  histList->Add(hTheXyVtxDif0);
  histList->Add(hTheXyVtxDif1);
  histList->Add(hTheXyVtxDcaDif0);
  histList->Add(hTheXyVtxDcaDif1);
  histList->Add(hEff0);
  histList->Add(hEff1);
  histList->Add(hMuPt0);
  histList->Add(hMuEta0);

  Float_t pMax = 50., ptMax = 5., theMin = 167., theMax = 180.;

  TH1F *hMftThe = new TH1F("hMftThe","Theta of all MFT tracks",200,theMin,theMax); 
  TH1F *hMftPhi = new TH1F("hMftPhi","Phi of all MFT tracks",100,0.,360.); 

  histList->Add(hMftThe);
  histList->Add(hMftPhi);

  //
  // references
  //
  TH1F *hRefPtot = new TH1F("hRefPtot","hRefPtot",200,0.,pMax);
  TH1F *hRefPt = new TH1F("hRefPt","hRefPt",200,0.,ptMax);
  TH1F *hRefThe = new TH1F("hRefThe","hRefThe",200,theMin,theMax);
  TH1F *hRefPhi = new TH1F("hRefPhi","hRefPhi",100,0.,360.);
  TH1F *hRefRap = new TH1F("hRefRap","hRefRap",500,-4.5,-0.0);
  TH1F *hRefEta = new TH1F("hRefEta","hRefEta",500,-6.0,-2.0);

  TH1F *hRefPtot1 = new TH1F("hRefPtot1","hRefPtot1",200,0.,pMax);
  TH1F *hRefPt1 = new TH1F("hRefPt1","hRefPt1",200,0.,ptMax);
  TH1F *hRefThe1 = new TH1F("hRefThe1","hRefThe1",200,theMin,theMax);
  TH1F *hRefPhi1 = new TH1F("hRefPhi1","hRefPhi1",100,0.,360.);
  TH1F *hRefRap1 = new TH1F("hRefRap1","hRefRap1",500,-4.5,-0.0);
  TH1F *hRefEta1 = new TH1F("hRefEta1","hRefEta1",500,-6.0,-2.0);

  TList *histList1 = new TList();
  histList1->Add(hRefPtot);
  histList1->Add(hRefPt);
  histList1->Add(hRefThe);
  histList1->Add(hRefPhi);
  histList1->Add(hRefRap);
  histList1->Add(hRefEta);
  histList1->Add(hRefPtot1);
  histList1->Add(hRefPt1);
  histList1->Add(hRefThe1);
  histList1->Add(hRefPhi1);
  histList1->Add(hRefRap1);
  histList1->Add(hRefEta1);

  double etaMin = -3.7, etaMax = -2.4;
  double phiMin = 0., phiMax = TMath::TwoPi();
  double thegMin, thegMax;

  thegMin = TMath::Pi()+2.*TMath::ATan(-TMath::Exp(etaMax));
  thegMax = TMath::Pi()+2.*TMath::ATan(-TMath::Exp(etaMin));

  // x = eta, y = phi
  int npx = 500, npy = 500;

  TH2F *hEtaPhi = new TH2F("hEtaPhi","",npx,phiMin*TMath::RadToDeg(),phiMax*TMath::RadToDeg(),npy,etaMin,etaMax);
  TH2F *hThePhi = new TH2F("hThePhi","",npx,phiMin*TMath::RadToDeg(),phiMax*TMath::RadToDeg(),npy,thegMin*TMath::RadToDeg(),thegMax*TMath::RadToDeg());

  TH2F *hEtaPhi1 = new TH2F("hEtaPhi1","",npx,phiMin*TMath::RadToDeg(),phiMax*TMath::RadToDeg(),npy,etaMin,etaMax);
  TH2F *hThePhi1 = new TH2F("hThePhi1","",npx,phiMin*TMath::RadToDeg(),phiMax*TMath::RadToDeg(),npy,thegMin*TMath::RadToDeg(),thegMax*TMath::RadToDeg());

  histList1->Add(hEtaPhi);
  histList1->Add(hThePhi);
  histList1->Add(hEtaPhi1);
  histList1->Add(hThePhi1);

  TH2F *hMftEtaPhi = new TH2F("hMftEtaPhi","",npx,phiMin*TMath::RadToDeg(),phiMax*TMath::RadToDeg(),npy,etaMin,etaMax);
  TH2F *hMftThePhi = new TH2F("hMftThePhi","",npx,phiMin*TMath::RadToDeg(),phiMax*TMath::RadToDeg(),npy,thegMin*TMath::RadToDeg(),thegMax*TMath::RadToDeg());

  histList->Add(hMftEtaPhi);
  histList->Add(hMftThePhi);

  TTree *esdTree;
  TTree *mftTree;
  TTree *cellTree;
  TTree *muonTree;
  TTree *eventTree;
  TClonesArray *mftTracks = new TClonesArray("AliMFTCATrack");
  TClonesArray *mftCells = new TClonesArray("AliMFTCACell");
  TClonesArray *muonTracks = new TClonesArray("AliMuonForwardTrack");
  AliMUONTrack *muonTrack = 0x0;
  AliMUONTrackParam trackParamMM, trackParamMM0vtx, trackParamMM0dca, trackParamM;
  TParticle *particle;
  AliStack *theStack;
  AliESDVertex *Vertex;
  AliMuonForwardTrack *mmTrack;
  AliMFTCATrack *mftTrack;

  //

  if (!gGeoManager) {
    TGeoManager::Import(Form("%s/geometry.root",baseDir.Data()));
    if (!gGeoManager) {
      printf("getting geometry from file failed");
      return;
    }
  }
  AliCDBManager::Instance()->SetDefaultStorage("local:///opt/alice/AliRootOCDB/OCDB");
  AliCDBManager::Instance()->SetSpecificStorage("GRP/GRP/Data",
			     Form("local://%s",baseDir.Data()));
  AliCDBManager::Instance()->SetSpecificStorage("MUON/Calib/RecoParam", "local:///home/vulpescu/alice/aliroot/OCDB");
  AliCDBManager::Instance()->SetSpecificStorage("MFT/Calib/RecoParam",  "local:///home/vulpescu/alice/aliroot/OCDB");
  AliCDBManager::Instance()->SetRun(runNumber);
  if (!AliMUONCDB::LoadField()) return;
  AliMUONTrackExtrap::SetField();

  AliRunLoader *runLoader;
  TFile *esdFile, *mftFile, *refFile;

  // batch loop

  ifstream in;
  //in.open(Form("ListSimsNewSet%02d_short.txt",listNr),ios::in);

  Char_t line[1024];
  Int_t iFile = 0;

  for (Int_t i = 1; i <= 10; i++) {
  /*
  while (in.getline(line,1024)) {

  baseDir = TString(line);
  */
  baseDir = TString(Form("/home/vulpescu/alice/aliroot/work/rec/%02d/%d",listNr,i));
  if (onlyRef) {
    baseDir = TString(Form("/home/vulpescu/alice/aliroot/work/sim/%02d/%d",listNr,i));
  }

  printf("Open %s \n",baseDir.Data());
  iFile++;

  if (test01) if (iFile == 2) break;

  if (onlyRef) {
    refFile = TFile::Open(Form("%s/readmfthits.root",baseDir.Data()));
    //refFile = TFile::Open(Form("%s/readmfthits_muons.root",baseDir.Data()));
    hRefPtot->Add((TH1F*)refFile->Get("hLPtot"));
    hRefPtot->Add((TH1F*)refFile->Get("hHPtot"));
    hRefPtot1->Add((TH1F*)refFile->Get("hLPtot1"));
    hRefPtot1->Add((TH1F*)refFile->Get("hHPtot1"));
    hRefPt->Add((TH1F*)refFile->Get("hLPt"));
    hRefPt->Add((TH1F*)refFile->Get("hHPt"));
    hRefPt1->Add((TH1F*)refFile->Get("hLPt1"));
    hRefPt1->Add((TH1F*)refFile->Get("hHPt1"));
    hRefThe->Add((TH1F*)refFile->Get("hLThe"));
    hRefThe->Add((TH1F*)refFile->Get("hHThe"));
    hRefThe1->Add((TH1F*)refFile->Get("hLThe1"));
    hRefThe1->Add((TH1F*)refFile->Get("hHThe1"));
    hRefPhi->Add((TH1F*)refFile->Get("hLPhi"));
    hRefPhi->Add((TH1F*)refFile->Get("hHPhi"));
    hRefPhi1->Add((TH1F*)refFile->Get("hLPhi1"));
    hRefPhi1->Add((TH1F*)refFile->Get("hHPhi1"));
    hRefRap->Add((TH1F*)refFile->Get("hLRap"));
    hRefRap->Add((TH1F*)refFile->Get("hHRap"));
    hRefRap1->Add((TH1F*)refFile->Get("hLRap1"));
    hRefRap1->Add((TH1F*)refFile->Get("hHRap1"));
    hRefEta->Add((TH1F*)refFile->Get("hLEta"));
    hRefEta->Add((TH1F*)refFile->Get("hHEta"));
    hRefEta1->Add((TH1F*)refFile->Get("hLEta1"));
    hRefEta1->Add((TH1F*)refFile->Get("hHEta1"));
    hEtaPhi->Add((TH2F*)refFile->Get("hEtaPhi"));
    hThePhi->Add((TH2F*)refFile->Get("hThePhi"));
    hEtaPhi1->Add((TH2F*)refFile->Get("hEtaPhi1"));
    hThePhi1->Add((TH2F*)refFile->Get("hThePhi1"));
    continue;

  }

  if (useStack) {
    runLoader = AliRunLoader::Open(Form("%s/generated/galice.root",baseDir.Data()));
    if (!runLoader) {
      printf("Getting run loader from file failed\n");
      return;
    }
    runLoader->LoadgAlice();
    gAlice = runLoader->GetAliRun();
    if (!gAlice) {
      printf("no galice object found\n");
      return;
    }
    runLoader->LoadHeader();
    if (runNumber != runLoader->GetHeader()->GetRun()) {
      printf("Mismatch between run number from ESD and from runLoader\n");
      return;
    }
    Int_t nevents = runLoader->GetNumberOfEvents();
    runLoader->LoadKinematics("READ");
  }

  //

  esdFile = TFile::Open(Form("%s/AliESDs.root",baseDir.Data()));
  mftFile = TFile::Open(Form("%s/MFT.Tracks.root",baseDir.Data()));
  esdTree = (TTree*) esdFile->Get("esdTree");
  AliESDEvent * esd = new AliESDEvent;
  esd->ReadFromTree(esdTree);
  printf("ESD Entries %d \n",(Int_t)esdTree->GetEntries());
  nEvents = esdTree->GetEntries();
  //nEvents = 1;
  iEventS = 0;
  iEventE = nEvents;

  // special case
  //iEventS = (iFile-1)*10;
  //iEventE = iEventS+10;
  printf("Event range: %d to %d \n",iEventS,(iEventE-1));

  // start loop over events
  //
  for (Int_t iEvent = iEventS; iEvent < iEventE; iEvent++) {

    esdTree->GetEvent(iEvent);
    
    if (useStack) {
      runLoader->GetEvent(iEvent);
      theStack = runLoader->Stack();
      nparticles = (Int_t)runLoader->TreeK()->GetEntries();
      nprimarypart = theStack->GetNprimary();
      ntracks = theStack->GetNtrack();
      printf("Event %d stack: %d %d %d \n",iEvent,nparticles,nprimarypart,ntracks);
      //continue;
    } else {
      printf("Event %d (no kine)\n",iEvent);
    }

    // needs SPD,SDD,SSD in the list of detectors
    //
    Vertex = (AliESDVertex*) esd->GetVertex();
    fXVertex = fYVertex = fZVertex = 0.0;
    if (Vertex->GetNContributors()) {
      fZVertex = Vertex->GetZ();
      fYVertex = Vertex->GetY();
      fXVertex = Vertex->GetX();      
      errXVtx = Vertex->GetXRes();
      errYVtx = Vertex->GetYRes();
      xVertex = fXVertex;
      yVertex = fYVertex;
      zVertex = fZVertex;
      useVertexFromESD = kTRUE;
    }
    //printf("ESD vertex: %10.5f %10.5f %10.5f , %10.5f %10.5f \n",fXVertex,fYVertex,fZVertex,errXVtx,errYVtx);
    //continue;
    
    nTracks = esd->GetNumberOfMuonTracks();
    
    // start loop over ESDMuon tracks
    //
    for(Int_t iTrack = 0; iTrack < nTracks; iTrack++) {

      const AliESDMuonTrack* esdTrack = esd->GetMuonTrack(iTrack);
      trackQ = esdTrack->Charge();

      // extract MUONTrack from ESDMuonTrack
      //
      muonTrack = new AliMUONTrack();
      AliMUONESDInterface::ESDToMUON(*esdTrack, *muonTrack, kFALSE);

      // has tracker?
      if (!muonTrack->GetTrackParamAtCluster()->First()) {
	continue;
      }
      // has trigger?
      if (!muonTrack->GetMatchTrigger()) {
	continue;
      }

      // all MUON tracks which trigger Apt
      //
      trackParamM = (*((AliMUONTrackParam*)(muonTrack->GetTrackParamAtCluster()->First())));
      //printf("TP0:  %8.3f   %8.3f   %8.3f \n",trackParamM.GetNonBendingCoor(),trackParamM.GetBendingCoor(),trackParamM.GetZ());

      muonTrackPxyz0[0] = trackParamM.Px();
      muonTrackPxyz0[1] = trackParamM.Py();
      muonTrackPxyz0[2] = trackParamM.Pz();
      muonTrackPt0 = TMath::Sqrt(muonTrackPxyz0[0]*muonTrackPxyz0[0]+muonTrackPxyz0[1]*muonTrackPxyz0[1]);

      hPt0->Fill(muonTrackPt0);

      esdTrackPxyzVtx[0] = esdTrack->Px();
      esdTrackPxyzVtx[1] = esdTrack->Py();
      esdTrackPxyzVtx[2] = esdTrack->Pz();
      esdTrackPtVtx = TMath::Sqrt(esdTrackPxyzVtx[0]*esdTrackPxyzVtx[0]+esdTrackPxyzVtx[1]*esdTrackPxyzVtx[1]);

      hVtxPt0->Fill(esdTrackPtVtx);

      esdTrackPxyzDca[0] = esdTrack->PxAtDCA();
      esdTrackPxyzDca[1] = esdTrack->PyAtDCA();
      esdTrackPxyzDca[2] = esdTrack->PzAtDCA();
      esdTrackPtDca = TMath::Sqrt(esdTrackPxyzDca[0]*esdTrackPxyzDca[0]+esdTrackPxyzDca[1]*esdTrackPxyzDca[1]);
	    
      hDcaPt0->Fill(esdTrackPtDca);

      hEta0->Fill(esdTrack->Eta());

      esdTrackV.SetX(esdTrackPxyzVtx[0]);
      esdTrackV.SetY(esdTrackPxyzVtx[1]);
      esdTrackV.SetZ(esdTrackPxyzVtx[2]);
      esdTrackTheVtx = esdTrackV.Theta()*TMath::RadToDeg();	      
      esdTrackPhiVtx = 180.+esdTrackV.Phi()*TMath::RadToDeg();
      hEff0->Fill(esdTrackPtVtx,esdTrackTheVtx,esdTrackPhiVtx);
	
      hChi2muon->Fill(esdTrack->GetNormalizedChi2());

      if (useStack) {
	if (esdTrack->GetLabel() >= 0 && esdTrack->GetLabel() < nparticles) {
	  particle = theStack->Particle(esdTrack->GetLabel());
	  if (TMath::Abs(particle->GetPdgCode()) == 13) {
	    hMuPt0->Fill(muonTrackPt0);
	    hMuEta0->Fill(esdTrack->Eta());
	  }
	}
      }

      delete muonTrack;

    } // end loop ESDMuon tracks
    
    // change event directory; go get MFT tracks
    //
    if (mftFile->cd(Form("Event%d",iEvent))) {
      
      // - a MUONTrack contained the MUON+MFT track in the old simulations
      // - in the new simulations this is MuonForwardTrack
      //
      muonTree = (TTree*)gDirectory->Get("MuonForwardTracks");
      muonTree->SetBranchAddress("tracks",&muonTracks);
      muonTree->GetEvent(0);
      
      mftTree = (TTree*)gDirectory->Get("MFTTracks");
      mftTree->SetBranchAddress("tracks",&mftTracks);
      mftTree->GetEvent(0);

      cellTree = (TTree*)gDirectory->Get("MFTCells");
      cellTree->SetBranchAddress("cells",&mftCells);
      cellTree->GetEvent(0);

      eventTree = (TTree*)gDirectory->Get("Events");
      eventTree->SetBranchAddress("fXVertexMC", &xVertexMC);
      eventTree->SetBranchAddress("fYVertexMC", &yVertexMC);
      eventTree->SetBranchAddress("fZVertexMC", &zVertexMC);
      eventTree->GetEvent(0);
	  
      //printf("ESD vertex MC: %10.5f %10.5f %10.5f \n",xVertexMC,yVertexMC,zVertexMC);

      if (!useVertexFromESD) {
	xVertex = xVertexMC;
	yVertex = yVertexMC;
	zVertex = zVertexMC;
      }

      nMUONTracks = muonTracks->GetEntries();
      nMFTTracks = mftTracks->GetEntries();
      nMFTCells = mftCells->GetEntries();

      // loop over all MFT tracks and fill statistics
      //
      for(Int_t iMFTTrack = 0; iMFTTrack < nMFTTracks; iMFTTrack++) {
	
	mftTrack = (AliMFTCATrack*)mftTracks->At(iMFTTrack);
	
	// statistics for all MFT tracks
	if (useStack) {
	  if (mftTrack->GetMCflag() == 1) {
	    if (mftTrack->GetMCindex() >= 0 && 
		mftTrack->GetMCindex() < nparticles) {
	      particle = theStack->Particle(mftTrack->GetMCindex());
	    }
	    phi = particle->Phi()*TMath::RadToDeg();
	    the = particle->Theta()*TMath::RadToDeg();
	    eta = particle->Eta();
	    hMftThe->Fill(the);
	    hMftPhi->Fill(phi);
	    hMftEtaPhi->Fill(phi,eta);
	    hMftThePhi->Fill(phi,the);
	  }
	} else {
	  phi = mftTrack->GetPhi();
	  the = 180.-mftTrack->GetTheta();
	  eta = -TMath::Log(TMath::Tan(the*TMath::DegToRad()/2.));
	  hMftThe->Fill(the);
	  hMftPhi->Fill(phi);
	  hMftEtaPhi->Fill(phi,eta);
	  hMftThePhi->Fill(phi,the);
	}

      }
      continue;

      //printf("Event %d tracks: ESD %5d MUON+MFT %5d MFT %5d  VertexMC:  %8.3f  %8.3f  %8.3f  %d \n",iEvent,nTracks,nMUONTracks,nMFTTracks,xVertex,yVertex,zVertex,nMFTCells);
      //continue;

      muonMCLabel = -1;
      foundMCLabelMFT = kFALSE;
      esdToMuonTrack = kFALSE;
      theXyDifMin = +9999.;

      // loop over MUON+MFT tracks
      //
      for(Int_t iMUONTrack = 0; iMUONTrack < nMUONTracks; iMUONTrack++) {

	mmTrack = (AliMuonForwardTrack*)muonTracks->At(iMUONTrack);

	// has trigger?
	if (!mmTrack->GetMatchTrigger()) {
	  continue;
	}

	//_______________________________________________________________

	// make some stat when changing the MUON track
	//
	if (mmTrack->GetMCLabel() != muonMCLabel) {

	  if (muonMCLabel >= 0) { // ... not the first time ...
	    
	    // this is to check that there is one couple MUON + MFT
	    // tracks which passed the slopes cut
	    if (esdToMuonTrack) {
	      
	      esdToMuonTrack = kFALSE;
	      /*  
	      if (kFALSE || (theXyDifMinMuonMCLabel != theXyDifMinMftMCLabel)) {
		
	        printf("Chi2 MM %8.3f Labels:  %5d   %5d   %5d   %5d   %5d   %5d   %1d\n",chi2min,mmTrack->GetMCLabel(),muonMCLabel,chi2minMuonMCLabel,chi2minMftMCLabel,pdgCodeMUON,pdgCodeMFT,foundMCLabelMFT);
		//printf("... %8.3f  %8.3f \n",mmTrackPtDca1,stackPtVtx1);
		printf("... %8.3f  %8.3f  %8.3f \n",mmTrackPxyzDca1[0],mmTrackPxyzDca1[1],mmTrackPxyzDca1[2]);
		printf("... %8.3f  %8.3f  %8.3f \n",stackPxyzVtx1[0],stackPxyzVtx1[1],stackPxyzVtx1[2]);
		printf("... %8.3f  %8.3f  %8.3f \n",esdTrackPxyzDca1[0],esdTrackPxyzDca1[1],esdTrackPxyzDca1[2]);
		printf("TheXyDifMin: %8.3f   %5d   %5d \n",theXyDifMin,theXyDifMinMuonMCLabel,theXyDifMinMftMCLabel);
		
	      }
	      */
	      if (chi2minMuonMCLabel != chi2minMftMCLabel) {
		/*
		Int_t fstm;
		particle = theStack->Particle(chi2minMuonMCLabel);
		fstm = particle->GetFirstMother();
		printf("----- mother of MUON %5d:  ",chi2minMuonMCLabel);
		while (fstm >= 0) {
		  printf("   %5d",fstm);
		  particle = theStack->Particle(fstm);
		  fstm = particle->GetFirstMother();
		}
		printf("\n");
		particle = theStack->Particle(chi2minMftMCLabel);
		fstm = particle->GetFirstMother();
		printf("----- mother of MUON %5d:  ",chi2minMftMCLabel);
		while (fstm >= 0) {
		  printf("   %5d",fstm);
		  particle = theStack->Particle(fstm);
		  fstm = particle->GetFirstMother();
		}
		printf("\n");
		*/		  
		hPt8->Fill(muonTrackPt1);
		hDcaPt8->Fill(esdTrackPtDca1);
		hVtxPt8->Fill(esdTrackPtVtx1);
	      } else { // end diff MUON and MFT MC labels for the min chi2

		// existing match in MC label and it corresponds 
		// also to the minimum chi-square
		
		hPt3->Fill(muonTrackPt1);
		hDcaPt3->Fill(esdTrackPtDca1);
		hVtxPt3->Fill(esdTrackPtVtx1);
		
		if (chi2minMuonQ == chi2minMftQ) {
		  
		  hPt5->Fill(muonTrackPt1);
		  hDcaPt5->Fill(esdTrackPtDca1);
		  hVtxPt5->Fill(esdTrackPtVtx1);
		  
		}
		
	      }
	      
	      if (theXyDifMinMuonMCLabel == theXyDifMinMftMCLabel) {
		
		// existing match in MC label and it corresponds 
		// also to the minimum difference in slopes
		
		hPt6->Fill(muonTrackPt1);
		hDcaPt6->Fill(esdTrackPtDca1);
		hVtxPt6->Fill(esdTrackPtVtx1);
		
		if (chi2minMuonMCLabel == chi2minMftMCLabel) {
		  
		  // ... and also to the minimum chi-square

		  hPt7->Fill(muonTrackPt1);
		  hDcaPt7->Fill(esdTrackPtDca1);
		  hVtxPt7->Fill(esdTrackPtVtx1);
		  
		}
	
	
	      } else {

		hPt9->Fill(muonTrackPt1);
		hDcaPt9->Fill(esdTrackPtDca1);
		hVtxPt9->Fill(esdTrackPtVtx1);
		
	      }

	      hEff1->Fill(esdTrackPtVtx1,esdTrackTheVtx1,esdTrackPhiVtx1);

	    } else { // esdToMuonTrack = kFALSE
	    
	      // due to the cut in slopes ...  
	      printf("No MUON+MFT was retained...\n");
	    
	    } // end check esdToMuonTrack set

	    chi2min = mmTrack->GetNormalizedChi2();
	    chi2minMuonMCLabel = mmTrack->GetMCLabel();
	    chi2minMftMCLabel = mmTrack->GetTrackMCId();
	    chi2minUpdate = kTRUE;

	    foundMCLabelMFT = kFALSE;

	  } else { // ... only the first time ...

	    chi2min = mmTrack->GetNormalizedChi2();
	    chi2minMuonMCLabel = mmTrack->GetMCLabel();
	    chi2minMftMCLabel = mmTrack->GetTrackMCId();

	    chi2minUpdate = kTRUE;

	  }

	  muonMCLabel = mmTrack->GetMCLabel();

	  theXyDifMin = +9999;

	} else { // same MUON track matched with the next MFT track

	  chi2minUpdate = kFALSE;

	  if (chi2min > mmTrack->GetNormalizedChi2()) {

	    chi2min = mmTrack->GetNormalizedChi2();
	    chi2minMuonMCLabel = mmTrack->GetMCLabel();
	    chi2minMftMCLabel = mmTrack->GetTrackMCId();

	    chi2minUpdate = kTRUE;

	  }

	}

	//_______________________________________________________________

	skipit = kFALSE;

	// additional chi-square cut
	//if (mmTrack->GetNormalizedChi2() > mmChi2Cut) continue;

	trackParamMM = (*((AliMUONTrackParam*)(mmTrack->GetTrackParamAtCluster()->First())));

	muonQ = (Short_t)TMath::Sign(1.,trackParamMM.GetInverseBendingMomentum());

	if (chi2minUpdate) chi2minMuonQ = muonQ;

	AliMUONTrackExtrap::ExtrapToZCov(&trackParamMM,zVertex);
	mmTrackPxyzDca[0] = trackParamMM.Px();
	mmTrackPxyzDca[1] = trackParamMM.Py();
	mmTrackPxyzDca[2] = trackParamMM.Pz();
	mmTrackXyzDca[0] = trackParamMM.GetNonBendingCoor();
	mmTrackXyzDca[1] = trackParamMM.GetBendingCoor();
	mmTrackXyzDca[2] = trackParamMM.GetZ();
	mmTrackPDca = trackParamMM.P();
	mmTrackPtDca = TMath::Sqrt(mmTrackPxyzDca[0]*mmTrackPxyzDca[0]+mmTrackPxyzDca[1]*mmTrackPxyzDca[1]);

	// among all matches, there is one, according to the MC label, 
	// which is the correct match between the MUON track and the MFT track
	//
        if (mmTrack->GetMCLabel() == mmTrack->GetTrackMCId()) {
	  foundMCLabelMFT = kTRUE;
	  matchMM_MC = kTRUE;
	  mmTrackPxyzDca1[0] = mmTrackPxyzDca[0];
	  mmTrackPxyzDca1[1] = mmTrackPxyzDca[1];
	  mmTrackPxyzDca1[2] = mmTrackPxyzDca[2];
	  mmTrackXyzDca1[0] = trackParamMM.GetNonBendingCoor();
	  mmTrackXyzDca1[1] = trackParamMM.GetBendingCoor();
	  mmTrackXyzDca1[2] = trackParamMM.GetZ();
	  mmTrackPDca1 = mmTrackPDca;
	  mmTrackPtDca1 = mmTrackPtDca;
	} else {
	  matchMM_MC = kFALSE;
	}
	/*	 
	for (Int_t i = 0; i < mmTrack->GetTrackParamAtCluster()->GetEntries(); i++) {
	  AliMUONTrackParam tp(*((AliMUONTrackParam*)(mmTrack->GetTrackParamAtCluster()->At(i))));
	  printf("TP1: %d   %8.3f   %8.3f   %8.3f \n",i,tp.GetNonBendingCoor(),tp.GetBendingCoor(),tp.GetZ());
	}
	*/
	if (prn) {
	  printf("-------------------------------------------------------\n");
	  printf("MUON+MFT Xyz:  %8.3f  %8.3f  %8.3f \n",mmTrackXyzDca[0],mmTrackXyzDca[1],mmTrackXyzDca[2]);
	  printf("MUON+MFT Pxyz:  %8.3f  %8.3f  %8.3f  %8.3f \n",mmTrackPxyzDca[0],mmTrackPxyzDca[1],mmTrackPxyzDca[2],mmTrackPDca);
	}
	
	// go get the MFT track which matches the MUON track
	//
	for(Int_t iMFTTrack = 0; iMFTTrack < nMFTTracks; iMFTTrack++) {
	  
	  mftTrack = (AliMFTCATrack*)mftTracks->At(iMFTTrack);
	  
	  if (mftTrack->GetMCindex() == mmTrack->GetTrackMCId()) {
	    
	    mftQ = (Short_t)mftTrack->GetChargeSign();
	    //printf("Q muon , mft: %d , %d \n",muonQ,mftQ);
	    
	    if (chi2minUpdate) chi2minMftQ = mftQ;

	    // estimation of the charge sign
	    if (muonQ != mftQ) {
	      matchMM_Q = kFALSE;
	    } else {
	      matchMM_Q = kTRUE;
	    }
	    if (matchMM_MC) {
	      hPt1s0->Fill(mmTrackPtDca);
	      if (!matchMM_Q) {
		nMisQMuonMft++;
	      } else {
		hPt1s1->Fill(mmTrackPtDca);
	      }
	    }
	    
	    mftTrackV.SetPtThetaPhi(1.0,mftTrack->GetTheta()*TMath::DegToRad(),mftTrack->GetPhi()*TMath::DegToRad());
	    mftTrackTheX = TMath::ATan(mftTrackV.X()/mftTrackV.Z())*TMath::RadToDeg();
	    mftTrackTheY = TMath::ATan(mftTrackV.Y()/mftTrackV.Z())*TMath::RadToDeg();

	    // loop over the cells of the MFT track
	    //
	    for (Int_t iCell = 0; iCell < mftTrack->GetNcells(); iCell++) {
	      //printf("%d %d \n",iCell,mftTrack->GetCellGID(iCell));
	      
	      AliMFTCACell *caCell = (AliMFTCACell*)mftCells->At(mftTrack->GetCellGID(iCell));;

	      // keep the last MFT track cluster (towards the absorber)
	      //
	      if (iCell == 0) {
		//printf("C: %d   %8.3f   %8.3f   %8.3f \n",iCell,caCell->GetHit2()[0],caCell->GetHit2()[1],caCell->GetHit2()[2]);
		mftTrackEnd[0] = caCell->GetHit2()[0];
		mftTrackEnd[1] = caCell->GetHit2()[1];
		mftTrackEnd[2] = caCell->GetHit2()[2];
		break;
	      }
	      //printf("C: %d   %8.3f   %8.3f   %8.3f \n",iCell,caCell->GetHit1()[0],caCell->GetHit1()[1],caCell->GetHit1()[2]);
	      
	    }
	    
	  } // end match MFT label
	  
	} // end search for the matching MFT CA track
	
	// go get the ESD MUON track which was built from the MUON track
	// before adding the MFT!
	// 
	for(Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
	  
	  const AliESDMuonTrack* esdTrack = esd->GetMuonTrack(iTrack);
	  if (!esdTrack->ContainTrackerData()) continue;
	  trackQ = esdTrack->Charge();
	  
	  if (esdTrack->GetLabel() == mmTrack->GetMCLabel()) {
	    
	    //printf("Q muon , esd: %d , %d \n",muonQ,trackQ);
	    
	    muonTrack = new AliMUONTrack();
	    AliMUONESDInterface::ESDToMUON(*esdTrack, *muonTrack, kFALSE);

	    trackParamMM0vtx = (*((AliMUONTrackParam*)(muonTrack->GetTrackParamAtCluster()->First())));
	    trackParamMM0dca = (*((AliMUONTrackParam*)(muonTrack->GetTrackParamAtCluster()->First())));

	    // at vertex with Branson correction
	    AliMUONTrackExtrap::ExtrapToVertex(&trackParamMM0vtx,xVertex,yVertex,zVertex,errXVtx,errYVtx);

	    // at DCA (without Branson correction)
	    AliMUONTrackExtrap::ExtrapToVertexWithoutBranson(&trackParamMM0dca,zVertex);

	    // go to the end of the MFT track
	    AliMUONTrackExtrap::ExtrapToZCov(&trackParamMM0vtx,mftTrackEnd[2]);
	    AliMUONTrackExtrap::ExtrapToZCov(&trackParamMM0dca,mftTrackEnd[2]);

	    // MUON track at exit of MFT (towards the absorber)

	    // Dca
	    esdTrackXyzDcaMft[0] = trackParamMM0dca.GetNonBendingCoor();
	    esdTrackXyzDcaMft[1] = trackParamMM0dca.GetBendingCoor();
	    esdTrackXyzDcaMft[2] = trackParamMM0dca.GetZ();
	    covDcaMFT = trackParamMM0dca.GetCovariances();
	    esdTrackDcaSqErrX = covDcaMFT(0,0);
	    esdTrackDcaSqErrY = covDcaMFT(2,2);
	    esdTrackDcaTheX = trackParamMM0dca.GetNonBendingSlope()*TMath::RadToDeg();
	    esdTrackDcaTheY = trackParamMM0dca.GetBendingSlope()*TMath::RadToDeg();

	    rxDcaMFT = esdTrackXyzDcaMft[0]-mftTrackEnd[0];
	    ryDcaMFT = esdTrackXyzDcaMft[1]-mftTrackEnd[1];
	    rxyDcaMFT = TMath::Sqrt(rxDcaMFT*rxDcaMFT+ryDcaMFT*ryDcaMFT);
	    errxyDcaMFT = TMath::Sqrt(esdTrackDcaSqErrX+esdTrackDcaSqErrY);
	    
	    // Vtx
	    esdTrackXyzVtxMft[0] = trackParamMM0vtx.GetNonBendingCoor();
	    esdTrackXyzVtxMft[1] = trackParamMM0vtx.GetBendingCoor();
	    esdTrackXyzVtxMft[2] = trackParamMM0vtx.GetZ();
	    covVtxMFT = trackParamMM0vtx.GetCovariances();
	    esdTrackVtxSqErrX = covVtxMFT(0,0);
	    esdTrackVtxSqErrY = covVtxMFT(2,2);
	    esdTrackVtxTheX = trackParamMM0vtx.GetNonBendingSlope()*TMath::RadToDeg();
	    esdTrackVtxTheY = trackParamMM0vtx.GetBendingSlope()*TMath::RadToDeg();

	    rxVtxMFT = esdTrackXyzVtxMft[0]-mftTrackEnd[0];
	    ryVtxMFT = esdTrackXyzVtxMft[1]-mftTrackEnd[1];
	    rxyVtxMFT = TMath::Sqrt(rxVtxMFT*rxVtxMFT+ryVtxMFT*ryVtxMFT);
	    errxyVtxMFT = TMath::Sqrt(esdTrackVtxSqErrX+esdTrackVtxSqErrY);

	    // use DCA (no constraint)
	    theXyDifX = esdTrackDcaTheX-mftTrackTheX;
	    theXyDifY = esdTrackDcaTheY-mftTrackTheY;
	    if ((TMath::Abs(theXyDifX) < mmDcaTheXcut) &&
		(TMath::Abs(theXyDifY) < mmDcaTheYcut)) {
	      matchMM_DcaTheXY = kTRUE;
	    } else {
	      matchMM_DcaTheXY = kFALSE;
	    }

	    // use vertex constraint
	    theXyDifX = esdTrackVtxTheX-mftTrackTheX;
	    theXyDifY = esdTrackVtxTheY-mftTrackTheY;
	    theXyDif = TMath::Sqrt(theXyDifX*theXyDifX+theXyDifY*theXyDifY);
	    if ((TMath::Abs(theXyDifX) < mmVtxTheXcut) &&
		(TMath::Abs(theXyDifY) < mmVtxTheYcut)) {
	      matchMM_VtxTheXY = kTRUE;
	    } else {
	      matchMM_VtxTheXY = kFALSE;
	    }
	    /*
	    if (!matchMM_VtxTheXY) {
	      skipit = kTRUE;
	      break;
	    }
	    
	    if (!matchMM_DcaTheXY) {
	      skipit = kTRUE;
	      break;
	    }
	    */
	    if (!matchMM_VtxTheXY || !matchMM_DcaTheXY) {
	      skipit = kTRUE;
	      break;
	    }
	    // if break here this, estToMuonTrack will no be set

	    if (theXyDifMin > theXyDif) {
	      theXyDifMin = theXyDif;
	      theXyDifMinMuonMCLabel = mmTrack->GetMCLabel();
	      theXyDifMinMftMCLabel = mmTrack->GetTrackMCId();
	    }
	    /*
	    if ((frxy*errxyDcaMFT) >= rxyDcaMFT) {
	      matchMM_Rxy = kTRUE; 
	    } else {
	      matchMM_Rxy = kFALSE; 
	    }

	    if (TMath::Abs(rxDcaMFT) < 6.0 && TMath::Abs(ryDcaMFT) < 6.0) {
	      matchMM_Rxy = kTRUE; 
	    } else {
	      matchMM_Rxy = kFALSE; 
	    }
	    */
	    if (matchMM_MC) {
	      //printf("ESD Xyz:       %8.3f  %8.3f  %8.3f   %8.3f   %8.3f   %8.3f   %8.3f (MFT)\n",esdTrackXyzMft[0],esdTrackXyzMft[1],esdTrackXyzMft[2],esdTrackSqErrX,esdTrackSqErrY,rxyMFT,errxyMFT);
	      hErrxyMFT->Fill(frxy*errxyDcaMFT-rxyDcaMFT);
	      hErrxyDifMFT1->Fill(rxyDcaMFT,errxyDcaMFT);
	      if ((kTRUE || matchMM_Q) && (kTRUE || matchMM_Rxy)) {
		hDxyDcaMft1->Fill(esdTrackXyzDcaMft[0]-mftTrackEnd[0],esdTrackXyzDcaMft[1]-mftTrackEnd[1]);
		hDxyVtxMft1->Fill(esdTrackXyzVtxMft[0]-mftTrackEnd[0],esdTrackXyzVtxMft[1]-mftTrackEnd[1]);
		hTheXyVtxDif1->Fill(esdTrackVtxTheX-mftTrackTheX,esdTrackVtxTheY-mftTrackTheY);
		hTheXyDcaDif1->Fill(esdTrackDcaTheX-mftTrackTheX,esdTrackDcaTheY-mftTrackTheY);
		hTheXyVtxDcaDif1->Fill(esdTrackVtxTheX-mftTrackTheX,esdTrackDcaTheX-mftTrackTheX);
	      }
	    } else {
	      hErrxyDifMFT0->Fill(rxyDcaMFT,errxyDcaMFT);
	      if ((kTRUE || matchMM_Q) && (kTRUE || matchMM_Rxy)) {
		hDxyDcaMft0->Fill(esdTrackXyzDcaMft[0]-mftTrackEnd[0],esdTrackXyzDcaMft[1]-mftTrackEnd[1]);
		hDxyVtxMft0->Fill(esdTrackXyzVtxMft[0]-mftTrackEnd[0],esdTrackXyzVtxMft[1]-mftTrackEnd[1]);
		hTheXyVtxDif0->Fill(esdTrackVtxTheX-mftTrackTheX,esdTrackVtxTheY-mftTrackTheY);
		hTheXyDcaDif0->Fill(esdTrackDcaTheX-mftTrackTheX,esdTrackDcaTheY-mftTrackTheY);
		hTheXyVtxDcaDif0->Fill(esdTrackVtxTheX-mftTrackTheX,esdTrackDcaTheX-mftTrackTheX);
	      }
	    }
	    
	    // ESDMuonTrack at vertex and DCA

	    esdTrackPxyzVtx[0] = esdTrack->Px();
	    esdTrackPxyzVtx[1] = esdTrack->Py();
	    esdTrackPxyzVtx[2] = esdTrack->Pz();
	    esdTrackPVtx = esdTrack->P();
	    esdTrackPtVtx = TMath::Sqrt(esdTrackPxyzVtx[0]*esdTrackPxyzVtx[0]+esdTrackPxyzVtx[1]*esdTrackPxyzVtx[1]);
	    
	    esdTrackPxyzDca[0] = esdTrack->PxAtDCA();
	    esdTrackPxyzDca[1] = esdTrack->PyAtDCA();
	    esdTrackPxyzDca[2] = esdTrack->PzAtDCA();
	    esdTrackPDca = esdTrack->PAtDCA();
	    esdTrackPtDca = TMath::Sqrt(esdTrackPxyzDca[0]*esdTrackPxyzDca[0]+esdTrackPxyzDca[1]*esdTrackPxyzDca[1]);
	    
	    esdTrackPxyzUnc[0] = esdTrack->PxUncorrected();
	    esdTrackPxyzUnc[1] = esdTrack->PyUncorrected();
	    esdTrackPxyzUnc[2] = esdTrack->PzUncorrected();
	    esdTrackPUnc = esdTrack->PUncorrected();
	    
	    esdTrackXyzVtx[0] = esdTrack->GetNonBendingCoor();
	    esdTrackXyzVtx[1] = esdTrack->GetBendingCoor();
	    esdTrackXyzVtx[2] = esdTrack->GetZ();
	    
	    esdTrackXyzDca[0] = esdTrack->GetNonBendingCoorAtDCA();
	    esdTrackXyzDca[1] = esdTrack->GetBendingCoorAtDCA();
	    esdTrackXyzDca[2] = zVertex;
	    
	    esdTrackXyzUnc[0] = esdTrack->GetNonBendingCoorUncorrected();
	    esdTrackXyzUnc[1] = esdTrack->GetBendingCoorUncorrected();
	    esdTrackXyzUnc[2] = esdTrack->GetZUncorrected();
	    
	    AliMUONTrackParam tp(*((AliMUONTrackParam*)(muonTrack->GetTrackParamAtCluster()->First())));
	    //printf("TP2: %8.3f   %8.3f   %8.3f \n",tp.GetNonBendingCoor(),tp.GetBendingCoor(),tp.GetZ());
	      
	    if (matchMM_MC) {

	      hEta1->Fill(esdTrack->Eta());

	      esdTrackPxyzDca1[0] = esdTrackPxyzDca[0];
	      esdTrackPxyzDca1[1] = esdTrackPxyzDca[1];
	      esdTrackPxyzDca1[2] = esdTrackPxyzDca[2];

	      hDxyVtx0->Fill(esdTrackXyzDca[0]-xVertex,esdTrackXyzDca[1]-yVertex);
	      hPtM->Fill(esdTrackPtVtx);
	    
	      if (prn) {
		printf("ESD Xyz:       %8.3f  %8.3f  %8.3f  (Vtx)\n",esdTrackXyzVtx[0],esdTrackXyzVtx[1],esdTrackXyzVtx[2]);
		printf("ESD Xyz:       %8.3f  %8.3f  %8.3f  (Dca)\n",esdTrackXyzDca[0],esdTrackXyzDca[1],esdTrackXyzDca[2]);
		printf("ESD Xyz:       %8.3f  %8.3f  %8.3f  (Unc)\n",esdTrackXyzUnc[0],esdTrackXyzUnc[1],esdTrackXyzUnc[2]);
		
		printf("ESD Pxyz:       %8.3f  %8.3f  %8.3f  %8.3f  (Vtx)\n",esdTrackPxyzVtx[0],esdTrackPxyzVtx[1],esdTrackPxyzVtx[2],esdTrackPVtx);
		printf("ESD Pxyz:       %8.3f  %8.3f  %8.3f  %8.3f  (Dca)\n",esdTrackPxyzDca[0],esdTrackPxyzDca[1],esdTrackPxyzDca[2],esdTrackPDca);
		
		printf("ESD Pxyz:       %8.3f  %8.3f  %8.3f  %8.3f  (Unc)\n",esdTrackPxyzUnc[0],esdTrackPxyzUnc[1],esdTrackPxyzUnc[2],esdTrackPUnc);
	      }
	      
	      if (muonQ != trackQ) {
		printf("Different MUON+MFT and ESD charge!\n");;
	      }

	      // MUON+MFT tracks with same MC label

	      muonTrackPxyz2[0] = tp.Px();
	      muonTrackPxyz2[1] = tp.Py();
	      muonTrackPxyz2[2] = tp.Pz();
	      muonTrackPt2 = TMath::Sqrt(muonTrackPxyz2[0]*muonTrackPxyz2[0]+muonTrackPxyz2[1]*muonTrackPxyz2[1]);
	      
	      hPt2->Fill(muonTrackPt2);
	      hDcaPt2->Fill(esdTrackPtDca);
	      hVtxPt2->Fill(esdTrackPtVtx);

	      if (matchMM_Q) {
		hPt4->Fill(muonTrackPt2);
		hDcaPt4->Fill(esdTrackPtDca);
		hVtxPt4->Fill(esdTrackPtVtx);
	      }

	    } // end matchMM_MC

	    if (!esdToMuonTrack) {

	      muonTrackPxyz1[0] = tp.Px();
	      muonTrackPxyz1[1] = tp.Py();
	      muonTrackPxyz1[2] = tp.Pz();
	      muonTrackPt1 = TMath::Sqrt(muonTrackPxyz1[0]*muonTrackPxyz1[0]+muonTrackPxyz1[1]*muonTrackPxyz1[1]);
	      esdTrackV.SetX(esdTrackPxyzVtx[0]);
	      esdTrackV.SetY(esdTrackPxyzVtx[1]);
	      esdTrackV.SetZ(esdTrackPxyzVtx[2]);
	      esdTrackTheVtx1 = esdTrackV.Theta()*TMath::RadToDeg();	      
	      esdTrackPhiVtx1 = 180.+esdTrackV.Phi()*TMath::RadToDeg();	      

	      hPt1->Fill(muonTrackPt1);
	      hDcaPt1->Fill(esdTrackPtDca);
	      hVtxPt1->Fill(esdTrackPtVtx);

	      esdTrackPtDca1 = esdTrackPtDca;
	      esdTrackPtVtx1 = esdTrackPtVtx;
	      
	      esdToMuonTrack = kTRUE;

	    }
	
	    delete muonTrack;

	  } // end MUON/ESD MC label match
	  
	} // end loop ESD tracks

	if (skipit) continue;
	
	if (matchMM_MC) {
	  if (matchMM_VtxTheXY) {
	    hChi2m1a->Fill(mmTrack->GetNormalizedChi2());
	  } else {
	    hChi2m1b->Fill(mmTrack->GetNormalizedChi2());
	  }
	  hDxyVtx1->Fill(mmTrackXyzDca[0]-xVertex,mmTrackXyzDca[1]-yVertex);
	  hPtMM->Fill(mmTrackPtDca);
	} else {
	  if (matchMM_VtxTheXY) {
	    hChi2m0a->Fill(mmTrack->GetNormalizedChi2());
	  } else {
	    hChi2m0b->Fill(mmTrack->GetNormalizedChi2());
	  }
	}

	if (useStack) {

	  if (mmTrack->GetMCLabel() >= 0) {
	    particle = theStack->Particle(mmTrack->GetMCLabel());
	    stackPxyzVtx[0] = particle->Px();
	    stackPxyzVtx[1] = particle->Py();
	    stackPxyzVtx[2] = particle->Pz();
	    stackPVtx = particle->P();
	    stackPtVtx = particle->Pt();
	    stackXyzVtx[0] = particle->Vx();	
	    stackXyzVtx[1] = particle->Vy();	
	    stackXyzVtx[2] = particle->Vz();	
	    
	    //particle = theStack->Particle(mmTrack->GetMCLabel());
	    if (chi2minMuonMCLabel >= 0) {
	      particle = theStack->Particle(chi2minMuonMCLabel);
	      pdgCodeMUON = particle->GetPdgCode();
	    } else {
	      pdgCodeMUON = -1;
	    }
	    
	    //particle = theStack->Particle(mmTrack->GetMFTMCLabel());
	    if (chi2minMftMCLabel >= 0) {
	      particle = theStack->Particle(chi2minMftMCLabel);
	      pdgCodeMFT = particle->GetPdgCode();
	    } else {
	      pdgCodeMFT = -1;
	    }
	    
	    if (matchMM_MC) {
	      
	      stackPxyzVtx1[0] = particle->Px();
	      stackPxyzVtx1[1] = particle->Py();
	      stackPxyzVtx1[2] = particle->Pz();
	      stackPVtx1 = particle->P();
	      stackPtVtx1 = particle->Pt();
	      stackXyzVtx1[0] = particle->Vx();	
	      stackXyzVtx1[1] = particle->Vy();	
	      stackXyzVtx1[2] = particle->Vz();	
	      
	      if (prn) {
		printf("Kine Pxyz:      %8.3f  %8.3f  %8.3f  %8.3f \n",stackPxyzVtx[0],stackPxyzVtx[1],stackPxyzVtx[2],stackPVtx);
	      }
	      //printf("Kine XYZ:      %10.5f  %10.5f  %10.5f  \n",stackXyzVtx1[0],stackXyzVtx1[1],stackXyzVtx1[2]);
	      
	      hDP0->Fill(esdTrackPDca-stackPVtx);
	      hDP1->Fill(mmTrackPDca-stackPVtx);
	      hDxy1->Fill(mmTrackXyzDca[0]-stackXyzVtx1[0],mmTrackXyzDca[1]-stackXyzVtx1[1]);
	      
	    } else {

	      hDxy0->Fill(mmTrackXyzDca[0]-stackXyzVtx[0],mmTrackXyzDca[1]-stackXyzVtx[1]);

	    }

	  } // end label >= 0

	} // end useStack

      } // end loop MUON+MFT tracks

    } // change event folder

    muonTracks->Clear("C");
    mftTracks->Clear("C");
    mftCells->Clear("C");
      
  } // end loop events
    
  printf("next file...\n");
  if (useStack) {
    runLoader->UnloadgAlice();
    runLoader->UnloadHeader();
    runLoader->UnloadKinematics();
    delete runLoader;
  }

  esdFile->Close();
  mftFile->Close();

  } // end batch loop

  //in.close();

  printf("MisQMuonMft: %d \n",nMisQMuonMft);
  /*
  if (hChi2m0->GetSumOfWeights() > 0.) 
    hChi2m0->Scale(1./hChi2m0->GetSumOfWeights());
  if (hChi2m1->GetSumOfWeights() > 0.) 
    hChi2m1->Scale(1./hChi2m1->GetSumOfWeights());
  if (hDxyVtx0->GetSumOfWeights() > 0.) 
    hDxyVtx0->Scale(1./hDxyVtx0->GetSumOfWeights());
  if (hDxyVtx1->GetSumOfWeights() > 0.) 
    hDxyVtx1->Scale(1./hDxyVtx1->GetSumOfWeights());
  */
  TFile *fout;

  if (onlyRef) {
    fout = new TFile(Form("readtracks05ref_%02d.root",listNr),"RECREATE"); 
    //fout = new TFile(Form("readtracks05ref_muons_%02d.root",listNr),"RECREATE"); 
    for (Int_t i = 0; i < histList1->GetSize(); i++) {
      histList1->At(i)->Write();
    }
    fout->Close();
    return;
  }

  fout = new TFile(Form("readtracks05_%02d.root",listNr),"RECREATE"); 
  for (Int_t i = 0; i < histList->GetSize(); i++) {
    histList->At(i)->Write();
  }
  fout->Close();

}
