#include <iostream>
#include <map>

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TText.h"
#include "TStyle.h"
#include "TLegend.h"

#include "StPicoEvent/StPicoDstReader.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoFcsHit.h"
#include "StPicoEvent/StPicoFcsCluster.h"
#include "StPicoEvent/StPicoFwdTrack.h"
#include "StFcsDbMaker/StFcsDbMaker.h"
#include "StFcsDbMaker/StFcsDb.h"

static const int mDebug = 0;
int mMaxEvent = 2e4;

// Get run number (shared utility, same as mip_ana.C)
int getrun(const char* fname) {
  TString f(fname);
  int idxst = f.Index("st_");
  int idxbr = f.Index("_", idxst+3);
  int idxar = f.Index("_", idxbr+1);
  TString r(f(idxbr+1, idxar-idxbr));
  int run = r.Atoi();
  std::cout << Form("RunNumber=%d found from %s (st=%d b=%d e=%d l=%d)",
                    run, fname, idxst, idxbr, idxar, idxar-idxbr) << std::endl;
  if (run == 0) {
    run = 23055058;
    std::cout << Form("Overwrite RunNumber = %d", run) << std::endl;
  }
  return run;
}

// FWD QA function
void fwd_qa(const Char_t *inFile = "infiles.lis") {
  int run = getrun(inFile);

  StFcsDbMaker *fcsDbMk = new StFcsDbMaker();
  fcsDbMk->Init();
  StFcsDb* fcsDb = dynamic_cast<StFcsDb*>(fcsDbMk->GetDataSet("fcsDb"));
  fcsDb->setDbAccess(0);
  fcsDb->InitRun(run);
  std::cout << "FcsDb initialized for run " << run << std::endl;
  for (int det = 0; det < 4; det++) {
    StThreeVectorD off = fcsDb->getDetectorOffset(det);
    printf("det=%d offset=%6.3f %6.3f %6.3f\n", det, off.x(), off.y(), off.z());
  }

  StPicoDstReader* picoReader = new StPicoDstReader(inFile);
  picoReader->Init();

  std::cout << "Explicit read status for some branches" << std::endl;
  picoReader->SetStatus("*", 0);
  picoReader->SetStatus("Event", 1);
  picoReader->SetStatus("FcsHits", 1);
  picoReader->SetStatus("FcsCluster", 1);
  picoReader->SetStatus("FwdTracks", 1);
  std::cout << "Status has been set" << std::endl;

  if (!picoReader->chain()) {
    std::cout << "No chain has been found." << std::endl;
  }
  Long64_t eventsInTree = picoReader->tree()->GetEntries();
  std::cout << "Events in Tree:  " << eventsInTree << std::endl;
  Long64_t events2read = picoReader->chain()->GetEntries();
  std::cout << "Events in chain: " << events2read << std::endl;
  if (events2read > mMaxEvent) {
    events2read = mMaxEvent;
    std::cout << "Limit # of event to read to " << events2read << std::endl;
  }

  // -------------------------------------------------------------------------
  // Define Histograms
  // -------------------------------------------------------------------------

  // --- Event-level ---

  // Primary vertex status
  TH1 *hVtxStatus = new TH1D("hVtxStatus", "Primary Vertex Status", 3, 0.5, 3.5);
  hVtxStatus->GetXaxis()->SetBinLabel(1, "No primary vertex");
  hVtxStatus->GetXaxis()->SetBinLabel(2, "Ranking < 0");
  hVtxStatus->GetXaxis()->SetBinLabel(3, "Ranking > 0");
  hVtxStatus->GetYaxis()->SetTitle("Number of Events");
  hVtxStatus->GetYaxis()->CenterTitle();
  hVtxStatus->SetMinimum(0);
  hVtxStatus->SetLineColor(kBlue); hVtxStatus->SetLineWidth(2);

  // VPD vertex status
  TH1 *hVpdVtxStatus = new TH1D("hVpdVtxStatus", "VPD Vertex Status", 2, 0.5, 2.5);
  hVpdVtxStatus->GetXaxis()->SetBinLabel(1, "No VPD vertex");
  hVpdVtxStatus->GetXaxis()->SetBinLabel(2, "VPD vertex present");
  hVpdVtxStatus->GetYaxis()->SetTitle("Number of Events");
  hVpdVtxStatus->GetYaxis()->CenterTitle();
  hVpdVtxStatus->SetMinimum(0);
  hVpdVtxStatus->SetLineColor(kBlue); hVpdVtxStatus->SetLineWidth(2);

  // Number of triggers per event
  TH1 *hNTriggers = new TH1D("hNTriggers", "Number of Triggers per Event", 6, -0.5, 5.5);
  hNTriggers->GetXaxis()->SetTitle("Number of Triggers");
  hNTriggers->GetXaxis()->CenterTitle();
  hNTriggers->GetYaxis()->SetTitle("Number of Events");
  hNTriggers->GetYaxis()->CenterTitle();
  hNTriggers->SetLineColor(kBlue); hNTriggers->SetLineWidth(2);

  // Small trigger IDs (< 100)
  TH1 *hTrigIdsSmall = new TH1D("hTrigIdsSmall", "Trigger IDs (small)", 20, 0, 100);
  hTrigIdsSmall->GetXaxis()->SetTitle("Trigger ID");
  hTrigIdsSmall->GetXaxis()->CenterTitle();
  hTrigIdsSmall->GetYaxis()->SetTitle("Number of Events");
  hTrigIdsSmall->GetYaxis()->CenterTitle();
  hTrigIdsSmall->SetLineColor(kBlue); hTrigIdsSmall->SetLineWidth(2);

  // Large trigger IDs (~890000), plotted as Trigger ID % 1000
  TH1 *hTrigIdsLarge = new TH1D("hTrigIdsLarge", "Trigger IDs (large)", 200, 0, 1000);
  hTrigIdsLarge->GetXaxis()->SetTitle("Trigger ID % 1000");
  hTrigIdsLarge->GetXaxis()->CenterTitle();
  hTrigIdsLarge->GetYaxis()->SetTitle("Number of Events");
  hTrigIdsLarge->GetYaxis()->CenterTitle();
  hTrigIdsLarge->SetLineColor(kBlue); hTrigIdsLarge->SetLineWidth(2);

  // --- Track-level ---

  const int nTrackTypes = 4;
  const char* trackTypeNames[nTrackTypes] = {
    "Global", "Beamline Const.", "Primary", "Fwd Vtx Const."
  };

  // Number of tracks per type
  TH1 *hTrackType = new TH1D("hTrackType", "Track Type", nTrackTypes, 0.5, nTrackTypes+0.5);
  for (int i = 0; i < nTrackTypes; i++)
    hTrackType->GetXaxis()->SetBinLabel(i+1, trackTypeNames[i]);
  hTrackType->GetYaxis()->SetTitle("Number of Tracks");
  hTrackType->GetYaxis()->CenterTitle();
  hTrackType->SetMinimum(0);
  hTrackType->SetLineColor(kBlue); hTrackType->SetLineWidth(2);

  // Number of tracks per event, one histogram per track type
  TH1 *hNTracksPerEvt[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hNTracksPerEvt[i] = new TH1D(Form("hNTracksPerEvt_%d", i), trackTypeNames[i], 100, 0, 100);
    hNTracksPerEvt[i]->GetXaxis()->SetTitle("Number of Tracks per Event");
    hNTracksPerEvt[i]->GetXaxis()->CenterTitle();
    hNTracksPerEvt[i]->GetYaxis()->SetTitle("Number of Events");
    hNTracksPerEvt[i]->GetYaxis()->CenterTitle();
    hNTracksPerEvt[i]->SetLineColor(kBlue); hNTracksPerEvt[i]->SetLineWidth(2);
  }

  // Fit convergence status, one histogram per track type
  TH1 *hFitStatus[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hFitStatus[i] = new TH1D(Form("hFitStatus_%d", i), trackTypeNames[i], 3, 0.5, 3.5);
    hFitStatus[i]->GetXaxis()->SetBinLabel(1, "No convergence");
    hFitStatus[i]->GetXaxis()->SetBinLabel(2, "Converged");
    hFitStatus[i]->GetXaxis()->SetBinLabel(3, "Fully converged");
    hFitStatus[i]->GetYaxis()->SetTitle("Number of Tracks");
    hFitStatus[i]->GetYaxis()->CenterTitle();
    hFitStatus[i]->SetMinimum(0);
    hFitStatus[i]->SetLineColor(kBlue); hFitStatus[i]->SetLineWidth(2);
  }

  // Vertex index used in track fit, one histogram per track type
  TH1 *hVtxIndex[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hVtxIndex[i] = new TH1D(Form("hVtxIndex_%d", i), trackTypeNames[i], 64, 0, 64);
    hVtxIndex[i]->GetXaxis()->SetTitle("Vertex Index");
    hVtxIndex[i]->GetXaxis()->CenterTitle();
    hVtxIndex[i]->GetYaxis()->SetTitle("Number of Tracks");
    hVtxIndex[i]->GetYaxis()->CenterTitle();
    hVtxIndex[i]->SetMinimum(0);
    hVtxIndex[i]->SetLineColor(kBlue); hVtxIndex[i]->SetLineWidth(2);
  }

  // Track charge, one histogram per track type
  TH1 *hCharge[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hCharge[i] = new TH1D(Form("hCharge_%d", i), trackTypeNames[i], 2, -1.5, 1.5);
    hCharge[i]->GetXaxis()->SetTitle("Charge");
    hCharge[i]->GetXaxis()->CenterTitle();
    hCharge[i]->GetYaxis()->SetTitle("Number of Tracks");
    hCharge[i]->GetYaxis()->CenterTitle();
    hCharge[i]->SetLineColor(kBlue); hCharge[i]->SetLineWidth(2);
  }

  // Number of fit points, one histogram per track type
  TH1 *hNFitPoints[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hNFitPoints[i] = new TH1D(Form("hNFitPoints_%d", i), trackTypeNames[i], 25, 0, 25);
    hNFitPoints[i]->GetXaxis()->SetTitle("Number of Fit Points");
    hNFitPoints[i]->GetXaxis()->CenterTitle();
    hNFitPoints[i]->GetYaxis()->SetTitle("Number of Tracks");
    hNFitPoints[i]->GetYaxis()->CenterTitle();
    hNFitPoints[i]->SetLineColor(kBlue); hNFitPoints[i]->SetLineWidth(2);
  }

  // Chi-square, one histogram per track type
  TH1 *hChi2[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hChi2[i] = new TH1D(Form("hChi2_%d", i), trackTypeNames[i], 101, -1, 100);
    hChi2[i]->GetXaxis()->SetTitle("Chi-Square");
    hChi2[i]->GetXaxis()->CenterTitle();
    hChi2[i]->GetYaxis()->SetTitle("Number of Tracks");
    hChi2[i]->GetYaxis()->CenterTitle();
    hChi2[i]->SetLineColor(kBlue); hChi2[i]->SetLineWidth(2);
  }

  // Track transverse momentum, one histogram per track type
  TH1 *hPT[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hPT[i] = new TH1D(Form("hPT_%d", i), trackTypeNames[i], 100, 0, 5);
    hPT[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hPT[i]->GetXaxis()->CenterTitle();
    hPT[i]->GetYaxis()->SetTitle("Number of Tracks");
    hPT[i]->GetYaxis()->CenterTitle();
    hPT[i]->SetLineColor(kBlue); hPT[i]->SetLineWidth(2);
  }

  // Track pseudo-rapidity, one histogram per track type
  TH1 *hEta[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hEta[i] = new TH1D(Form("hEta_%d", i), trackTypeNames[i], 100, 2.2, 4.5);
    hEta[i]->GetXaxis()->SetTitle("#eta");
    hEta[i]->GetXaxis()->CenterTitle();
    hEta[i]->GetYaxis()->SetTitle("Number of Tracks");
    hEta[i]->GetYaxis()->CenterTitle();
    hEta[i]->SetLineColor(kBlue); hEta[i]->SetLineWidth(2);
  }

  // Track phi, one histogram per track type (pT > 0 only)
  TH1 *hPhi[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hPhi[i] = new TH1D(Form("hPhi_%d", i), trackTypeNames[i], 100, -3.2, 3.2);
    hPhi[i]->GetXaxis()->SetTitle("#phi [rad]");
    hPhi[i]->GetXaxis()->CenterTitle();
    hPhi[i]->GetYaxis()->SetTitle("Number of Tracks");
    hPhi[i]->GetYaxis()->CenterTitle();
    hPhi[i]->SetLineColor(kBlue); hPhi[i]->SetLineWidth(2);
  }

  // DCA_Z (pT > 0 tracks only), one histogram per track type
  TH1 *hDcaZ[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hDcaZ[i] = new TH1D(Form("hDcaZ_%d", i), trackTypeNames[i], 110, -110, 110);
    hDcaZ[i]->GetXaxis()->SetTitle("DCA_{Z} [cm]");
    hDcaZ[i]->GetXaxis()->CenterTitle();
    hDcaZ[i]->GetYaxis()->SetTitle("Number of Tracks");
    hDcaZ[i]->GetYaxis()->CenterTitle();
    hDcaZ[i]->SetLineColor(kBlue); hDcaZ[i]->SetLineWidth(2);
  }

  // DCA_XY (pT > 0 tracks only), one histogram per track type
  TH1 *hDcaXY[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hDcaXY[i] = new TH1D(Form("hDcaXY_%d", i), trackTypeNames[i], 75, -0.5, 7);
    hDcaXY[i]->GetXaxis()->SetTitle("DCA_{XY} [cm]");
    hDcaXY[i]->GetXaxis()->CenterTitle();
    hDcaXY[i]->GetYaxis()->SetTitle("Number of Tracks");
    hDcaXY[i]->GetYaxis()->CenterTitle();
    hDcaXY[i]->SetLineColor(kBlue); hDcaXY[i]->SetLineWidth(2);
  }

  // q/pT, one histogram per track type (pT > 0 tracks only)
  TH1 *hQOverPt[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hQOverPt[i] = new TH1D(Form("hQOverPt_%d", i), trackTypeNames[i], 200, -10, 10);
    hQOverPt[i]->GetXaxis()->SetTitle("q/p_{T} [c/GeV]");
    hQOverPt[i]->GetXaxis()->CenterTitle();
    hQOverPt[i]->GetYaxis()->SetTitle("Number of Tracks");
    hQOverPt[i]->GetYaxis()->CenterTitle();
    hQOverPt[i]->SetLineColor(kBlue); hQOverPt[i]->SetLineWidth(2);
  }

  // --- Cut level 2: pT > 0 + chi2 cut (0 < chi2 < 1000) ---

  TH1 *hPhi_chi2[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hPhi_chi2[i] = new TH1D(Form("hPhi_chi2_%d", i), trackTypeNames[i], 100, -3.2, 3.2);
    hPhi_chi2[i]->GetXaxis()->SetTitle("#phi [rad]");
    hPhi_chi2[i]->GetXaxis()->CenterTitle();
    hPhi_chi2[i]->GetYaxis()->SetTitle("Number of Tracks");
    hPhi_chi2[i]->GetYaxis()->CenterTitle();
    hPhi_chi2[i]->SetLineColor(kBlue); hPhi_chi2[i]->SetLineWidth(2);
  }

  TH1 *hEta_chi2[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hEta_chi2[i] = new TH1D(Form("hEta_chi2_%d", i), trackTypeNames[i], 100, 2.2, 4.5);
    hEta_chi2[i]->GetXaxis()->SetTitle("#eta");
    hEta_chi2[i]->GetXaxis()->CenterTitle();
    hEta_chi2[i]->GetYaxis()->SetTitle("Number of Tracks");
    hEta_chi2[i]->GetYaxis()->CenterTitle();
    hEta_chi2[i]->SetLineColor(kBlue); hEta_chi2[i]->SetLineWidth(2);
  }

  TH1 *hPT_chi2[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hPT_chi2[i] = new TH1D(Form("hPT_chi2_%d", i), trackTypeNames[i], 100, 0, 5);
    hPT_chi2[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hPT_chi2[i]->GetXaxis()->CenterTitle();
    hPT_chi2[i]->GetYaxis()->SetTitle("Number of Tracks");
    hPT_chi2[i]->GetYaxis()->CenterTitle();
    hPT_chi2[i]->SetLineColor(kBlue); hPT_chi2[i]->SetLineWidth(2);
  }

  TH1 *hDcaZ_chi2[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hDcaZ_chi2[i] = new TH1D(Form("hDcaZ_chi2_%d", i), trackTypeNames[i], 110, -110, 110);
    hDcaZ_chi2[i]->GetXaxis()->SetTitle("DCA_{Z} [cm]");
    hDcaZ_chi2[i]->GetXaxis()->CenterTitle();
    hDcaZ_chi2[i]->GetYaxis()->SetTitle("Number of Tracks");
    hDcaZ_chi2[i]->GetYaxis()->CenterTitle();
    hDcaZ_chi2[i]->SetLineColor(kBlue); hDcaZ_chi2[i]->SetLineWidth(2);
  }

  TH1 *hDcaXY_chi2[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hDcaXY_chi2[i] = new TH1D(Form("hDcaXY_chi2_%d", i), trackTypeNames[i], 75, -0.5, 7);
    hDcaXY_chi2[i]->GetXaxis()->SetTitle("DCA_{XY} [cm]");
    hDcaXY_chi2[i]->GetXaxis()->CenterTitle();
    hDcaXY_chi2[i]->GetYaxis()->SetTitle("Number of Tracks");
    hDcaXY_chi2[i]->GetYaxis()->CenterTitle();
    hDcaXY_chi2[i]->SetLineColor(kBlue); hDcaXY_chi2[i]->SetLineWidth(2);
  }

  TH1 *hQOverPt_chi2[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hQOverPt_chi2[i] = new TH1D(Form("hQOverPt_chi2_%d", i), trackTypeNames[i], 200, -10, 10);
    hQOverPt_chi2[i]->GetXaxis()->SetTitle("q/p_{T} [c/GeV]");
    hQOverPt_chi2[i]->GetXaxis()->CenterTitle();
    hQOverPt_chi2[i]->GetYaxis()->SetTitle("Number of Tracks");
    hQOverPt_chi2[i]->GetYaxis()->CenterTitle();
    hQOverPt_chi2[i]->SetLineColor(kBlue); hQOverPt_chi2[i]->SetLineWidth(2);
  }

  // --- Cut level 3: pT > 0 + chi2 cut + eta cut (2.5 < eta < 4.0) ---

  TH1 *hPhi_etaCut[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hPhi_etaCut[i] = new TH1D(Form("hPhi_etaCut_%d", i), trackTypeNames[i], 100, -3.2, 3.2);
    hPhi_etaCut[i]->GetXaxis()->SetTitle("#phi [rad]");
    hPhi_etaCut[i]->GetXaxis()->CenterTitle();
    hPhi_etaCut[i]->GetYaxis()->SetTitle("Number of Tracks");
    hPhi_etaCut[i]->GetYaxis()->CenterTitle();
    hPhi_etaCut[i]->SetLineColor(kBlue); hPhi_etaCut[i]->SetLineWidth(2);
  }

  TH1 *hEta_etaCut[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hEta_etaCut[i] = new TH1D(Form("hEta_etaCut_%d", i), trackTypeNames[i], 100, 2.2, 4.5);
    hEta_etaCut[i]->GetXaxis()->SetTitle("#eta");
    hEta_etaCut[i]->GetXaxis()->CenterTitle();
    hEta_etaCut[i]->GetYaxis()->SetTitle("Number of Tracks");
    hEta_etaCut[i]->GetYaxis()->CenterTitle();
    hEta_etaCut[i]->SetLineColor(kBlue); hEta_etaCut[i]->SetLineWidth(2);
  }

  TH1 *hPT_etaCut[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hPT_etaCut[i] = new TH1D(Form("hPT_etaCut_%d", i), trackTypeNames[i], 100, 0, 5);
    hPT_etaCut[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hPT_etaCut[i]->GetXaxis()->CenterTitle();
    hPT_etaCut[i]->GetYaxis()->SetTitle("Number of Tracks");
    hPT_etaCut[i]->GetYaxis()->CenterTitle();
    hPT_etaCut[i]->SetLineColor(kBlue); hPT_etaCut[i]->SetLineWidth(2);
  }

  TH1 *hDcaZ_etaCut[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hDcaZ_etaCut[i] = new TH1D(Form("hDcaZ_etaCut_%d", i), trackTypeNames[i], 110, -110, 110);
    hDcaZ_etaCut[i]->GetXaxis()->SetTitle("DCA_{Z} [cm]");
    hDcaZ_etaCut[i]->GetXaxis()->CenterTitle();
    hDcaZ_etaCut[i]->GetYaxis()->SetTitle("Number of Tracks");
    hDcaZ_etaCut[i]->GetYaxis()->CenterTitle();
    hDcaZ_etaCut[i]->SetLineColor(kBlue); hDcaZ_etaCut[i]->SetLineWidth(2);
  }

  TH1 *hDcaXY_etaCut[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hDcaXY_etaCut[i] = new TH1D(Form("hDcaXY_etaCut_%d", i), trackTypeNames[i], 75, -0.5, 7);
    hDcaXY_etaCut[i]->GetXaxis()->SetTitle("DCA_{XY} [cm]");
    hDcaXY_etaCut[i]->GetXaxis()->CenterTitle();
    hDcaXY_etaCut[i]->GetYaxis()->SetTitle("Number of Tracks");
    hDcaXY_etaCut[i]->GetYaxis()->CenterTitle();
    hDcaXY_etaCut[i]->SetLineColor(kBlue); hDcaXY_etaCut[i]->SetLineWidth(2);
  }

  TH1 *hQOverPt_etaCut[nTrackTypes];
  for (int i = 0; i < nTrackTypes; i++) {
    hQOverPt_etaCut[i] = new TH1D(Form("hQOverPt_etaCut_%d", i), trackTypeNames[i], 200, -10, 10);
    hQOverPt_etaCut[i]->GetXaxis()->SetTitle("q/p_{T} [c/GeV]");
    hQOverPt_etaCut[i]->GetXaxis()->CenterTitle();
    hQOverPt_etaCut[i]->GetYaxis()->SetTitle("Number of Tracks");
    hQOverPt_etaCut[i]->GetYaxis()->CenterTitle();
    hQOverPt_etaCut[i]->SetLineColor(kBlue); hQOverPt_etaCut[i]->SetLineWidth(2);
  }

  // -------------------------------------------------------------------------
  // Calorimeter hit histograms
  // -------------------------------------------------------------------------

  const int nFcsDets = 4;
  const char* fcsDetNames[nFcsDets] = {
    "ECal North (det=0)", "ECal South (det=1)",
    "HCal North (det=2)", "HCal South (det=3)"
  };

  // Number of hits per event, one histogram per detector
  TH1 *hFcsNHitsPerEvt[nFcsDets];
  for (int i = 0; i < nFcsDets; i++) {
    hFcsNHitsPerEvt[i] = new TH1D(Form("hFcsNHitsPerEvt_%d", i), fcsDetNames[i], 200, 0, 200);
    hFcsNHitsPerEvt[i]->GetXaxis()->SetTitle("Number of Hits per Event");
    hFcsNHitsPerEvt[i]->GetXaxis()->CenterTitle();
    hFcsNHitsPerEvt[i]->GetYaxis()->SetTitle("Number of Events");
    hFcsNHitsPerEvt[i]->GetYaxis()->CenterTitle();
    hFcsNHitsPerEvt[i]->SetLineColor(kBlue); hFcsNHitsPerEvt[i]->SetLineWidth(2);
  }

  // Tower energy spectrum, one histogram per detector
  TH1 *hFcsEnergy[nFcsDets];
  for (int i = 0; i < nFcsDets; i++) {
    hFcsEnergy[i] = new TH1D(Form("hFcsEnergy_%d", i), fcsDetNames[i], 200, 0, 3);
    hFcsEnergy[i]->GetXaxis()->SetTitle("Tower Energy [GeV]");
    hFcsEnergy[i]->GetXaxis()->CenterTitle();
    hFcsEnergy[i]->GetYaxis()->SetTitle("Number of Towers");
    hFcsEnergy[i]->GetYaxis()->CenterTitle();
    hFcsEnergy[i]->SetLineColor(kBlue); hFcsEnergy[i]->SetLineWidth(2);
  }

  // Total energy per event, one histogram per detector
  TH1 *hFcsTotalEnergyPerEvt[nFcsDets];
  for (int i = 0; i < nFcsDets; i++) {
    hFcsTotalEnergyPerEvt[i] = new TH1D(Form("hFcsTotalEnergyPerEvt_%d", i), fcsDetNames[i], 200, 0, 50);
    hFcsTotalEnergyPerEvt[i]->GetXaxis()->SetTitle("Total Energy per Event [GeV]");
    hFcsTotalEnergyPerEvt[i]->GetXaxis()->CenterTitle();
    hFcsTotalEnergyPerEvt[i]->GetYaxis()->SetTitle("Number of Events");
    hFcsTotalEnergyPerEvt[i]->GetYaxis()->CenterTitle();
    hFcsTotalEnergyPerEvt[i]->SetLineColor(kBlue); hFcsTotalEnergyPerEvt[i]->SetLineWidth(2);
  }
  for (Long64_t iEvent = 0; iEvent < events2read; iEvent++) {
    if (iEvent%1000 == 0) std::cout << "Working on event #[" << iEvent << "/" << events2read << "]" << std::endl;

    Bool_t readEvent = picoReader->readPicoEvent(iEvent);
    if (!readEvent) {
      std::cout << "Something went wrong! Nothing to analyze..." << std::endl;
      break;
    }

    // Retrieve picoDst
    StPicoDst *dst = picoReader->picoDst();

    // Retrieve event information
    StPicoEvent *event = dst->event();
    if (!event) {
      std::cout << "Something went wrong! Event is hiding from me..." << std::endl;
      break;
    }

    // --- Event-level QA ---

    Float_t evtRanking = event->ranking();

    if      (fabs(evtRanking + 999) < 0.5) hVtxStatus->Fill(1);
    else if (evtRanking < 0)               hVtxStatus->Fill(2);
    else                                   hVtxStatus->Fill(3);

    Float_t vzVpd = event->vzVpd();
    if (fabs(vzVpd + 999) < 0.5) hVpdVtxStatus->Fill(1);
    else                          hVpdVtxStatus->Fill(2);

    std::vector<unsigned int> triggerIds = event->triggerIds();
    hNTriggers->Fill(triggerIds.size());
    for (Size_t iTrg = 0; iTrg < triggerIds.size(); iTrg++) {
      unsigned int tid = triggerIds[iTrg];
      if      (tid < 100)   hTrigIdsSmall->Fill(tid);
      else if (tid > 1000)  hTrigIdsLarge->Fill(tid % 1000);
    }

    // --- Track-level QA ---
    int nTrksByType[nTrackTypes] = {0, 0, 0, 0};
    int nTrks = dst->numberOfFwdTracks();
    for (int iTrk = 0; iTrk < nTrks; iTrk++) {
      StPicoFwdTrack* t = dst->fwdTrack(iTrk);

      int typeIdx = -1;
      if      (t->isGlobalTrack())               typeIdx = 0;
      else if (t->isBeamLineConstrainedTrack())   typeIdx = 1;
      else if (t->isPrimaryTrack())              typeIdx = 2;
      else if (t->isFwdVertexConstrainedTrack()) typeIdx = 3;

      if (typeIdx < 0) continue;

      nTrksByType[typeIdx]++;
      hTrackType->Fill(typeIdx + 1);
      hFitStatus[typeIdx]->Fill(t->status() + 1); // status 0,1,2 -> bins 1,2,3
      hVtxIndex[typeIdx]->Fill(t->vertexIndex());
      hCharge[typeIdx]->Fill(t->charge());
      hNFitPoints[typeIdx]->Fill(abs(t->numberOfFitPoints()));
      hChi2[typeIdx]->Fill(t->chi2());

      TVector3 mom = t->momentum();

      if (mom.Perp() > 0.1) {
        // Cut level 1: pT > 0
        hPT[typeIdx]->Fill(mom.Perp());
        hEta[typeIdx]->Fill(mom.Eta());
        hPhi[typeIdx]->Fill(mom.Phi());
        hDcaZ[typeIdx]->Fill(t->dcaZ());
        hDcaXY[typeIdx]->Fill(t->dcaXY());
        hQOverPt[typeIdx]->Fill(t->charge() / mom.Perp());

        // Cut level 2: pT > 0 + chi2 cut
        bool passChi2 = (t->chi2() > 0 && t->chi2() < 1000);
        if (passChi2) {
          hPhi_chi2[typeIdx]->Fill(mom.Phi());
          hEta_chi2[typeIdx]->Fill(mom.Eta());
          hPT_chi2[typeIdx]->Fill(mom.Perp());
          hDcaZ_chi2[typeIdx]->Fill(t->dcaZ());
          hDcaXY_chi2[typeIdx]->Fill(t->dcaXY());
          hQOverPt_chi2[typeIdx]->Fill(t->charge() / mom.Perp());

          // Cut level 3: pT > 0 + chi2 cut + eta cut
          bool passEta = (mom.Eta() > 2.5 && mom.Eta() < 4.0);
          if (passEta) {
            hPhi_etaCut[typeIdx]->Fill(mom.Phi());
            hEta_etaCut[typeIdx]->Fill(mom.Eta());
            hPT_etaCut[typeIdx]->Fill(mom.Perp());
            hDcaZ_etaCut[typeIdx]->Fill(t->dcaZ());
            hDcaXY_etaCut[typeIdx]->Fill(t->dcaXY());
            hQOverPt_etaCut[typeIdx]->Fill(t->charge() / mom.Perp());
          }
        }
      }
    } //track loop

    // Fill tracks-per-event histograms after the track loop
    for (int i = 0; i < nTrackTypes; i++)
      hNTracksPerEvt[i]->Fill(nTrksByType[i]);

    // --- Calorimeter hit QA ---
    int   nHitsByDet[nFcsDets]        = {0, 0, 0, 0};
    float totalEnergyByDet[nFcsDets]  = {0, 0, 0, 0};
    int nHits = dst->numberOfFcsHits();
    for (int ihit = 0; ihit < nHits; ++ihit) {
      StPicoFcsHit* h = dst->fcsHit(ihit);
      int det_id = h->detectorId();
      if (det_id < 0 || det_id >= nFcsDets) continue;

      float energy = h->energy();
      if (energy < 0.0001) continue; // Skip zero-energy hits

      nHitsByDet[det_id]++;
      totalEnergyByDet[det_id] += energy;
      hFcsEnergy[det_id]->Fill(energy);
    } //FCS hit loop

    // Fill per-event histograms after the hit loop
    for (int i = 0; i < nFcsDets; i++) {
      hFcsNHitsPerEvt[i]->Fill(nHitsByDet[i]);
      hFcsTotalEnergyPerEvt[i]->Fill(totalEnergyByDet[i]);
    }

    // --- Calorimeter cluster QA ---
    // (to be added)

    // --- Track-calorimeter matching QA ---
    // (to be added)

  } //event loop
  picoReader->Finish();

  // -------------------------------------------------------------------------
  // Plotting
  // -------------------------------------------------------------------------
  gStyle->SetOptStat(0);

  // --- Event-level ---
  TCanvas *cVtxStatus = new TCanvas("cVtxStatus");
  hVtxStatus->Draw();

  TCanvas *cVpdVtxStatus = new TCanvas("cVpdVtxStatus");
  hVpdVtxStatus->Draw();

  TCanvas *cNTriggers = new TCanvas("cNTriggers");
  hNTriggers->Draw();

  TCanvas *cTrigIdsSmall = new TCanvas("cTrigIdsSmall");
  hTrigIdsSmall->Draw();

  TCanvas *cTrigIdsLarge = new TCanvas("cTrigIdsLarge");
  hTrigIdsLarge->Draw();

  // --- Track-level ---
  TCanvas *cTrackType = new TCanvas("cTrackType");
  hTrackType->Draw();

  TCanvas *cNTracksPerEvt = new TCanvas("cNTracksPerEvt", "Number of Tracks per Event", 800, 600);
  cNTracksPerEvt->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cNTracksPerEvt->cd(i + 1);
    gPad->SetLogy();
    hNTracksPerEvt[i]->Draw();
  }

  TCanvas *cFitStatus = new TCanvas("cFitStatus", "Track Fit Convergence Status", 800, 600);
  cFitStatus->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cFitStatus->cd(i + 1);
    hFitStatus[i]->Draw();
  }

  TCanvas *cVtxIndex = new TCanvas("cVtxIndex", "Track Vertex Index", 800, 600);
  cVtxIndex->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cVtxIndex->cd(i + 1);
    hVtxIndex[i]->Draw();
  }

  TCanvas *cCharge = new TCanvas("cCharge", "Track Charge", 800, 600);
  cCharge->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cCharge->cd(i + 1);
    hCharge[i]->Draw();
  }

  TCanvas *cNFitPoints = new TCanvas("cNFitPoints", "Number of Fit Points", 800, 600);
  cNFitPoints->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cNFitPoints->cd(i + 1);
    hNFitPoints[i]->Draw();
  }

  TCanvas *cChi2 = new TCanvas("cChi2", "Track Chi-Square", 800, 600);
  cChi2->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cChi2->cd(i + 1);
    gPad->SetLogy();
    hChi2[i]->Draw();
  }

  TCanvas *cPT = new TCanvas("cPT", "Track Transverse Momentum", 800, 600);
  cPT->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cPT->cd(i + 1);
    gPad->SetLogy();
    hPT[i]->Draw();
  }

  TCanvas *cEta = new TCanvas("cEta", "Track Pseudo-Rapidity", 800, 600);
  cEta->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cEta->cd(i + 1);
    hEta[i]->Draw();
  }

  TCanvas *cPhi = new TCanvas("cPhi", "Track Phi (pT > 0)", 800, 600);
  cPhi->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cPhi->cd(i + 1);
    hPhi[i]->Draw();
  }

  TCanvas *cDcaZ = new TCanvas("cDcaZ", "Track DCA_{Z}", 800, 600);
  cDcaZ->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cDcaZ->cd(i + 1);
    hDcaZ[i]->Draw();
  }

  TCanvas *cDcaXY = new TCanvas("cDcaXY", "Track DCA_{XY}", 800, 600);
  cDcaXY->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cDcaXY->cd(i + 1);
    gPad->SetLogy();
    hDcaXY[i]->Draw();
  }

  TCanvas *cQOverPt = new TCanvas("cQOverPt", "Track q/p_{T}", 800, 600);
  cQOverPt->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cQOverPt->cd(i + 1);
    gPad->SetLogy();
    hQOverPt[i]->Draw();
  }

  // --- Cut level 2: pT > 0 + chi2 cut ---
  TCanvas *cPhi_chi2 = new TCanvas("cPhi_chi2", "Track #phi (pT>0, chi2 cut)", 800, 600);
  cPhi_chi2->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cPhi_chi2->cd(i + 1);
    hPhi_chi2[i]->Draw();
  }

  TCanvas *cEta_chi2 = new TCanvas("cEta_chi2", "Track #eta (pT>0, chi2 cut)", 800, 600);
  cEta_chi2->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cEta_chi2->cd(i + 1);
    hEta_chi2[i]->Draw();
  }

  TCanvas *cPT_chi2 = new TCanvas("cPT_chi2", "Track p_{T} (pT>0, chi2 cut)", 800, 600);
  cPT_chi2->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cPT_chi2->cd(i + 1);
    gPad->SetLogy();
    hPT_chi2[i]->Draw();
  }

  TCanvas *cDcaZ_chi2 = new TCanvas("cDcaZ_chi2", "Track DCA_{Z} (pT>0, chi2 cut)", 800, 600);
  cDcaZ_chi2->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cDcaZ_chi2->cd(i + 1);
    hDcaZ_chi2[i]->Draw();
  }

  TCanvas *cDcaXY_chi2 = new TCanvas("cDcaXY_chi2", "Track DCA_{XY} (pT>0, chi2 cut)", 800, 600);
  cDcaXY_chi2->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cDcaXY_chi2->cd(i + 1);
    gPad->SetLogy();
    hDcaXY_chi2[i]->Draw();
  }

  TCanvas *cQOverPt_chi2 = new TCanvas("cQOverPt_chi2", "Track q/p_{T} (pT>0, chi2 cut)", 800, 600);
  cQOverPt_chi2->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cQOverPt_chi2->cd(i + 1);
    gPad->SetLogy();
    hQOverPt_chi2[i]->Draw();
  }

  // --- Cut level 3: pT > 0 + chi2 cut + eta cut ---
  TCanvas *cPhi_etaCut = new TCanvas("cPhi_etaCut", "Track #phi (pT>0, chi2 cut, eta cut)", 800, 600);
  cPhi_etaCut->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cPhi_etaCut->cd(i + 1);
    hPhi_etaCut[i]->Draw();
  }

  TCanvas *cEta_etaCut = new TCanvas("cEta_etaCut", "Track #eta (pT>0, chi2 cut, eta cut)", 800, 600);
  cEta_etaCut->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cEta_etaCut->cd(i + 1);
    hEta_etaCut[i]->Draw();
  }

  TCanvas *cPT_etaCut = new TCanvas("cPT_etaCut", "Track p_{T} (pT>0, chi2 cut, eta cut)", 800, 600);
  cPT_etaCut->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cPT_etaCut->cd(i + 1);
    gPad->SetLogy();
    hPT_etaCut[i]->Draw();
  }

  TCanvas *cDcaZ_etaCut = new TCanvas("cDcaZ_etaCut", "Track DCA_{Z} (pT>0, chi2 cut, eta cut)", 800, 600);
  cDcaZ_etaCut->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cDcaZ_etaCut->cd(i + 1);
    hDcaZ_etaCut[i]->Draw();
  }

  TCanvas *cDcaXY_etaCut = new TCanvas("cDcaXY_etaCut", "Track DCA_{XY} (pT>0, chi2 cut, eta cut)", 800, 600);
  cDcaXY_etaCut->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cDcaXY_etaCut->cd(i + 1);
    gPad->SetLogy();
    hDcaXY_etaCut[i]->Draw();
  }

  TCanvas *cQOverPt_etaCut = new TCanvas("cQOverPt_etaCut", "Track q/p_{T} (pT>0, chi2 cut, eta cut)", 800, 600);
  cQOverPt_etaCut->Divide(2, 2);
  for (int i = 0; i < nTrackTypes; i++) {
    cQOverPt_etaCut->cd(i + 1);
    gPad->SetLogy();
    hQOverPt_etaCut[i]->Draw();
  }

  // --- Calorimeter hit level ---
  TCanvas *cFcsNHitsPerEvt = new TCanvas("cFcsNHitsPerEvt", "FCS Hits per Event", 800, 600);
  cFcsNHitsPerEvt->Divide(2, 2);
  for (int i = 0; i < nFcsDets; i++) {
    cFcsNHitsPerEvt->cd(i + 1);
    gPad->SetLogy();
    hFcsNHitsPerEvt[i]->Draw();
  }

  TCanvas *cFcsEnergy = new TCanvas("cFcsEnergy", "FCS Tower Energy", 800, 600);
  cFcsEnergy->Divide(2, 2);
  for (int i = 0; i < nFcsDets; i++) {
    cFcsEnergy->cd(i + 1);
    gPad->SetLogy();
    hFcsEnergy[i]->Draw();
  }

  TCanvas *cFcsTotalEnergyPerEvt = new TCanvas("cFcsTotalEnergyPerEvt", "FCS Total Energy per Event", 800, 600);
  cFcsTotalEnergyPerEvt->Divide(2, 2);
  for (int i = 0; i < nFcsDets; i++) {
    cFcsTotalEnergyPerEvt->cd(i + 1);
    gPad->SetLogy();
    hFcsTotalEnergyPerEvt[i]->Draw();
  }
  cVtxStatus->Print("plots/fwd_qa.pdf[");
  cVtxStatus->Print("plots/fwd_qa.pdf");
  cVpdVtxStatus->Print("plots/fwd_qa.pdf");
  cNTriggers->Print("plots/fwd_qa.pdf");
  cTrigIdsSmall->Print("plots/fwd_qa.pdf");
  cTrigIdsLarge->Print("plots/fwd_qa.pdf");
  cTrackType->Print("plots/fwd_qa.pdf");
  cNTracksPerEvt->Print("plots/fwd_qa.pdf");
  cFitStatus->Print("plots/fwd_qa.pdf");
  cVtxIndex->Print("plots/fwd_qa.pdf");
  cCharge->Print("plots/fwd_qa.pdf");
  cNFitPoints->Print("plots/fwd_qa.pdf");
  cChi2->Print("plots/fwd_qa.pdf");
  // phi: all three cut levels
  cPhi->Print("plots/fwd_qa.pdf");
  cPhi_chi2->Print("plots/fwd_qa.pdf");
  cPhi_etaCut->Print("plots/fwd_qa.pdf");
  // pT: all three cut levels
  cPT->Print("plots/fwd_qa.pdf");
  cPT_chi2->Print("plots/fwd_qa.pdf");
  cPT_etaCut->Print("plots/fwd_qa.pdf");
  // eta: all three cut levels
  cEta->Print("plots/fwd_qa.pdf");
  cEta_chi2->Print("plots/fwd_qa.pdf");
  cEta_etaCut->Print("plots/fwd_qa.pdf");
  // DCA_Z: all three cut levels
  cDcaZ->Print("plots/fwd_qa.pdf");
  cDcaZ_chi2->Print("plots/fwd_qa.pdf");
  cDcaZ_etaCut->Print("plots/fwd_qa.pdf");
  // DCA_XY: all three cut levels
  cDcaXY->Print("plots/fwd_qa.pdf");
  cDcaXY_chi2->Print("plots/fwd_qa.pdf");
  cDcaXY_etaCut->Print("plots/fwd_qa.pdf");
  // q/pT: all three cut levels
  cQOverPt->Print("plots/fwd_qa.pdf");
  cQOverPt_chi2->Print("plots/fwd_qa.pdf");
  cQOverPt_etaCut->Print("plots/fwd_qa.pdf");
  cFcsNHitsPerEvt->Print("plots/fwd_qa.pdf");
  cFcsEnergy->Print("plots/fwd_qa.pdf");
  cFcsTotalEnergyPerEvt->Print("plots/fwd_qa.pdf");
  cFcsTotalEnergyPerEvt->Print("plots/fwd_qa.pdf]");

}
