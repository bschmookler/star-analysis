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
int mMaxEvent = 1e6;

// Event cuts
float mVtxCutMax = 80;    // Cut on primary vertex z position (cm)
bool  mVpdVtxCut = false; // Require vpd vertex to be present

// Track cuts
float  mChi2CutMin  = 0;
float  mChi2CutMax  = 65;
double mNHitCut     = 4;
double mDcaZCutMax  = 80;  // Cut on |DCA_Z| (cm); negative = disabled
double mDcaXYCutMax = 0.5; // Cut on |DCA_XY| (cm); negative = disabled

// ECal cuts. Full width of ECal towers is ~5.57 cm.
double mEcalLooseTowCut  = 15;   // Loose cut: require exactly 1 track within this distance of tower (cm)
double mEcalTightTowCut  = -1;   // Tight cut: additionally require track within this distance (cm); negative = disabled
double mEcalFrac      = 0.05; // Max energy fraction of distance-1 adjacent towers in ECal
double mEcalFrac2     = 0.05; // Max energy fraction of distance-2 adjacent towers in ECal; negative = disabled
int    mEcalColumns   = 22;  // Number of columns in ECal
int    mEcalColumnMin = 0;   // Minimum ECal column number to consider

// HCal cuts. Full width of HCal towers is ~9.99 cm.
double mHcalLooseTowCut  = 20;   // Loose cut: require exactly 1 track within this distance of tower (cm)
double mHcalTightTowCut  = -1;   // Tight cut: additionally require track within this distance (cm); negative = disabled
double mHcalFrac      = 0.15; // Max energy fraction of distance-1 adjacent towers in HCal
double mHcalFrac2     = 0.15; // Max energy fraction of distance-2 adjacent towers in HCal; negative = disabled
int    mHcalColumns   = 13;  // Number of columns in HCal
int    mHcalColumnMin = 0;   // Minimum HCal column number to consider

// Struct to cache track projections for a single event
struct TrkProjections {
  double ecal_x, ecal_y;
  double hcal_x, hcal_y;
};

// Key type for the adjacency map: {detector id, row, column}
// Implemented as a struct with operator< so it works as a std::map key in C++03
struct TowerKey {
  int det, row, col;
  TowerKey(int d, int r, int c) : det(d), row(r), col(c) {}
  bool operator<(const TowerKey& o) const {
    if (det != o.det) return det < o.det;
    if (row != o.row) return row < o.row;
    return col < o.col;
  }
};

// Get run number
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

// MIP analysis function
void mip_ana(const Char_t *inFile = "infiles.lis") {
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

  // This is a way if you want to speed up IO
  std::cout << "Explicit read status for some branches" << std::endl;
  picoReader->SetStatus("*", 0);
  picoReader->SetStatus("Event", 1);
  picoReader->SetStatus("FcsHits", 1);
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

  // Define Histograms

  // Event Primary Vertex z position
  TH1 *he1a = new TH1D("he1a", "Event Primary Vertex Z position", 100, -150, 150);
  he1a->GetXaxis()->SetTitle("Vertex Z position [cm]"); he1a->GetXaxis()->CenterTitle();
  he1a->SetLineColor(kBlue); he1a->SetLineWidth(2);

  // Event Primary Vertex z position -- vertices with Positive Rankings
  TH1 *he1b = new TH1D("he1b", "Event Primary Vertex Z position", 100, -150, 150);
  he1b->GetXaxis()->SetTitle("Vertex Z position [cm]"); he1b->GetXaxis()->CenterTitle();
  he1b->SetLineColor(kGreen); he1b->SetLineWidth(2);

  // Event Primary Vertex z position -- vertices with Positive Rankings + events w/ vpd vertex present
  TH1 *he1c = new TH1D("he1c", "Event Primary Vertex Z position", 100, -150, 150);
  he1c->GetXaxis()->SetTitle("Vertex Z position [cm]"); he1c->GetXaxis()->CenterTitle();
  he1c->SetLineColor(kMagenta); he1c->SetLineWidth(2);

  // Tower energy spectrum -- FCS ECal
  TH1 *h1a = new TH1D("h1a", "FCS ECal", 200, 0, 1.5);
  h1a->GetXaxis()->SetTitle("Tower Energy [GeV]"); h1a->GetXaxis()->CenterTitle();
  h1a->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)"); h1a->GetYaxis()->CenterTitle();
  h1a->SetLineColor(kBlue); h1a->SetLineWidth(2);

  // Tower energy spectrum Zoomed in -- FCS ECal
  TH1 *h1b = new TH1D("h1b", "FCS ECal", 500, 0.1, 0.2);
  h1b->GetXaxis()->SetTitle("Tower Energy [GeV]"); h1b->GetXaxis()->CenterTitle();
  h1b->SetLineColor(kBlue); h1b->SetLineWidth(2);

  // Tower energy spectrum with track matching -- FCS ECal
  TH1 *h1c = new TH1D("h1c", "FCS ECal", 200, 0, 1.5);
  h1c->GetXaxis()->SetTitle("Tower Energy [GeV]"); h1c->GetXaxis()->CenterTitle();
  h1c->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)"); h1c->GetYaxis()->CenterTitle();
  h1c->SetLineColor(kGreen); h1c->SetLineWidth(2);

  // Tower energy spectrum with track matching and isolation, Scenario 1 (single tower) -- FCS ECal
  TH1 *h1d_s1 = new TH1D("h1d_s1", "FCS ECal", 200, 0, 1.5);
  h1d_s1->GetXaxis()->SetTitle("Tower Energy [GeV]"); h1d_s1->GetXaxis()->CenterTitle();
  h1d_s1->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)"); h1d_s1->GetYaxis()->CenterTitle();
  h1d_s1->SetLineColor(kOrange-1); h1d_s1->SetLineWidth(2);

  // Tower energy spectrum with track matching and isolation, Scenario 2 (two towers summed) -- FCS ECal
  TH1 *h1d_s2 = new TH1D("h1d_s2", "FCS ECal", 200, 0, 1.5);
  h1d_s2->GetXaxis()->SetTitle("Tower Energy [GeV]"); h1d_s2->GetXaxis()->CenterTitle();
  h1d_s2->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)"); h1d_s2->GetYaxis()->CenterTitle();
  h1d_s2->SetLineColor(kOrange+1); h1d_s2->SetLineWidth(2);

  // Tower energy spectrum with track matching and isolation, combined Scenarios 1+2 -- FCS ECal
  TH1 *h1d = new TH1D("h1d", "FCS ECal", 200, 0, 1.5);
  h1d->GetXaxis()->SetTitle("Tower Energy [GeV]"); h1d->GetXaxis()->CenterTitle();
  h1d->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)"); h1d->GetYaxis()->CenterTitle();
  h1d->SetLineColor(kRed); h1d->SetLineWidth(2);

  // Tower energy spectrum with isolation only, Scenario 1 (single tower) -- FCS ECal
  TH1 *h1e_s1 = new TH1D("h1e_s1", "FCS ECal", 200, 0, 1.5);
  h1e_s1->GetXaxis()->SetTitle("Tower Energy [GeV]"); h1e_s1->GetXaxis()->CenterTitle();
  h1e_s1->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)"); h1e_s1->GetYaxis()->CenterTitle();
  h1e_s1->SetLineColor(kViolet-1); h1e_s1->SetLineWidth(2);

  // Tower energy spectrum with isolation only, Scenario 2 (two towers summed) -- FCS ECal
  TH1 *h1e_s2 = new TH1D("h1e_s2", "FCS ECal", 200, 0, 1.5);
  h1e_s2->GetXaxis()->SetTitle("Tower Energy [GeV]"); h1e_s2->GetXaxis()->CenterTitle();
  h1e_s2->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)"); h1e_s2->GetYaxis()->CenterTitle();
  h1e_s2->SetLineColor(kViolet+1); h1e_s2->SetLineWidth(2);

  // Tower energy spectrum with isolation only, combined Scenarios 1+2 -- FCS ECal
  TH1 *h1e = new TH1D("h1e", "FCS ECal", 200, 0, 1.5);
  h1e->GetXaxis()->SetTitle("Tower Energy [GeV]"); h1e->GetXaxis()->CenterTitle();
  h1e->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)"); h1e->GetYaxis()->CenterTitle();
  h1e->SetLineColor(kMagenta); h1e->SetLineWidth(2);

  // Tower row vs. column spectrum -- FCS ECal
  TH2 *h1f = new TH2D("h1f", "FCS ECal", 50, -25, 25, 36, -1, 35);
  h1f->GetXaxis()->SetTitle("#pm (Column Number+2)"); h1f->GetXaxis()->CenterTitle();
  h1f->GetYaxis()->SetTitle("Row Number"); h1f->GetYaxis()->CenterTitle();

  // Tower energy spectrum -- FCS HCal
  TH1 *h2a = new TH1D("h2a", "FCS HCal", 200, 0, 3);
  h2a->GetXaxis()->SetTitle("Tower Energy [GeV]"); h2a->GetXaxis()->CenterTitle();
  h2a->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)"); h2a->GetYaxis()->CenterTitle();
  h2a->SetLineColor(kBlue); h2a->SetLineWidth(2);

  // Tower energy spectrum Zoomed in -- FCS HCal
  TH1 *h2b = new TH1D("h2b", "FCS HCal", 500, 0.1, 0.2);
  h2b->GetXaxis()->SetTitle("Tower Energy [GeV]"); h2b->GetXaxis()->CenterTitle();
  h2b->SetLineColor(kBlue); h2b->SetLineWidth(2);

  // Tower energy spectrum with track matching -- FCS HCal
  TH1 *h2c = new TH1D("h2c", "FCS HCal", 200, 0, 3);
  h2c->GetXaxis()->SetTitle("Tower Energy [GeV]"); h2c->GetXaxis()->CenterTitle();
  h2c->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)"); h2c->GetYaxis()->CenterTitle();
  h2c->SetLineColor(kGreen); h2c->SetLineWidth(2);

  // Tower energy spectrum with track matching and isolation, Scenario 1 (single tower) -- FCS HCal
  TH1 *h2d_s1 = new TH1D("h2d_s1", "FCS HCal", 200, 0, 3);
  h2d_s1->GetXaxis()->SetTitle("Tower Energy [GeV]"); h2d_s1->GetXaxis()->CenterTitle();
  h2d_s1->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)"); h2d_s1->GetYaxis()->CenterTitle();
  h2d_s1->SetLineColor(kOrange-1); h2d_s1->SetLineWidth(2);

  // Tower energy spectrum with track matching and isolation, Scenario 2 (two towers summed) -- FCS HCal
  TH1 *h2d_s2 = new TH1D("h2d_s2", "FCS HCal", 200, 0, 3);
  h2d_s2->GetXaxis()->SetTitle("Tower Energy [GeV]"); h2d_s2->GetXaxis()->CenterTitle();
  h2d_s2->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)"); h2d_s2->GetYaxis()->CenterTitle();
  h2d_s2->SetLineColor(kOrange+1); h2d_s2->SetLineWidth(2);

  // Tower energy spectrum with track matching and isolation, combined Scenarios 1+2 -- FCS HCal
  TH1 *h2d = new TH1D("h2d", "FCS HCal", 200, 0, 3);
  h2d->GetXaxis()->SetTitle("Tower Energy [GeV]"); h2d->GetXaxis()->CenterTitle();
  h2d->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)"); h2d->GetYaxis()->CenterTitle();
  h2d->SetLineColor(kRed); h2d->SetLineWidth(2);

  // Tower energy spectrum with isolation only, Scenario 1 (single tower) -- FCS HCal
  TH1 *h2e_s1 = new TH1D("h2e_s1", "FCS HCal", 200, 0, 3);
  h2e_s1->GetXaxis()->SetTitle("Tower Energy [GeV]"); h2e_s1->GetXaxis()->CenterTitle();
  h2e_s1->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)"); h2e_s1->GetYaxis()->CenterTitle();
  h2e_s1->SetLineColor(kViolet-1); h2e_s1->SetLineWidth(2);

  // Tower energy spectrum with isolation only, Scenario 2 (two towers summed) -- FCS HCal
  TH1 *h2e_s2 = new TH1D("h2e_s2", "FCS HCal", 200, 0, 3);
  h2e_s2->GetXaxis()->SetTitle("Tower Energy [GeV]"); h2e_s2->GetXaxis()->CenterTitle();
  h2e_s2->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)"); h2e_s2->GetYaxis()->CenterTitle();
  h2e_s2->SetLineColor(kViolet+1); h2e_s2->SetLineWidth(2);

  // Tower energy spectrum with isolation only, combined Scenarios 1+2 -- FCS HCal
  TH1 *h2e = new TH1D("h2e", "FCS HCal", 200, 0, 3);
  h2e->GetXaxis()->SetTitle("Tower Energy [GeV]"); h2e->GetXaxis()->CenterTitle();
  h2e->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)"); h2e->GetYaxis()->CenterTitle();
  h2e->SetLineColor(kMagenta); h2e->SetLineWidth(2);

  // Tower row vs. column spectrum -- FCS HCal
  TH2 *h2f = new TH2D("h2f", "FCS HCal", 30, -15, 15, 22, -1, 21);
  h2f->GetXaxis()->SetTitle("#pm (Column Number+1)"); h2f->GetXaxis()->CenterTitle();
  h2f->GetYaxis()->SetTitle("Row Number"); h2f->GetYaxis()->CenterTitle();

  // Track-tower residuals -- FCS ECal (Scenario 1+2, isolation only)
  TH1 *h1g_dx = new TH1D("h1g_dx", "FCS ECal track-tower residual", 100, -20, 20);
  h1g_dx->GetXaxis()->SetTitle("Tower x - Projection x [cm]"); h1g_dx->GetXaxis()->CenterTitle();
  h1g_dx->GetYaxis()->SetTitle("Number of towers"); h1g_dx->GetYaxis()->CenterTitle();
  h1g_dx->SetLineColor(kBlue); h1g_dx->SetLineWidth(2);

  TH1 *h1g_dy = new TH1D("h1g_dy", "FCS ECal track-tower residual", 100, -20, 20);
  h1g_dy->GetXaxis()->SetTitle("Tower y - Projection y [cm]"); h1g_dy->GetXaxis()->CenterTitle();
  h1g_dy->GetYaxis()->SetTitle("Number of towers"); h1g_dy->GetYaxis()->CenterTitle();
  h1g_dy->SetLineColor(kBlue); h1g_dy->SetLineWidth(2);

  TH2 *h1g_dxdy = new TH2D("h1g_dxdy", "FCS ECal track-tower residual", 100, -20, 20, 100, -20, 20);
  h1g_dxdy->GetXaxis()->SetTitle("Tower x - Projection x [cm]"); h1g_dxdy->GetXaxis()->CenterTitle();
  h1g_dxdy->GetYaxis()->SetTitle("Tower y - Projection y [cm]"); h1g_dxdy->GetYaxis()->CenterTitle();

  // Track-tower residuals -- FCS HCal (Scenario 1+2, isolation only)
  TH1 *h2g_dx = new TH1D("h2g_dx", "FCS HCal track-tower residual", 100, -25, 25);
  h2g_dx->GetXaxis()->SetTitle("Tower x - Projection x [cm]"); h2g_dx->GetXaxis()->CenterTitle();
  h2g_dx->GetYaxis()->SetTitle("Number of towers"); h2g_dx->GetYaxis()->CenterTitle();
  h2g_dx->SetLineColor(kBlue); h2g_dx->SetLineWidth(2);

  TH1 *h2g_dy = new TH1D("h2g_dy", "FCS HCal track-tower residual", 100, -25, 25);
  h2g_dy->GetXaxis()->SetTitle("Tower y - Projection y [cm]"); h2g_dy->GetXaxis()->CenterTitle();
  h2g_dy->GetYaxis()->SetTitle("Number of towers"); h2g_dy->GetYaxis()->CenterTitle();
  h2g_dy->SetLineColor(kBlue); h2g_dy->SetLineWidth(2);

  TH2 *h2g_dxdy = new TH2D("h2g_dxdy", "FCS HCal track-tower residual", 100, -25, 25, 100, -25, 25);
  h2g_dxdy->GetXaxis()->SetTitle("Tower x - Projection x [cm]"); h2g_dxdy->GetXaxis()->CenterTitle();
  h2g_dxdy->GetYaxis()->SetTitle("Tower y - Projection y [cm]"); h2g_dxdy->GetYaxis()->CenterTitle();

  // Track Chi-square distribution
  TH1 *h3a = new TH1D("h3a", Form("Fwd Primary tracks w/ at least %.0f hits", mNHitCut), 200, -5, 100);
  h3a->GetXaxis()->SetTitle("Track Chi-Square"); h3a->GetXaxis()->CenterTitle();
  h3a->SetLineColor(kBlue); h3a->SetLineWidth(2);

  // Track DCA-Z distribution
  TH1 *h3b = new TH1D("h3b", Form("Fwd Primary tracks w/ at least %.0f hits", mNHitCut), 200, -110, 110);
  h3b->GetXaxis()->SetTitle("Track DCA_{Z} [cm]"); h3b->GetXaxis()->CenterTitle();
  h3b->SetLineColor(kBlue); h3b->SetLineWidth(2);

  // Track DCA-XY distribution
  TH1 *h3c = new TH1D("h3c", Form("Fwd Primary tracks w/ at least %.0f hits", mNHitCut), 200, -0.5, 8);
  h3c->GetXaxis()->SetTitle("Track DCA_{XY} [cm]"); h3c->GetXaxis()->CenterTitle();
  h3c->SetLineColor(kBlue); h3c->SetLineWidth(2);

  // Track DCA-XY vs. DCA-Z
  TH2D *h3d = new TH2D("h3d",Form("Fwd Primary tracks w/ at least %.0f hits", mNHitCut), 200, -0.5, 8, 200, -110, 110);
  h3d->GetXaxis()->SetTitle("Track DCA_{XY} [cm]"); h3d->GetXaxis()->CenterTitle();
  h3d->GetYaxis()->SetTitle("Track DCA_{Z} [cm]"); h3d->GetYaxis()->CenterTitle();

  // Loop over events
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

    // Apply event-level cuts
    TVector3 evtVtx    = event->primaryVertex();
    Float_t evtRanking = event->ranking();

    he1a->Fill(evtVtx.Z());
    if (evtRanking > 0) he1b->Fill(evtVtx.Z());
    if (evtRanking > 0 && fabs(event->vzVpd()) < 500.) he1c->Fill(evtVtx.Z());

    // Cut on primary vertex ranking
    if (evtRanking < 0) continue;

    // Cut on primary vertex z position
    if (fabs(evtVtx.z()) > mVtxCutMax) continue;

    // Cut on presence of vpd vertex
    if (mVpdVtxCut && fabs(event->vzVpd()) > 500.) continue;

    // Loop over reconstructed Fwd tracks, fill QA histograms,
    // and cache projections for good tracks
    std::vector<TrkProjections> good_trk_proj;

    int nTrks = dst->numberOfFwdTracks();
    for (int iTrk = 0; iTrk < nTrks; iTrk++) {
      StPicoFwdTrack* t = dst->fwdTrack(iTrk);

      // Remove tracks where the fit did not fully converge
      if ( !(t->didFitConvergeFully()) ) continue;

      // Remove any tracks with zero PT
      TVector3 mom = t->momentum();
      if( fabs(mom.X()) < 0.001 && fabs(mom.Y()) < 0.001 ) continue;

      if (t->isPrimaryTrack() && abs(t->numberOfFitPoints()) >= mNHitCut) {

        h3a->Fill(t->chi2());
        h3b->Fill(t->dcaZ());
        h3c->Fill(t->dcaXY());
        h3d->Fill( t->dcaXY(),t->dcaZ() );

        bool passChi2  = (t->chi2()  > mChi2CutMin && t->chi2() < mChi2CutMax);
        bool passDcaZ  = (mDcaZCutMax  < 0 || fabs(t->dcaZ())  < mDcaZCutMax);
        bool passDcaXY = (mDcaXYCutMax < 0 || fabs(t->dcaXY()) < mDcaXYCutMax);

        if (passChi2 && passDcaZ && passDcaXY) {
          TVector3 pe = t->ecalProjection();
          TVector3 ph = t->hcalProjection();
          TrkProjections proj = {pe.X(), pe.Y(), ph.X(), ph.Y()};
          good_trk_proj.push_back(proj);
        }
      }
    } //track loop

    // Build adjacency map: {det, row, col} -> energy, for O(1) neighbor lookups
    int nHits = dst->numberOfFcsHits();
    std::map<TowerKey, float> tower_energy_map;
    for (int ihit = 0; ihit < nHits; ++ihit) {
      StPicoFcsHit* h = dst->fcsHit(ihit);
      int det_id = h->detectorId();
      int nCols  = (det_id == 0 || det_id == 1) ? mEcalColumns : mHcalColumns;
      int row    = int(h->id()/nCols);
      int col    = h->id()%nCols;
      tower_energy_map[TowerKey(det_id, row, col)] = h->energy();
    }

    // Loop over FCS hits
    for (int ihit = 0; ihit < nHits; ++ihit) {
      StPicoFcsHit* h = dst->fcsHit(ihit);
      int det_id = h->detectorId();

      // Determine detector type and load detector-specific parameters
      int    nColumns, columnMin, colOffset;
      double looseTowCut, tightTowCut, frac, frac2;
      TH1   *hAll, *hZoom, *hMatched, *hMatchedIso, *hMatchedIso_s1, *hMatchedIso_s2;
      TH1   *hIso, *hIso_s1, *hIso_s2;
      TH1   *hResX, *hResY;
      TH2   *hMap, *hResDxDy;

      if (det_id == 0 || det_id == 1) {
        nColumns       = mEcalColumns;    columnMin      = mEcalColumnMin;
        looseTowCut    = mEcalLooseTowCut; tightTowCut   = mEcalTightTowCut;
        frac           = mEcalFrac;       frac2          = mEcalFrac2;
        colOffset      = 2; // ECal map x-axis uses (column+2)
        hAll           = h1a;  hZoom = h1b;  hMatched = h1c;
        hMatchedIso    = h1d;  hMatchedIso_s1 = h1d_s1;  hMatchedIso_s2 = h1d_s2;
        hIso           = h1e;  hIso_s1        = h1e_s1;  hIso_s2        = h1e_s2;
        hMap           = h1f;
        hResX          = h1g_dx;  hResY = h1g_dy;  hResDxDy = h1g_dxdy;
      } else if (det_id == 2 || det_id == 3) {
        nColumns       = mHcalColumns;    columnMin      = mHcalColumnMin;
        looseTowCut    = mHcalLooseTowCut; tightTowCut   = mHcalTightTowCut;
        frac           = mHcalFrac;       frac2          = mHcalFrac2;
        colOffset      = 1; // HCal map x-axis uses (column+1)
        hAll           = h2a;  hZoom = h2b;  hMatched = h2c;
        hMatchedIso    = h2d;  hMatchedIso_s1 = h2d_s1;  hMatchedIso_s2 = h2d_s2;
        hIso           = h2e;  hIso_s1        = h2e_s1;  hIso_s2        = h2e_s2;
        hMap           = h2f;
        hResX          = h2g_dx;  hResY = h2g_dy;  hResDxDy = h2g_dxdy;
      } else {
        continue; // Unknown detector, skip
      }

      // Get tower information
      float tower_energy = h->energy();
      int tower_row      = int(h->id()/nColumns);
      int tower_column   = h->id()%nColumns;
      StThreeVectorD xyz = fcsDb->getStarXYZ(det_id, h->id());
      double tower_x     = xyz.x();
      double tower_y     = xyz.y();

      // Some FCS tower hits have energy = 0
      if (tower_energy < 0.0001) continue;

      // Only keep columns greater than min value
      if (tower_column < columnMin) continue;

      hAll->Fill(tower_energy);
      hZoom->Fill(tower_energy);

      // ECal north (det=0) and HCal north (det=2) are flipped positive, south negative
      int flip = (det_id == 0 || det_id == 2) ? 1 : -1;
      hMap->Fill(flip*(tower_column + colOffset), tower_row);

      // Track matching: require exactly 1 good track within the loose cut distance,
      // and optionally require it is also within the tight cut distance.
      int nTracksLoose = 0;
      bool passedTight = false;
      double matched_px = 0, matched_py = 0;
      for (Size_t iProj = 0; iProj < good_trk_proj.size(); iProj++) {
        const TrkProjections& proj = good_trk_proj[iProj];
        double px = (det_id == 0 || det_id == 1) ? proj.ecal_x : proj.hcal_x;
        double py = (det_id == 0 || det_id == 1) ? proj.ecal_y : proj.hcal_y;
        if (fabs(px - tower_x) < looseTowCut && fabs(py - tower_y) < looseTowCut) {
          nTracksLoose++;
          matched_px = px;
          matched_py = py;
          if (tightTowCut < 0 ||
              (fabs(px - tower_x) < tightTowCut && fabs(py - tower_y) < tightTowCut))
            passedTight = true;
        }
      }

      // Exactly 1 track in loose window, and tight cut satisfied (or disabled)
      bool matched_to_trk = (nTracksLoose == 1) && passedTight;

      if (matched_to_trk) hMatched->Fill(tower_energy);

      // Classify tower into Scenario 1 or 2 using neighbor energy deposits.
      // Track:
      //   nEdgeAbove   -- edge-sharing (dist=1, axis-aligned) neighbors above frac
      //   partnerEnergy -- energy of that one edge neighbor if nEdgeAbove == 1
      //   nOtherAbove  -- any other dist-1 (diagonal) or dist-2 neighbors above their threshold
      int   nEdgeAbove    = 0;
      float partnerEnergy = 0;
      int   nOtherAbove   = 0;
      int   maxDist       = (frac2 >= 0) ? 2 : 1;

      for (int drow = -maxDist; drow <= maxDist; ++drow) {
        for (int dcol = -maxDist; dcol <= maxDist; ++dcol) {
          int dist = std::max(abs(drow), abs(dcol)); // Chebyshev distance
          if (dist == 0) continue; // Skip the tower itself

          std::map<TowerKey,float>::iterator it = tower_energy_map.find(
            TowerKey(det_id, tower_row+drow, tower_column+dcol));
          if (it == tower_energy_map.end()) continue; // No hit in this cell

          float adj_energy = it->second;
          bool  isEdge     = (dist == 1) && (drow == 0 || dcol == 0); // axis-aligned only

          if (dist == 1 && (adj_energy/tower_energy) > frac) {
            if (isEdge) {
              nEdgeAbove++;
              partnerEnergy = adj_energy;
            } else {
              nOtherAbove++; // Diagonal dist-1 neighbor above threshold
            }
          }
          if (dist == 2 && (adj_energy/tower_energy) > frac2) {
            nOtherAbove++; // dist-2 neighbor above threshold
          }
        }
      } //Neighbor loop

      // Scenario 1: all dist-1 and dist-2 neighbors below threshold
      bool scenario1 = (nEdgeAbove == 0 && nOtherAbove == 0);

      // Scenario 2: exactly one edge neighbor above threshold (and below tower energy
      // to avoid double-counting), all other neighbors below threshold
      bool scenario2 = (nEdgeAbove == 1 && partnerEnergy < tower_energy && nOtherAbove == 0);

      float fill_energy = 0;
      if      (scenario1) fill_energy = tower_energy;
      else if (scenario2) fill_energy = tower_energy + partnerEnergy;

      if (scenario1 || scenario2) {
        if (matched_to_trk) {
          hMatchedIso->Fill(fill_energy);
          if (scenario1) hMatchedIso_s1->Fill(fill_energy);
          if (scenario2) hMatchedIso_s2->Fill(fill_energy);
        }
        hIso->Fill(fill_energy);
        if (scenario1) hIso_s1->Fill(fill_energy);
        if (scenario2) hIso_s2->Fill(fill_energy);

        // Fill residuals for isolated towers with exactly 1 track in loose window
        if (nTracksLoose == 1) {
          hResX->Fill(tower_x - matched_px);
          hResY->Fill(tower_y - matched_py);
          hResDxDy->Fill(tower_x - matched_px, tower_y - matched_py);
        }
      }

    } //FCS hit loop

  } //event loop
  picoReader->Finish();

  // Set Style
  gStyle->SetOptStat(0);

  // Make plots
  TCanvas *b1a = new TCanvas("b1a");
  he1a->Draw();
  he1b->Draw("same");

  TLegend *legb1 = new TLegend(0.55, 0.7, 0.85, 0.875);
  legb1->SetBorderSize(0); legb1->SetTextSize(0.0275);
  legb1->SetFillStyle(0);
  legb1->AddEntry(he1a, "All events");
  legb1->AddEntry(he1b, "Events w/ Prim. Vtx Ranking > 0");
  legb1->Draw();

  TCanvas *b1b = new TCanvas("b1b");
  he1a->Draw();
  he1b->Draw("same");
  he1c->Draw("same");

  TLegend *legb2 = new TLegend(0.55, 0.7, 0.85, 0.875);
  legb2->SetBorderSize(0); legb2->SetTextSize(0.0275);
  legb2->SetFillStyle(0);
  legb2->AddEntry(he1a, "All events");
  legb2->AddEntry(he1b, "Events w/ Prim. Vtx Ranking > 0");
  legb2->AddEntry(he1c, "+Events w/ VPD Vtx present");
  legb2->Draw();

  TCanvas *b1c = new TCanvas("b1c");
  b1c->SetLogy();
  h3a->Draw();

  TCanvas *b1d = new TCanvas("b1d");
  h3b->Draw();

  TCanvas *b1e = new TCanvas("b1e");
  h3c->Draw();

  TCanvas *b1f = new TCanvas("b1f");
  h3d->Draw("colz");

  TCanvas *c1a = new TCanvas("c1a");
  c1a->SetLogy();
  h1a->Draw();

  TCanvas *c1b = new TCanvas("c1b");
  h1b->Draw();

  TCanvas *c1c = new TCanvas("c1c");
  c1c->SetLogy();
  h1a->SetMinimum(10);
  h1a->Draw();
  h1c->Draw("same");
  h1d->Draw("same");

  TLegend *leg1 = new TLegend(0.55, 0.625, 0.85, 0.875);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetFillStyle(0);
  leg1->AddEntry(h1a, "All FCS ECal tower hits");
  leg1->AddEntry(h1c, "+ good track isolation");
  leg1->AddEntry(h1d, "+ hit isolation (S1+S2)");
  leg1->Draw();

  TCanvas *c1d = new TCanvas("c1d");
  c1d->SetLogy();
  h1a->Draw();
  h1e->Draw("same");
  h1d->Draw("same");

  TLegend *leg1a = new TLegend(0.55, 0.625, 0.85, 0.875);
  leg1a->SetBorderSize(0);
  leg1a->SetTextSize(0.03);
  leg1a->SetFillStyle(0);
  leg1a->AddEntry(h1a, "All FCS ECal tower hits");
  leg1a->AddEntry(h1e, "+ hit isolation (S1+S2)");
  leg1a->AddEntry(h1d, "+ good track isolation");
  leg1a->Draw();

  // ECal diagnostic: Scenario 1 vs Scenario 2 breakdown (matched + isolation)
  TCanvas *c1c_scen = new TCanvas("c1c_scen");
  c1c_scen->SetLogy();
  h1d->Draw();
  h1d_s1->Draw("same");
  h1d_s2->Draw("same");

  TLegend *leg1b = new TLegend(0.55, 0.625, 0.85, 0.875);
  leg1b->SetBorderSize(0);
  leg1b->SetTextSize(0.03);
  leg1b->SetFillStyle(0);
  leg1b->AddEntry(h1d,    "Track + hit isolation (S1+S2)");
  leg1b->AddEntry(h1d_s1, "Scenario 1 (single tower)");
  leg1b->AddEntry(h1d_s2, "Scenario 2 (two towers summed)");
  leg1b->Draw();

  // ECal diagnostic: Scenario 1 vs Scenario 2 breakdown (isolation only)
  TCanvas *c1d_scen = new TCanvas("c1d_scen");
  c1d_scen->SetLogy();
  h1e->Draw();
  h1e_s1->Draw("same");
  h1e_s2->Draw("same");

  TLegend *leg1c = new TLegend(0.55, 0.625, 0.85, 0.875);
  leg1c->SetBorderSize(0);
  leg1c->SetTextSize(0.03);
  leg1c->SetFillStyle(0);
  leg1c->AddEntry(h1e,    "Hit Isolation (S1+S2)");
  leg1c->AddEntry(h1e_s1, "Scenario 1 (single tower)");
  leg1c->AddEntry(h1e_s2, "Scenario 2 (two towers summed)");
  leg1c->Draw();

  TCanvas *c2a = new TCanvas("c2a");
  c2a->SetLogy();
  h2a->Draw();

  TCanvas *c2b = new TCanvas("c2b");
  h2b->Draw();

  TCanvas *c2c = new TCanvas("c2c");
  c2c->SetLogy();
  h2a->SetMinimum(10);
  h2a->Draw();
  h2c->Draw("same");
  h2d->Draw("same");

  TLegend *leg2 = new TLegend(0.55, 0.625, 0.85, 0.875);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.03);
  leg2->SetFillStyle(0);
  leg2->AddEntry(h2a, "All FCS HCal tower hits");
  leg2->AddEntry(h2c, "+ good track isolation");
  leg2->AddEntry(h2d, "+ hit isolation (S1+S2)");
  leg2->Draw();

  TCanvas *c2d = new TCanvas("c2d");
  c2d->SetLogy();
  h2a->Draw();
  h2e->Draw("same");
  h2d->Draw("same");

  TLegend *leg2a = new TLegend(0.55, 0.625, 0.85, 0.875);
  leg2a->SetBorderSize(0);
  leg2a->SetTextSize(0.03);
  leg2a->SetFillStyle(0);
  leg2a->AddEntry(h2a, "All FCS HCal tower hits");
  leg2a->AddEntry(h2e, "+ hit isolation (S1+S2)");
  leg2a->AddEntry(h2d, "+ good track isolation");
  leg2a->Draw();

  // HCal diagnostic: Scenario 1 vs Scenario 2 breakdown (matched + isolation)
  TCanvas *c2c_scen = new TCanvas("c2c_scen");
  c2c_scen->SetLogy();
  h2d->Draw();
  h2d_s1->Draw("same");
  h2d_s2->Draw("same");

  TLegend *leg2b = new TLegend(0.55, 0.625, 0.85, 0.875);
  leg2b->SetBorderSize(0);
  leg2b->SetTextSize(0.03);
  leg2b->SetFillStyle(0);
  leg2b->AddEntry(h2d,    "Track + hit isolation (S1+S2)");
  leg2b->AddEntry(h2d_s1, "Scenario 1 (single tower)");
  leg2b->AddEntry(h2d_s2, "Scenario 2 (two towers summed)");
  leg2b->Draw();

  // HCal diagnostic: Scenario 1 vs Scenario 2 breakdown (isolation only)
  TCanvas *c2d_scen = new TCanvas("c2d_scen");
  c2d_scen->SetLogy();
  h2e->Draw();
  h2e_s1->Draw("same");
  h2e_s2->Draw("same");

  TLegend *leg2c = new TLegend(0.55, 0.625, 0.85, 0.875);
  leg2c->SetBorderSize(0);
  leg2c->SetTextSize(0.03);
  leg2c->SetFillStyle(0);
  leg2c->AddEntry(h2e,    "Hit Isolation (S1+S2)");
  leg2c->AddEntry(h2e_s1, "Scenario 1 (single tower)");
  leg2c->AddEntry(h2e_s2, "Scenario 2 (two towers summed)");
  leg2c->Draw();

  TCanvas *d1a = new TCanvas("d1a");
  h1f->Draw("colz");

  TCanvas *d2a = new TCanvas("d2a");
  h2f->Draw("colz");

  // ECal residuals
  TCanvas *e1a = new TCanvas("e1a");
  h1g_dx->Draw();

  TCanvas *e1b = new TCanvas("e1b");
  h1g_dy->Draw();

  TCanvas *e1c = new TCanvas("e1c");
  h1g_dxdy->Draw("colz");

  // HCal residuals
  TCanvas *e2a = new TCanvas("e2a");
  h2g_dx->Draw();

  TCanvas *e2b = new TCanvas("e2b");
  h2g_dy->Draw();

  TCanvas *e2c = new TCanvas("e2c");
  h2g_dxdy->Draw("colz");

  // Print plots to file
  b1a->Print("plots/mip_ana.pdf[");
  b1a->Print("plots/mip_ana.pdf");
  b1b->Print("plots/mip_ana.pdf");
  b1c->Print("plots/mip_ana.pdf");
  b1d->Print("plots/mip_ana.pdf");
  b1e->Print("plots/mip_ana.pdf");
  b1f->Print("plots/mip_ana.pdf");
  c1a->Print("plots/mip_ana.pdf");
  c1b->Print("plots/mip_ana.pdf");
  c1c->Print("plots/mip_ana.pdf");
  c1d->Print("plots/mip_ana.pdf");
  c1c_scen->Print("plots/mip_ana.pdf");
  c1d_scen->Print("plots/mip_ana.pdf");
  c2a->Print("plots/mip_ana.pdf");
  c2b->Print("plots/mip_ana.pdf");
  c2c->Print("plots/mip_ana.pdf");
  c2d->Print("plots/mip_ana.pdf");
  c2c_scen->Print("plots/mip_ana.pdf");
  c2d_scen->Print("plots/mip_ana.pdf");
  d1a->Print("plots/mip_ana.pdf");
  d2a->Print("plots/mip_ana.pdf");
  e1a->Print("plots/mip_ana.pdf");
  e1b->Print("plots/mip_ana.pdf");
  e1c->Print("plots/mip_ana.pdf");
  e2a->Print("plots/mip_ana.pdf");
  e2b->Print("plots/mip_ana.pdf");
  e2c->Print("plots/mip_ana.pdf");
  e2c->Print("plots/mip_ana.pdf]");

}
