#include <iostream>

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

static const int mDebug=0;
int mMaxEvent=2e7;

// Event Cuts
float mVtxCutMax = 80; // Cut on primary vertex z position (cm)
bool mVpdVtxCut = false; // Require vpd vertex to be present

// Track cuts
float mChi2CutMin = 0;
float mChi2CutMax = 65;
double mNHitCut = 4;

double mEcalTowCut = 3; //Cut on track projection to Ecal (cm)
double mEcalFrac = 0.1; //Energy fraction of adjacent towers in ECal
int mEcalColumns = 22;  //Number of columns in ECal
int mEcalColumnMin = 0; //Minimum Ecal column number to consider 

double mHcalTowCut = 5; //Cut on track projection to Hcal (cm)
double mHcalFrac = 0.2; //Energy fraction of adjacent towers in HCal
int mHcalColumns = 13;  //Number of columns in HCal
int mHcalColumnMin = 0; //Minimum Hcal column number to consider

// Get run number
int getrun(const char* fname){
  TString f(fname);
  int idxst = f.Index("st_");
  int idxbr = f.Index("_",idxst+3);
  int idxar = f.Index("_",idxbr+1);
  TString r(f(idxbr+1,idxar-idxbr));
  int run=r.Atoi();
  std::cout << Form("RunNumber=%d found from %s (st=%d b=%d e=%d l=%d)",
		    run,fname,idxst,idxbr,idxar,idxar-idxbr) << std::endl;
  if(run==0){
    run=23055058;
    std::cout << Form("Overwrite RunNumber = %d",run) << std::endl; 
  }
  return run;
}

// MIP analysis function
void mip_ana(const Char_t *inFile = "infiles.lis") {
  int run=getrun(inFile);

  StFcsDbMaker *fcsDbMk = new StFcsDbMaker();
  fcsDbMk->Init();
  StFcsDb* fcsDb = dynamic_cast<StFcsDb*>(fcsDbMk->GetDataSet("fcsDb"));
  fcsDb->setDbAccess(0);
  fcsDb->InitRun(run);
  std::cout << "FcsDb initialized for run " << run << std::endl;
  for(int det=0; det<4; det++){
    StThreeVectorD off=fcsDb->getDetectorOffset(det);
    printf("det=%d offset=%6.3f %6.3f %6.3f\n",det,off.x(),off.y(),off.z());
  }

  StPicoDstReader* picoReader = new StPicoDstReader(inFile);
  picoReader->Init();

  // This is a way if you want to spead up IO
  std::cout << "Explicit read status for some branches" << std::endl;
  picoReader->SetStatus("*",0);
  picoReader->SetStatus("Event",1);  
  picoReader->SetStatus("FcsHits",1);
  picoReader->SetStatus("FwdTracks",1);
  std::cout << "Status has been set" << std::endl;

  if( !picoReader->chain() ) {
    std::cout << "No chain has been found." << std::endl;
  }
  Long64_t eventsInTree = picoReader->tree()->GetEntries();
  std::cout << "Events in Tree:  "  << eventsInTree << std::endl;
  Long64_t events2read = picoReader->chain()->GetEntries();
  std::cout << "Events in chain: " << events2read << std::endl;
  if(events2read>mMaxEvent){
    events2read=mMaxEvent;
    std::cout << "Limit # of event to read to " << events2read << std::endl;
  }
  
  // Define Histograms

  // Event Primary Vertex z position
  TH1 *he1a = new TH1D("he1a","Event Primary Vertex Z position",100,-150,150);
  he1a->GetXaxis()->SetTitle("Vertex Z position [cm]");he1a->GetXaxis()->CenterTitle();
  he1a->SetLineColor(kBlue);he1a->SetLineWidth(2);

  // Event Primary Vertex z position -- vertices with Positive Rankings
  TH1 *he1b = new TH1D("he1b","Event Primary Vertex Z position",100,-150,150);
  he1b->GetXaxis()->SetTitle("Vertex Z position [cm]");he1b->GetXaxis()->CenterTitle();
  he1b->SetLineColor(kGreen);he1b->SetLineWidth(2);

  // Event Primary Vertex z position -- vertices with Positive Rankings + events w/ vpd vertex present
  TH1 *he1c = new TH1D("he1c","Event Primary Vertex Z position",100,-150,150);
  he1c->GetXaxis()->SetTitle("Vertex Z position [cm]");he1c->GetXaxis()->CenterTitle();
  he1c->SetLineColor(kMagenta);he1c->SetLineWidth(2);

  // Tower energy spectrum -- FCS ECal
  TH1 *h1a = new TH1D("h1a","FCS ECal",200,0,3);
  h1a->GetXaxis()->SetTitle("Tower Energy [GeV]");h1a->GetXaxis()->CenterTitle();
  h1a->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)");h1a->GetYaxis()->CenterTitle();
  h1a->SetLineColor(kBlue);h1a->SetLineWidth(2);

  // Tower energy spectrum Zoomed in -- FCS ECal
  TH1 *h1b = new TH1D("h1b","FCS ECal",500,0.1,0.2);
  h1b->GetXaxis()->SetTitle("Tower Energy [GeV]");h1b->GetXaxis()->CenterTitle();
  h1b->SetLineColor(kBlue);h1b->SetLineWidth(2);
 
  // Tower energy spectrum with track matching -- FCS ECal
  TH1 *h1c = new TH1D("h1c","FCS ECal",200,0,3);
  h1c->GetXaxis()->SetTitle("Tower Energy [GeV]");h1c->GetXaxis()->CenterTitle();
  h1c->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)");h1c->GetYaxis()->CenterTitle();
  h1c->SetLineColor(kGreen);h1c->SetLineWidth(2);

  // Tower energy spectrum with track matching and isolation -- FCS ECal
  TH1 *h1d = new TH1D("h1d","FCS ECal",200,0,3);
  h1d->GetXaxis()->SetTitle("Tower Energy [GeV]");h1d->GetXaxis()->CenterTitle();
  h1d->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)");h1d->GetYaxis()->CenterTitle();
  h1d->SetLineColor(kRed);h1d->SetLineWidth(2);

  // Tower energy spectrum with isolation but no track matching -- FCS ECal
  TH1 *h1e = new TH1D("h1e","FCS ECal",200,0,3);
  h1e->GetXaxis()->SetTitle("Tower Energy [GeV]");h1e->GetXaxis()->CenterTitle();
  h1e->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)");h1e->GetYaxis()->CenterTitle();
  h1e->SetLineColor(kMagenta);h1e->SetLineWidth(2);

  // Tower row vs. column spectrum -- FCS ECal
  TH2 *h1f = new TH2D("h1f","FCS ECal",50,-25,25,36,-1,35);
  h1f->GetXaxis()->SetTitle("#pm (Column Number+2)");h1f->GetXaxis()->CenterTitle();
  h1f->GetYaxis()->SetTitle("Row Number");h1f->GetYaxis()->CenterTitle();

  // Tower energy spectrum -- FCS HCal
  TH1 *h2a = new TH1D("h2a","FCS HCal",200,0,3);
  h2a->GetXaxis()->SetTitle("Tower Energy [GeV]");h2a->GetXaxis()->CenterTitle();
  h2a->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)");h2a->GetYaxis()->CenterTitle();
  h2a->SetLineColor(kBlue);h2a->SetLineWidth(2);

  // Tower energy spectrum Zoomed in -- FCS HCal
  TH1 *h2b = new TH1D("h2b","FCS HCal",500,0.1,0.2);
  h2b->GetXaxis()->SetTitle("Tower Energy [GeV]");h2b->GetXaxis()->CenterTitle();
  h2b->SetLineColor(kBlue);h2b->SetLineWidth(2);

  // Tower energy spectrum with track matching -- FCS HCal
  TH1 *h2c = new TH1D("h2c","FCS HCal",200,0,3);
  h2c->GetXaxis()->SetTitle("Tower Energy [GeV]");h2c->GetXaxis()->CenterTitle();
  h2c->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)");h2c->GetYaxis()->CenterTitle();
  h2c->SetLineColor(kGreen);h2c->SetLineWidth(2);

  // Tower energy spectrum with track matching and isolation -- FCS HCal
  TH1 *h2d = new TH1D("h2d","FCS HCal",200,0,3);
  h2d->GetXaxis()->SetTitle("Tower Energy [GeV]");h2d->GetXaxis()->CenterTitle();
  h2d->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)");h2d->GetYaxis()->CenterTitle();
  h2d->SetLineColor(kRed);h2d->SetLineWidth(2); 
 
  // Tower energy spectrum with isolation but no track matching -- FCS HCal
  TH1 *h2e = new TH1D("h2e","FCS HCal",200,0,3);
  h2e->GetXaxis()->SetTitle("Tower Energy [GeV]");h2e->GetXaxis()->CenterTitle();
  h2e->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)");h2e->GetYaxis()->CenterTitle();
  h2e->SetLineColor(kMagenta);h2e->SetLineWidth(2);

  // Tower row vs. column spectrum -- FCS HCal
  TH2 *h2f = new TH2D("h2f","FCS HCal",30,-15,15,22,-1,21);
  h2f->GetXaxis()->SetTitle("#pm (Column Number+1)");h2f->GetXaxis()->CenterTitle();
  h2f->GetYaxis()->SetTitle("Row Number");h2f->GetYaxis()->CenterTitle();

  // Track Chi-square distrbution
  TH1 *h3a = new TH1D("h3a",Form("Fwd Primary tracks w/ at least %.0f hits",mNHitCut),200,0,100);
  h3a->GetXaxis()->SetTitle("Track Chi-Square");h3a->GetXaxis()->CenterTitle();
  h3a->SetLineColor(kBlue);h3a->SetLineWidth(2);

  // Tower energy spectrum with track matching and isolation -- FCS ECal North tower with (Column, Row) = (1,18)
  TH1 *h4e1 = new TH1D("h4e1","FCS ECal North tower: (Column, Row) = (1,18)",200,0,3);
  h4e1->GetXaxis()->SetTitle("Tower Energy [GeV]");h4e1->GetXaxis()->CenterTitle();
  h4e1->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)");h4e1->GetYaxis()->CenterTitle();
  h4e1->SetLineColor(kRed);h4e1->SetLineWidth(2);

  // Tower energy spectrum with track matching and isolation -- FCS HCal North tower with (Column, Row) = (1,12)
  TH1 *h4h1 = new TH1D("h4h1","FCS HCal North tower: (Column, Row) = (1,12)",200,0,3);
  h4h1->GetXaxis()->SetTitle("Tower Energy [GeV]");h4h1->GetXaxis()->CenterTitle();
  h4h1->GetYaxis()->SetTitle("Number of towers (/ 15 MeV)");h4h1->GetYaxis()->CenterTitle();
  h4h1->SetLineColor(kRed);h4h1->SetLineWidth(2);

  // Loop over events
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {
    if(iEvent%1000==0) std::cout << "Working on event #[" << iEvent<< "/" << events2read << "]" << std::endl;

    Bool_t readEvent = picoReader->readPicoEvent(iEvent);
    if( !readEvent ) {
      std::cout << "Something went wrong, Master! Nothing to analyze..."<< std::endl;
      break;
    }

    // Retrieve picoDst
    StPicoDst *dst = picoReader->picoDst();

    // Retrieve event information
    StPicoEvent *event = dst->event();
    if( !event ) {
      std::cout << "Something went wrong, Master! Event is hiding from me..."<< std::endl;
      break;
    }

    // Apply event-level cuts
    TVector3 evtVtx = event->primaryVertex();
    Float_t evtRanking = event->ranking();

    he1a->Fill(evtVtx.Z());
    if( evtRanking > 0 ) he1b->Fill(evtVtx.Z());
    if( evtRanking > 0 && fabs(event->vzVpd()) < 500. ) he1c->Fill(evtVtx.Z());

    // Cut on primary vertex
    if( fabs(evtVtx.z()) > mVtxCutMax && evtRanking < 0) 
	continue;

    // Cut on presence of vpd vertex
    if( mVpdVtxCut && fabs(event->vzVpd()) > 500. )
	continue;

    // Vector to hold indices of good tracks
    std::vector<int> good_trk_index;

    // Loop over reconstructed Fwd tracks and find set of 'good' tracks
    int nTrks = dst->numberOfFwdTracks();
    for(int iTrk = 0; iTrk < nTrks; iTrk++){
	StPicoFwdTrack* t = dst->fwdTrack(iTrk);

      	if( t->isPrimaryTrack() && abs(t->numberOfFitPoints()) >= mNHitCut){ // Primary tracks

	    h3a->Fill(t->chi2());

	    if( t->chi2() > mChi2CutMin && t->chi2() < mChi2CutMax ) // Chi-square cut
		good_trk_index.push_back(iTrk);
	}
    } //track loop

    //loop over FCS hits
    int nHits = dst->numberOfFcsHits();
    for (int ihit = 0; ihit < nHits; ++ihit){

	StPicoFcsHit* h = dst->fcsHit(ihit);

	//FCS Ecal
	if( h->detectorId() == 0 || h->detectorId() == 1 ){
		
           // Get tower information
           float tower_energy = h->energy();
           int tower_det = h->detectorId();
           int tower_row = int(h->id()/mEcalColumns);
           int tower_column = h->id()%mEcalColumns;
	   StThreeVectorD xyz = fcsDb->getStarXYZ(tower_det,h->id());
           auto tower_x = xyz.x();
	   auto tower_y = xyz.y();

	   // Some FCS tower hit have energy = 0.
	   if( tower_energy < 0.00001 ) continue;

           // Only keep columns greater than min value
	   if( tower_column < mEcalColumnMin ) continue;

	   h1a->Fill(tower_energy);
           h1b->Fill(tower_energy); 

	   int flip = (tower_det == 0) ? 1 : -1;
           h1f->Fill( flip*(tower_column+2) , tower_row );

	   // Require hit to be close to a track projection
	   // Full width of ECal towers ~5.57cm
	   bool matched_to_trk = false;
	   for(Size_t iTrk = 0; iTrk < good_trk_index.size(); iTrk++){

	      StPicoFwdTrack* t = dst->fwdTrack( good_trk_index[iTrk] ); 
	      TVector3 pe = t->ecalProjection();
	      double pex=pe.X();
      	      double pey=pe.Y();

	      if( fabs(pex-tower_x) < mEcalTowCut &&
	          fabs(pey-tower_y) < mEcalTowCut ){
			matched_to_trk = true;
			break;
		}
	   }

	   if(matched_to_trk) h1c->Fill(tower_energy);

	   // Require small relative energy deposit in adjacent towers
	   bool keep_tower = true;
	   for (int iAdj = 0; iAdj < nHits; ++iAdj){

	     // Don't use the same tower
	     if(iAdj == ihit) continue;

	     StPicoFcsHit* ha = dst->fcsHit(iAdj);

	     float adj_energy = ha->energy();
             int adj_det = ha->detectorId();
             int adj_row = int(ha->id()/mEcalColumns);
             int adj_column = ha->id()%mEcalColumns;

	     if( adj_det == tower_det &&
	         abs(adj_row - tower_row) <=1 &&
	         abs(adj_column - tower_column) <=1 &&
                 (adj_energy / tower_energy) > mEcalFrac ){
			keep_tower = false;
			break;
	     }
	   } //Adjacent tower loop

           if(matched_to_trk && keep_tower) h1d->Fill(tower_energy);
	   if(keep_tower) h1e->Fill(tower_energy);

	   if(tower_det == 0 && tower_column == 1 && tower_row == 18 && matched_to_trk && keep_tower)
	       h4e1->Fill(tower_energy);

	} //ECal if statement

	//FCS Hcal
	if( h->detectorId() == 2 || h->detectorId() == 3 ){
	   // Get tower information
           float tower_energy = h->energy();
           int tower_det = h->detectorId();
           int tower_row = int(h->id()/mHcalColumns);
           int tower_column = h->id()%mHcalColumns;
           StThreeVectorD xyz = fcsDb->getStarXYZ(tower_det,h->id());
           auto tower_x = xyz.x();
           auto tower_y = xyz.y();

           // Some FCS tower hit have energy = 0.
           if( tower_energy < 0.00001 ) continue;

	   // Only keep columns greater than min value
	   if( tower_column < mHcalColumnMin ) continue;

	   h2a->Fill(tower_energy);
           h2b->Fill(tower_energy);

	   int flip = (tower_det == 2) ? 1 : -1;
           h2f->Fill( flip*(tower_column+1) , tower_row );

	   // Require hit to be close to a track projection
	   // Full width of HCal towers ~9.99cm
	   bool matched_to_trk = false;
           for(Size_t iTrk = 0; iTrk < good_trk_index.size(); iTrk++){

              StPicoFwdTrack* t = dst->fwdTrack( good_trk_index[iTrk] );
              TVector3 ph = t->hcalProjection();
              double phx=ph.X();
              double phy=ph.Y();

              if( fabs(phx-tower_x) < mHcalTowCut &&
                  fabs(phy-tower_y) < mHcalTowCut ){
                        matched_to_trk = true;
                        break;
                }
           }	

	   if(matched_to_trk) h2c->Fill(tower_energy);

	   // Require small relative energy deposit in adjacent towers
           bool keep_tower = true;
           for (int iAdj = 0; iAdj < nHits; ++iAdj){

             // Don't use the same tower
             if(iAdj == ihit) continue;

             StPicoFcsHit* ha = dst->fcsHit(iAdj);

             float adj_energy = ha->energy();
             int adj_det = ha->detectorId();
             int adj_row = int(ha->id()/mHcalColumns);
             int adj_column = ha->id()%mHcalColumns;

             if( adj_det == tower_det &&
                 abs(adj_row - tower_row) <=1 &&
                 abs(adj_column - tower_column) <=1 &&
                 (adj_energy / tower_energy) > mHcalFrac )
	     {
                 keep_tower = false;
                 break;
             }
           } //Adjacent tower loop

           if(matched_to_trk && keep_tower) h2d->Fill(tower_energy);
           if(keep_tower) h2e->Fill(tower_energy);
	   if(tower_det == 2 && tower_column == 1 && tower_row == 12 && matched_to_trk && keep_tower)
               h4h1->Fill(tower_energy);

	} //HCal if statement
	
    } //FCS hit loop

  } //event loop
  picoReader->Finish();

  // Set Style
  gStyle->SetOptStat(0);

  // Make plots 
  TCanvas *b1a = new TCanvas("b1a");
  he1a->Draw();
  he1b->Draw("same");

  TLegend *legb1 = new TLegend(0.55,0.7,0.85,0.875);
  legb1->SetBorderSize(0); legb1->SetTextSize(0.0275);
  legb1->AddEntry(he1a,"All events");
  legb1->AddEntry(he1b,"Events w/ Prim. Vtx Ranking > 0");
  legb1->Draw();

  TCanvas *b1b = new TCanvas("b1b");
  he1a->Draw();
  he1b->Draw("same");
  he1c->Draw("same");

  TLegend *legb2 = new TLegend(0.55,0.7,0.85,0.875);
  legb2->SetBorderSize(0); legb2->SetTextSize(0.0275);
  legb2->AddEntry(he1a,"All events");
  legb2->AddEntry(he1b,"Events w/ Prim. Vtx Ranking > 0");
  legb2->AddEntry(he1c,"+Events w/ VPD Vtx present");
  legb2->Draw();

  TCanvas *b1c = new TCanvas("b1c");
  b1c->SetLogy();
  h3a->Draw();

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

  TLegend *leg1 = new TLegend(0.4,0.6,0.75,0.85);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->AddEntry(h1a,"All FCS ECal tower hits");
  leg1->AddEntry(h1c,"Hits matched to good track");
  leg1->AddEntry(h1d,"Hits matched to good track + isolation");
  leg1->Draw();

  TCanvas *c1d = new TCanvas("c1d");
  c1d->SetLogy();
  h1a->Draw();
  h1e->Draw("same");

  TLegend *leg1a = new TLegend(0.4,0.6,0.75,0.85);
  leg1a->SetBorderSize(0);
  leg1a->SetTextSize(0.03);
  leg1a->AddEntry(h1a,"All FCS ECal tower hits");
  leg1a->AddEntry(h1e,"Hits w/ isolation requirement only");
  leg1a->Draw();

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

  TLegend *leg2 = new TLegend(0.4,0.6,0.75,0.85);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.03);
  leg2->AddEntry(h2a,"All FCS HCal tower hits");
  leg2->AddEntry(h2c,"Hits matched to good track");
  leg2->AddEntry(h2d,"Hits matched to good track + isolation");
  leg2->Draw();

  TCanvas *c2d = new TCanvas("c2d");
  c2d->SetLogy();
  h2a->Draw();
  h2e->Draw("same");

  TLegend *leg2a = new TLegend(0.4,0.6,0.75,0.85);
  leg2a->SetBorderSize(0);
  leg2a->SetTextSize(0.03);
  leg2a->AddEntry(h1a,"All FCS HCal tower hits");
  leg2a->AddEntry(h1e,"Hits w/ isolation requirement only");
  leg2a->Draw();
 
  TCanvas *d1a = new TCanvas("d1a");
  h1f->Draw("colz");

  TCanvas *d2a = new TCanvas("d2a");
  h2f->Draw("colz");

  TCanvas *e1a = new TCanvas("e1a");
  h4e1->Draw();

  TCanvas *e2a = new TCanvas("e2a");
  h4h1->Draw();

  // Print plots to file
  b1a->Print("plots/mip_ana.pdf[");
  b1a->Print("plots/mip_ana.pdf");
  b1b->Print("plots/mip_ana.pdf");
  b1c->Print("plots/mip_ana.pdf");
  c1a->Print("plots/mip_ana.pdf");
  c1b->Print("plots/mip_ana.pdf");
  c1c->Print("plots/mip_ana.pdf");
  c1d->Print("plots/mip_ana.pdf");
  c2a->Print("plots/mip_ana.pdf");
  c2b->Print("plots/mip_ana.pdf");
  c2c->Print("plots/mip_ana.pdf");
  c2d->Print("plots/mip_ana.pdf");
  d1a->Print("plots/mip_ana.pdf");
  d2a->Print("plots/mip_ana.pdf");
  e1a->Print("plots/mip_ana.pdf");
  e2a->Print("plots/mip_ana.pdf"); 
  e2a->Print("plots/mip_ana.pdf]");

}

