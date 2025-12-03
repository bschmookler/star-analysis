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

#include "StPicoEvent/StPicoDstReader.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoFcsHit.h"
#include "StPicoEvent/StPicoFcsCluster.h"
#include "StPicoEvent/StPicoFwdTrack.h"
#include "StFcsDbMaker/StFcsDbMaker.h"
#include "StFcsDbMaker/StFcsDb.h"

static const int mDebug=0;
static const int mNCut=3;
static const int mMaxTrack=100;
double mChi2Cut=10000;
double mNHitCut=4;
double mEcalXminCut=35;
double mEcalXmaxCut=150;
double mHcalXminCut=40;
double mHcalXmaxCut=150;
double mYmaxCut=100;
double mDCut=25;           //dy cut for dx, dx cut for dy
int mMaxEvent=400000;

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

double PI=TMath::Pi();
double TWOPI=2.0*PI;
double PIO2=PI/2.0;
double THREEPIO2=PI/2.0;
double ANG(double phi){
  while(phi < -PIO2) phi+=TWOPI;
  while(phi >= THREEPIO2) phi-=TWOPI;
  return phi;
}
double DANG(double phi){
  while(phi < -PI) phi+=TWOPI;
  while(phi >= PI) phi-=TWOPI;
  return phi;
}
double R(double x, double y){return sqrt(x*x + y*y);}
double PHI(double x, double y){return ANG(atan2(y,x));}
double DPHI(double phi1, double phi2){return DANG(phi1-phi2);}

void match(const Char_t *inFile = "alive.lis") {
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
  picoReader->SetStatus("FcsClusters",1);
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
  const char *EH[2]={"Ecal","Hcal"};
  const char *NSTB[4][2]={{"North","South"},{"Top","Bottom"},{"North","South"},{"R<70","R>70"}};
  const char *CUT[4][mNCut+1]={{"Same Event","Mixed Event","Same-Mixed North","Same-Mixed South"},
			       {"Same Event","Mixed Event","Same-Mixed Top","Same-Mixed Bottom"},
			       {"Same Event","Mixed Event","Same-Mixed North","Same-Mixed South"},
			       {"Same Event","Mixed Event","Same-Mixed R<70","Same-Mixed R>70"}};
  const int  CUTCOL[mNCut+1]={kBlack,kGreen,kBlue,kRed};
  TH1F *hBunchId = new TH1F("bunchId","bunchId",120,0.0,120.0);
  TH2F *hXY[4];
  hXY[0] = new TH2F("xyEcal","Ecal; x; y",100,-150,150,100,-100,100);
  hXY[1] = new TH2F("xyHcal","Hcal; x; y",100,-150,150,100,-100,100);
  hXY[2] = new TH2F("xyETrk","Trk Ecal Projection; x; y",100,-150,150,100,-100,100);
  hXY[3] = new TH2F("xyHTrk","Trk Hcal Projection; x; y",100,-150,150,100,-100,100);
  TH1F *hE[3];
  TH1F *hET[3];
  TH1F *hdx[2][2][3], *hdy[2][2][3], *hdp[2][2][3], *hdr[2][2][3];
  TH2F *hxdx[2][3], *hydx[2][3], *hpdx[2][3], *hrdx[2][3];
  TH2F *hxdy[2][3], *hydy[2][3], *hpdy[2][3], *hrdy[2][3];
  TH2F *hxdp[2][3], *hydp[2][3], *hpdp[2][3], *hrdp[2][3];
  TH2F *hxdr[2][3], *hydr[2][3], *hpdr[2][3], *hrdr[2][3];
  for(int eh=0; eh<2; eh++){
    hE[eh]=new TH1F(Form("E%s",EH[eh]),Form("E %s; E[GeV]",EH[eh]),100,0.0,40.0);
    hET[eh]=new TH1F(Form("ET%s",EH[eh]),Form("ET %s; ET[GeV]",EH[eh]),100,0.0,5.0);
    for(int cut=0; cut<mNCut; cut++){
      hxdx[eh][cut]=new TH2F(Form("xdx%s%s",EH[eh],CUT[0][cut]),Form("dX(fcsX) %s; fcsX[cm]; dx[cm]",      EH[eh]),100,-150,150,100,-100,100);
      hydx[eh][cut]=new TH2F(Form("ydx%s%s",EH[eh],CUT[0][cut]),Form("dX(fcsY) %s; fcsY[cm]; dx[cm]",      EH[eh]),100,-100,100,100,-100,100);
      hpdx[eh][cut]=new TH2F(Form("pdx%s%s",EH[eh],CUT[0][cut]),Form("dX(fcsR*fcsPhi) %s; fcsR*Phi[cd]; dx[cm]",EH[eh]),100,-150,150,100,-100,100);
      hrdx[eh][cut]=new TH2F(Form("rdx%s%s",EH[eh],CUT[0][cut]),Form("dX(fcsR) %s; fcsR[cm]; dx[cm]",      EH[eh]),100,   0,180,100,-100,100);
      hxdy[eh][cut]=new TH2F(Form("xdy%s%s",EH[eh],CUT[1][cut]),Form("dY(fcsX) %s; fcsX[cm]; dy[cm]",      EH[eh]),100,-150,150,100,-100,100);
      hydy[eh][cut]=new TH2F(Form("ydy%s%s",EH[eh],CUT[1][cut]),Form("dY(fcsY) %s; fcsY[cm]; dy[cm]",      EH[eh]),100,-100,100,100,-100,100);
      hpdy[eh][cut]=new TH2F(Form("pdy%s%s",EH[eh],CUT[1][cut]),Form("dY(fcsR*fcsPhi) %s; fcsR*Phi[cm]; dy[cm]",EH[eh]),100,-150,150,100,-100,100);
      hrdy[eh][cut]=new TH2F(Form("rdy%s%s",EH[eh],CUT[1][cut]),Form("dY(fcsR) %s; fcsR[cm]; dy[cm]",      EH[eh]),100,   0,180,100,-100,100);
      hxdp[eh][cut]=new TH2F(Form("xdp%s%s",EH[eh],CUT[2][cut]),Form("fcsR*dPhi(fcsX) %s; fcsX[cm]; fcsr*rdphi[cm]",     EH[eh]),100,-150,150,100,-100,100);
      hydp[eh][cut]=new TH2F(Form("ydp%s%s",EH[eh],CUT[2][cut]),Form("fcsR*dPhi(fcsY) %s; fcsY[cm]; fcsr*dphi[cm]",      EH[eh]),100,-100,100,100,-100,100);
      hpdp[eh][cut]=new TH2F(Form("pdp%s%s",EH[eh],CUT[2][cut]),Form("fcsR*dPhi(fcsR*fcsPhi) %s; fcsR*Phi[cm]; fcsr*dphi[cm]",EH[eh]),100,-150,150,100,-100,100);
      hrdp[eh][cut]=new TH2F(Form("rdp%s%s",EH[eh],CUT[2][cut]),Form("fcsR*dPhi(fcsR) %s; fcsR[cm]; fcsr*dphi[cm]",      EH[eh]),100,   0,180,100,-100,100);
      hxdr[eh][cut]=new TH2F(Form("xdr%s%s",EH[eh],CUT[3][cut]),Form("dR(fcsX) %s; fcsX[cm]; dr[cm]",      EH[eh]),100,-150,150,100,-100,100);
      hydr[eh][cut]=new TH2F(Form("ydr%s%s",EH[eh],CUT[3][cut]),Form("dR(fcsY) %s; fcsY[cm]; dr[cm]",      EH[eh]),100,-100,100,100,-100,100);
      hpdr[eh][cut]=new TH2F(Form("pdr%s%s",EH[eh],CUT[3][cut]),Form("dR(fcsR*fcsPhi) %s; fcsR*Phi[cm]; dr[cm]",EH[eh]),100,-150,150,100,-100,100);
      hrdr[eh][cut]=new TH2F(Form("rdr%s%s",EH[eh],CUT[3][cut]),Form("dR(fcdR) %s; fcsR[cm]; dr[cm]",      EH[eh]),100,   0,180,100,-100,100);
      for(int nstb=0; nstb<2; nstb++){
	hdx[eh][nstb][cut] = new TH1F(Form("%s%sdx%s",EH[eh],NSTB[0][nstb],CUT[0][cut]),
				      Form("%s%s-Trk dX %s; dX[cm]",EH[eh],NSTB[0][nstb],CUT[0][cut]),100,-200.0,200.0);
	hdy[eh][nstb][cut] = new TH1F(Form("%s%s%sdy%s",EH[eh],NSTB[1][nstb],CUT[1][cut]),
				      Form("%s%s-Trk dY %s; dY[cm]",EH[eh],NSTB[1][nstb],CUT[1][cut]),100,-200.0,200.0);      
	hdp[eh][nstb][cut] = new TH1F(Form("%s%sdp%s",EH[eh],NSTB[2][nstb],CUT[2][cut]),
				      Form("%s%s-Trk R*dPhi %s; fcsR*dPhi[cm]",EH[eh],NSTB[0][nstb],CUT[2][cut]),100,-200.0,200.0);
	hdr[eh][nstb][cut] = new TH1F(Form("%s%s%sdr%s",EH[eh],NSTB[3][nstb],CUT[3][cut]),
				      Form("%s%s-Trk dR %s; dR[cm]",EH[eh],NSTB[1][nstb],CUT[3][cut]),100,-200.0,200.0);      
	hdx[eh][nstb][cut]->SetLineColor(CUTCOL[cut]);
	hdy[eh][nstb][cut]->SetLineColor(CUTCOL[cut]);
	hdp[eh][nstb][cut]->SetLineColor(CUTCOL[cut]);
	hdr[eh][nstb][cut]->SetLineColor(CUTCOL[cut]);
      }
    }
  }
  hE[2]=new TH1F(Form("ETrk"),Form("E Trk; p[GeV]"),100,0.0,40.0);
  hET[2]=new TH1F(Form("PTTrk"),Form("PT Trk; pT[GeV]"),100,0.0,5.0);

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
    hBunchId->Fill( event->bunchId() );

    //Select good tracks and keep it for this and previous event
    int same=iEvent%2;
    int ntrk[2];
    double tx[2][2][mMaxTrack],ty[2][2][mMaxTrack];    
    int nTracks = dst->numberOfFwdTracks();
    ntrk[same]=0;
    for(int it=0; it<nTracks; it++) {
      StPicoFwdTrack* t=dst->fwdTrack(it);
      TVector3 pe = t->ecalProjection();	
      TVector3 ph = t->hcalProjection();	
      double x=ph.X();
      double y=ph.Y();
      double r=sqrt(x*x + y*y);
      if(t->chi2()<mChi2Cut && 
	 t->didFitConvergeFully() && 
	 abs(t->numberOfFitPoints())>=mNHitCut &&
	 fabs(x)>mHcalXminCut &&
	 fabs(x)<mHcalXmaxCut &&
	 fabs(y)<mYmaxCut ){
	TVector3 p=t->momentum();
	hXY[2]->Fill(pe.X(),pe.Y());
	hXY[3]->Fill(ph.X(),ph.Y());
	hE[2]->Fill(p.Mag());
	hET[2]->Fill(p.Perp());
	tx[same][0][ntrk[same]]=x;
	ty[same][0][ntrk[same]]=y;
	tx[same][1][ntrk[same]]=ph.X();
	ty[same][1][ntrk[same]]=ph.Y();
	ntrk[same]++;
	if(mDebug>0) printf("ev=%6lld n=%3d trk=%3d Exy=%6.3f %6.3f Hxy=%6.3f %6.3f\n",
			    iEvent,nTracks,it,pe.X(),pe.Y(),ph.X(),ph.Y());
      }
    }

    //loop over FCS clusters
    int nClusters = dst->numberOfFcsClusters();
    for(int ic=0; ic<nClusters; ic++) {
      StPicoFcsCluster* c = dst->fcsCluster(ic);
      int det=c->detectorId();
      int ehp=fcsDb->ecalHcalPres(det);
      StThreeVectorD xyz = fcsDb->getStarXYZfromColumnRow(det,c->x(),c->y());
      double fcsx=xyz.x();
      double fcsy=xyz.y();
      double fcse=c->energy();
      double fcset=c->fourMomentum().Et();
      if(ehp==0 && (fabs(fcsx)<mEcalXminCut || fabs(fcsx)>mEcalXmaxCut)) continue;
      if(ehp==1 && (fabs(fcsx)<mHcalXminCut || fabs(fcsx)>mHcalXmaxCut)) continue;
      double fcsp=PHI(fcsx,fcsy);
      double fcsr=R(fcsx,fcsy);
      hXY[ehp]->Fill(fcsx,fcsy);
      hE[ehp]->Fill(fcse);
      hET[ehp]->Fill(fcset);
      if(mDebug>0) printf("ev=%6lld nc=%3d nt=%3d cl=%4d det=%1d ehp=%1d cr=%6.3f %6.3f xyz=%6.3f %6.3f %6.3f E=%6.3f\n",
			 iEvent,nClusters,nTracks,ic,det,ehp,c->x(),c->y(),xyz.x(),xyz.y(),xyz.z(),c->energy());
      for(int i=0; i<2; i++){ //same event or mixed event
	int evt=(same+i)%2;
	if(iEvent==0 && i==1) continue; //no mixed event for 1st event
	for(int it=0; it<ntrk[evt]; it++) {
	  if(fcsx * tx[evt][ehp][it] > 0){ //same side
	    //if(ehp==0 && fcse>0.5) continue;  //select MIP for ecal-trk
	    //if(ehp==1 && fcse<1.0) continue; //cut off low evenry for hcal-trk
	    double trkx=tx[evt][ehp][it];
	    double trky=ty[evt][ehp][it];
	    double trkp=PHI(trkx,trky);
	    double trkr=R(trkx,trky);
	    double fcsrp=fcsr*fcsp;
	    double dx=fcsx - trkx;
	    double dy=fcsy - trky;
	    double rdp=fcsr*DPHI(fcsp,trkp);
	    double dr=fcsr - trkr;
	    int nstb=0;
	    if(fcsx>0) {nstb=0;} else {nstb=1;}
	    if(fabs(dy)<mDCut) {
	      hdx[ehp][nstb][i]->Fill(dx);  
	      hxdx[ehp][i]->Fill(fcsx,dx); 
	      hydx[ehp][i]->Fill(fcsy,dx);
	      hpdx[ehp][i]->Fill(fcsrp,dx);
	      hrdx[ehp][i]->Fill(fcsr,dx);
	    }
	    if(fcsy>0) {nstb=0;} else {nstb=1;}
	    if(fabs(dx)<mDCut){
	      hdy[ehp][nstb][i]->Fill(dy);  
	      hxdy[ehp][i]->Fill(fcsx,dy); 
	      hydy[ehp][i]->Fill(fcsy,dy);
	      hpdy[ehp][i]->Fill(fcsrp,dy); 
	      hrdy[ehp][i]->Fill(fcsr,dy);
	    }
	    if(fcsx>0) {nstb=0;} else {nstb=1;}
	    if(fabs(dr)<mDCut) {
	      hdp[ehp][nstb][i]->Fill(rdp);  
	      hxdp[ehp][i]->Fill(fcsx,rdp); 
	      hydp[ehp][i]->Fill(fcsy,rdp);
	      hpdp[ehp][i]->Fill(fcsrp,rdp); 
	      hrdp[ehp][i]->Fill(fcsr,rdp);
	    }
	    if(fcsr<70) {nstb=0;} else {nstb=1;}
	    if(fabs(rdp)<mDCut) {
	      hdr[ehp][nstb][i]->Fill(dr);  
	      hxdr[ehp][i]->Fill(fcsx,dr); 
	      hydr[ehp][i]->Fill(fcsy,dr);
	      hpdr[ehp][i]->Fill(fcsrp,dr); 
	      hrdr[ehp][i]->Fill(fcsr,dr);
	    }
	  }//same side
	}//fwd track loop
      }//same event or mixed
    }//fcs cluster loop    
  } //event loop
  picoReader->Finish();
						
  for(int eh=0; eh<2; eh++){
    for(int nstb=0; nstb<2; nstb++){
      hdx[eh][nstb][2]= (TH1F*)hdx[eh][nstb][0]->Clone(); hdx[eh][nstb][2]->Add(hdx[eh][nstb][1],-1); hdx[eh][nstb][2]->SetLineColor(CUTCOL[2+nstb]);
      hdy[eh][nstb][2]= (TH1F*)hdy[eh][nstb][0]->Clone(); hdy[eh][nstb][2]->Add(hdy[eh][nstb][1],-1); hdy[eh][nstb][2]->SetLineColor(CUTCOL[2+nstb]);
      hdp[eh][nstb][2]= (TH1F*)hdp[eh][nstb][0]->Clone(); hdp[eh][nstb][2]->Add(hdp[eh][nstb][1],-1); hdp[eh][nstb][2]->SetLineColor(CUTCOL[2+nstb]);
      hdr[eh][nstb][2]= (TH1F*)hdr[eh][nstb][0]->Clone(); hdr[eh][nstb][2]->Add(hdr[eh][nstb][1],-1); hdr[eh][nstb][2]->SetLineColor(CUTCOL[2+nstb]);
    }
    hdx[eh][0][0]->Add(hdx[eh][1][0]);
    hdy[eh][0][0]->Add(hdy[eh][1][0]);
    hdp[eh][0][0]->Add(hdp[eh][1][0]);
    hdr[eh][0][0]->Add(hdr[eh][1][0]);
    hdx[eh][0][1]->Add(hdx[eh][1][1]);
    hdy[eh][0][1]->Add(hdy[eh][1][1]);
    hdp[eh][0][1]->Add(hdp[eh][1][1]);
    hdr[eh][0][1]->Add(hdr[eh][1][1]);
    hxdx[eh][2]= (TH2F*)hxdx[eh][0]->Clone(); hxdx[eh][2]->Add(hxdx[eh][1],-1);
    hydx[eh][2]= (TH2F*)hydx[eh][0]->Clone(); hydx[eh][2]->Add(hydx[eh][1],-1);
    hpdx[eh][2]= (TH2F*)hpdx[eh][0]->Clone(); hpdx[eh][2]->Add(hpdx[eh][1],-1);
    hrdx[eh][2]= (TH2F*)hrdx[eh][0]->Clone(); hrdx[eh][2]->Add(hrdx[eh][1],-1);
    hxdy[eh][2]= (TH2F*)hxdy[eh][0]->Clone(); hxdy[eh][2]->Add(hxdy[eh][1],-1);
    hydy[eh][2]= (TH2F*)hydy[eh][0]->Clone(); hydy[eh][2]->Add(hydy[eh][1],-1);
    hpdy[eh][2]= (TH2F*)hpdy[eh][0]->Clone(); hpdy[eh][2]->Add(hpdy[eh][1],-1);
    hrdy[eh][2]= (TH2F*)hrdy[eh][0]->Clone(); hrdy[eh][2]->Add(hrdy[eh][1],-1);
    hxdp[eh][2]= (TH2F*)hxdp[eh][0]->Clone(); hxdp[eh][2]->Add(hxdp[eh][1],-1);
    hydp[eh][2]= (TH2F*)hydp[eh][0]->Clone(); hydp[eh][2]->Add(hydp[eh][1],-1);
    hpdp[eh][2]= (TH2F*)hpdp[eh][0]->Clone(); hpdp[eh][2]->Add(hpdp[eh][1],-1);
    hrdp[eh][2]= (TH2F*)hrdp[eh][0]->Clone(); hrdp[eh][2]->Add(hrdp[eh][1],-1);
    hxdr[eh][2]= (TH2F*)hxdr[eh][0]->Clone(); hxdr[eh][2]->Add(hxdr[eh][1],-1);
    hydr[eh][2]= (TH2F*)hydr[eh][0]->Clone(); hydr[eh][2]->Add(hydr[eh][1],-1);
    hpdr[eh][2]= (TH2F*)hpdr[eh][0]->Clone(); hpdr[eh][2]->Add(hpdr[eh][1],-1);
    hrdr[eh][2]= (TH2F*)hrdr[eh][0]->Clone(); hrdr[eh][2]->Add(hrdr[eh][1],-1);
  }

  TCanvas *c1=new TCanvas("FCSTRK", "FCSTRK",50,0,1200,1200);
  TText* t;
  float lx=0.15,ly0=0.8,ldy=0.04, ly;

  c1->Divide(2,2);
  c1->cd(1); hdx[0][0][0]->Draw(); hdx[0][0][1]->Draw("same"); hdx[0][0][2]->Draw("same"); hdx[0][1][2]->Draw("same"); 
  c1->cd(2); hdy[0][0][0]->Draw(); hdy[0][0][1]->Draw("same"); hdy[0][0][2]->Draw("same"); hdy[0][1][2]->Draw("same");
  c1->cd(3); hdp[0][0][0]->Draw(); hdp[0][0][1]->Draw("same"); hdp[0][0][2]->Draw("same"); hdp[0][1][2]->Draw("same"); 
  c1->cd(4); hdr[0][0][0]->Draw(); hdr[0][0][1]->Draw("same"); hdr[0][0][2]->Draw("same"); hdr[0][1][2]->Draw("same");
  for(int j=0; j<4; j++){
    c1->cd(j+1);
    ly=ly0;
    for(int i=0; i<mNCut+1; i++){
      t=new TText(lx,ly,CUT[j][i]); 
      t->SetNDC(); t->SetTextColor(CUTCOL[i]); t->SetTextSize(0.03); 
      t->Draw();
      ly-=ldy;
    }
  }
  c1->SaveAs("fcsTrkMatchEcal.png");

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hdx[1][0][0]->Draw(); hdx[1][0][1]->Draw("same"); hdx[1][0][2]->Draw("same"); hdx[1][1][2]->Draw("same"); 
  c1->cd(2); hdy[1][0][0]->Draw(); hdy[1][0][1]->Draw("same"); hdy[1][0][2]->Draw("same"); hdy[1][1][2]->Draw("same");
  c1->cd(3); hdp[1][0][0]->Draw(); hdp[1][0][1]->Draw("same"); hdp[1][0][2]->Draw("same"); hdp[1][1][2]->Draw("same"); 
  c1->cd(4); hdr[1][0][0]->Draw(); hdr[1][0][1]->Draw("same"); hdr[1][0][2]->Draw("same"); hdr[1][1][2]->Draw("same");
  ly=ly0;
  for(int j=0; j<4; j++){
    c1->cd(j+1);
    ly=ly0;
    for(int i=0; i<mNCut+1; i++){      
      t=new TText(lx,ly,CUT[j][i]); 
      t->SetNDC(); t->SetTextColor(CUTCOL[i]); t->SetTextSize(0.03); 
      t->Draw();
      ly-=ldy;
    }
  }
  c1->SaveAs("fcsTrkMatchHcal.png");

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hXY[0]->Draw("colz");
  c1->cd(2); hXY[1]->Draw("colz");
  c1->cd(3); hXY[2]->Draw("colz");
  c1->cd(4); hXY[3]->Draw("colz");
  c1->SaveAs("fcsTrkxy.png");

  c1->Clear();
  c1->Divide(3,2);
  c1->cd(1); hET[0]->Draw();
  c1->cd(2); hET[1]->Draw();
  c1->cd(3); hET[2]->Draw();
  c1->cd(4); hE[0]->Draw();
  c1->cd(5); hE[1]->Draw();
  c1->cd(6); hE[2]->Draw();
  c1->SaveAs("fcsTrkEt.png");

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hxdx[0][2]->Draw("colz");
  c1->cd(2); hydx[0][2]->Draw("colz");
  c1->cd(3); hpdx[0][2]->Draw("colz");
  c1->cd(4); hrdx[0][2]->Draw("colz");
  c1->SaveAs("fcsTrkEcaldx.png");
  c1->Clear();

  c1->Divide(2,2);
  c1->cd(1); hxdy[0][2]->Draw("colz");
  c1->cd(2); hydy[0][2]->Draw("colz");
  c1->cd(3); hpdy[0][2]->Draw("colz");
  c1->cd(4); hrdy[0][2]->Draw("colz");
  c1->SaveAs("fcsTrkEcaldy.png");

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hxdp[0][2]->Draw("colz");
  c1->cd(2); hydp[0][2]->Draw("colz");
  c1->cd(3); hpdp[0][2]->Draw("colz");
  c1->cd(4); hrdp[0][2]->Draw("colz");
  c1->SaveAs("fcsTrkEcaldp.png");

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hxdr[0][2]->Draw("colz");
  c1->cd(2); hydr[0][2]->Draw("colz");
  c1->cd(3); hpdr[0][2]->Draw("colz");
  c1->cd(4); hrdr[0][2]->Draw("colz");
  c1->SaveAs("fcsTrkEcaldr.png");

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hxdx[1][2]->Draw("colz");
  c1->cd(2); hydx[1][2]->Draw("colz");
  c1->cd(3); hpdx[1][2]->Draw("colz");
  c1->cd(4); hrdx[1][2]->Draw("colz");
  c1->SaveAs("fcsTrkHcaldx.png");
  c1->Clear();

  c1->Divide(2,2);
  c1->cd(1); hxdy[1][2]->Draw("colz");
  c1->cd(2); hydy[1][2]->Draw("colz");
  c1->cd(3); hpdy[1][2]->Draw("colz");
  c1->cd(4); hrdy[1][2]->Draw("colz");
  c1->SaveAs("fcsTrkHcaldy.png");

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hxdp[1][2]->Draw("colz");
  c1->cd(2); hydp[1][2]->Draw("colz");
  c1->cd(3); hpdp[1][2]->Draw("colz");
  c1->cd(4); hrdp[1][2]->Draw("colz");
  c1->SaveAs("fcsTrkHcaldp.png");

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hxdr[1][2]->Draw("colz");
  c1->cd(2); hydr[1][2]->Draw("colz");
  c1->cd(3); hpdr[1][2]->Draw("colz");
  c1->cd(4); hrdr[1][2]->Draw("colz");
  c1->SaveAs("fcsTrkHcaldr.png");
}

