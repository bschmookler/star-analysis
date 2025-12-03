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
#include "StPicoEvent/StPicoMcTrack.h"
#include "StPicoEvent/StPicoMcVertex.h"
#include "StPicoEvent/StPicoFwdTrack.h"
#include "StPicoEvent/StPicoFwdVertex.h"
#include "StPicoEvent/StPicoFcsHit.h"
#include "StPicoEvent/StPicoFcsCluster.h"

int testTrackType = 0; // 0=Global, 1=BLC, 2=Primary, 3=FwdVertex
int mMaxEvent=400000;

//-------------------
bool inAcc(StPicoMcTrack *track) {
    // Example acceptance criteria
    return (track->nHitsFts() > 2 || track->nHitsStg() > 4);
}

//-------------------
void format_eff_plot(TH1* hNum, TH1 * hDen, TString name, TString title, int color = kBlack) {
    TH1 * hEff = (TH1*)hNum->Clone(name);
    hEff->Divide(hNum,hDen, 1, 1, "B");
    hEff->SetTitle(title);
    hEff->SetLineColor(color);
    hEff->SetLineWidth(4);
}

//-------------------
void track_qa(const Char_t *inFile = "infiles.lis", 
              int tt = -1, // track type
              TString label = ""
	     ){
   
    StPicoDstReader* picoReader = new StPicoDstReader(inFile);
    picoReader->Init();

    // This is a way if you want to spead up IO
    std::cout << "Explicit read status for some branches" << std::endl;
    picoReader->SetStatus("*",0);
    picoReader->SetStatus("Event",1);  
    picoReader->SetStatus("McTrack",1);
    picoReader->SetStatus("McVertex",1);
    picoReader->SetStatus("FwdTracks",1);
    picoReader->SetStatus("FwdVertices",1);
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

    if (tt >= 0) {
        testTrackType = tt; // Set the track type if specified
    }

    // Create output ROOT file
    TFile *fOutput = new TFile( TString::Format( "efficiency_trackType%d_%s.root", testTrackType, label.Data() ) , "RECREATE");
    fOutput->cd();

    // Define histograms
    int nEtaBins = 200;
    float etaMin = -5.0;
    float etaMax = 5.0;

    int nPtBins = 100;
    float ptMin = 0.0;
    float ptMax = 5.0;

    int nPhiBins = 100;
    float phiMin = -3.24;
    float phiMax = 3.24;

    int nMultBins = 100;
    float multMin = 0.0;
    float multMax = 100.0;

    TH1 * hMcPrimaryEta = new TH1F("hMcPrimaryEta", "MC Primary Eta", nEtaBins, etaMin, etaMax);
    TH1 * hMcPrimaryEtaInAcc = new TH1F("hMcPrimaryEtaInAcc", "MC Primary (in acc) Eta", nEtaBins, etaMin, etaMax);

    TH1 * hMcPrimaryPt = new TH1F("hMcPrimaryPt", "MC Primary Pt", nPtBins, ptMin, ptMax);
    TH1 * hMcPrimaryPtInAcc = new TH1F("hMcPrimaryPtInAcc", "MC Primary (in acc) Pt", nPtBins, ptMin, ptMax);

    TH1 * hMcPrimaryPhi = new TH1F("hMcPrimaryPhi", "MC Primary Phi", nPhiBins, phiMin, phiMax);
    TH1 * hMcPrimaryPhiInAcc = new TH1F("hMcPrimaryPhiInAcc", "MC Primary (in acc) Phi", nPhiBins, phiMin, phiMax);
    
    TH1 * hMcPrimaryMult = new TH1F("hMcPrimaryMult", "MC Primary Mult", nMultBins, multMin, multMax);
    TH1 * hMcPrimaryMultInAcc = new TH1F("hMcPrimaryMultInAcc", "MC Primary (in acc) Mult", nMultBins, multMin, multMax);

    TH1 * hMatchedMcEta = new TH1F("hMatchedMcEta", "Matched MC Eta", nEtaBins, etaMin, etaMax);
    TH1 * hMatchedMcPt = new TH1F("hMatchedMcPt", "Matched MC Pt", nPtBins, ptMin, ptMax);
    TH1 * hMatchedMcPhi = new TH1F("hMatchedMcPhi", "Matched MC Phi", nPhiBins, phiMin, phiMax);

    TH1 * hMatchedMcPrimaryEta = new TH1F("hMatchedMcPrimaryEta", "Matched MC Primary Eta", nEtaBins, etaMin, etaMax);
    TH1 * hMatchedMcPrimaryPt = new TH1F("hMatchedMcPrimaryPt", "Matched MC Primary Pt", nPtBins, ptMin, ptMax);
    TH1 * hMatchedMcPrimaryPhi = new TH1F("hMatchedMcPrimaryPhi", "Matched MC Primary Phi", nPhiBins, phiMin, phiMax);

    TH1 * hMatchedMcPrimaryMult = new TH1F("hMatchedMcPrimaryMult", "Matched MC Primary Mult", nMultBins, multMin, multMax);

    TH2 * hMatchedMcPrimaryEtaCharge = new TH2F("hMatchedMcPrimaryEtaCharge", "Matched MC Primary Eta vs Charge", nEtaBins, etaMin, etaMax, 3, -0.5, 2.5);
    TH2 * hMatchedMcPrimaryPtCharge = new TH2F("hMatchedMcPrimaryPtCharge", "Matched MC Primary Pt vs Charge", nPtBins, ptMin, ptMax, 3, -0.5, 2.5);
    TH2 * hMatchedMcPrimaryPhiCharge = new TH2F("hMatchedMcPrimaryPhiCharge", "Matched MC Primary Phi vs Charge", nPhiBins, phiMin, phiMax, 3, -0.5, 2.5);

    TH1 * hCurveResolution = new TH1F("hCurveResolution", "Curve Resolution", 100, -5.0, 5.0);
    TH2 * hCurveResolutionEta = new TH2F("hCurveResolutionEta", "Curve Resolution vs Eta", nEtaBins, etaMin, etaMax, 100, -5.0, 5.0);
    TH2 * hCurveResolutionPt = new TH2F("hCurveResolutionPt", "Curve Resolution vs Pt", nPtBins, ptMin, ptMax, 100, -5.0, 5.0);
    TH2 * hCurveResolutionPhi = new TH2F("hCurveResolutionPhi", "Curve Resolution vs Phi", nPhiBins, phiMin, phiMax, 100, -5.0, 5.0);

    TH1 * hQATruth = new TH1F("hQATruth", "QA Truth", 102, -1.0, 101.0);
    TH1 * hQATruthGoodChi2 = new TH1F("hQATruthGoodChi2", "QA Truth", 102, -1.0, 101.0);

    TH1 * hFwdTrackType = new TH1F("hFwdTrackType", "Forward Track Type", 4, -0.5, 3.5);

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

	// Reset counters
        size_t fwdMult = 0; // Forward multiplicity counter
        size_t mcMult = 0; // MC multiplicity counter
        size_t fwdMultReco = 0; // Reconstructed forward multiplicity counter

        // Loop over the McTracks
        int nMCPart = dst->numberOfMcTracks();
        for(int iMCPart = 0; iMCPart < nMCPart; iMCPart++){

            StPicoMcTrack* t1 = dst->mcTrack(iMCPart);
            if (!t1) continue;
	    
            TLorentzVector lv1 = t1->fourMomentum();
            
            StPicoMcVertex *vtx1 = dst->mcVertex(t1->idVtxStart() - 1);
            if (!vtx1) continue;

            if ( t1->idVtxStart() != 1 ){ // Only consider primary particles
                continue;
            }

            mcMult++; // Increment MC multiplicity counter
            hMcPrimaryEta->Fill(lv1.PseudoRapidity());
            hMcPrimaryPt->Fill(lv1.Pt());
            hMcPrimaryPhi->Fill(lv1.Phi());

            if (inAcc(t1)) {
                hMcPrimaryEtaInAcc->Fill(lv1.PseudoRapidity());
                hMcPrimaryPtInAcc->Fill(lv1.Pt());
                hMcPrimaryPhiInAcc->Fill(lv1.Phi());
                fwdMult++; // Increment Forward multiplicity counter
            }
        } //MCTracks Loop
 
        // Loop over the FwdTracks
        int nTrks = dst->numberOfFwdTracks();
        for(int iTrk = 0; iTrk < nTrks; iTrk++){

	    StPicoFwdTrack* fwdTrack = dst->fwdTrack(iTrk);
            if (!fwdTrack) continue;

            hFwdTrackType->Fill(fwdTrack->trackType());
            if ( fwdTrack->trackType() != testTrackType ) continue; // Only use tracks of chosen type

            // Check matched McTrack exists
            if ( fwdTrack->idTruth() < 1  || fwdTrack->idTruth() > nMCPart ) {
                // printf("FwdTrack with idTruth < 1 || idTruth > len(mcTracks): 
                // %d -> status = %d, chi2 = %f \n", fwdTrack->idTruth(), (int)fwdTrack->status(), fwdTrack->chi2());
                continue;
            }

            hQATruth->Fill(fwdTrack->qaTruth());

            // if ( fwdTrack->chi2() < 0.01 ) continue; // Skip tracks with chi2 < 0.01
            //FIXME -- Does not include chi2 cut
            hQATruthGoodChi2->Fill(fwdTrack->qaTruth());

            // Get matched McTrack
            StPicoMcTrack* mcTrack = dst->mcTrack(fwdTrack->idTruth() - 1);
            if (!mcTrack) { 
                printf("No matched MC track for FwdTrack with idTruth: %d\n", fwdTrack->idTruth());
                continue;
            }

            TLorentzVector lvMc = mcTrack->fourMomentum();
            hMatchedMcEta->Fill(lvMc.PseudoRapidity());
            hMatchedMcPt->Fill(lvMc.Pt());
            hMatchedMcPhi->Fill(lvMc.Phi());

            if ( mcTrack->idVtxStart() == 1 ){ // Require matched MCTrack to be primary particle
                hMatchedMcPrimaryEta->Fill(lvMc.PseudoRapidity());
                hMatchedMcPrimaryPt->Fill(lvMc.Pt());
                hMatchedMcPrimaryPhi->Fill(lvMc.Phi());

                int weightCharge = 1.0 - (int)(fwdTrack->charge() == mcTrack->charge());
                
                hMatchedMcPrimaryEtaCharge->Fill(lvMc.PseudoRapidity(), weightCharge);
                hMatchedMcPrimaryPtCharge->Fill(lvMc.Pt(), weightCharge);
                hMatchedMcPrimaryPhiCharge->Fill(lvMc.Phi(), weightCharge);

                float fwdCurve = 1.0 / (fwdTrack->momentum().Pt());
                float mcCurve = 1.0 / (lvMc.Pt());
                float curveRes = (fwdCurve - mcCurve) / mcCurve;
                hCurveResolution->Fill(curveRes);
                hCurveResolutionEta->Fill(lvMc.PseudoRapidity(), curveRes);
                hCurveResolutionPt->Fill(lvMc.Pt(), curveRes);
                hCurveResolutionPhi->Fill(lvMc.Phi(), curveRes);

                fwdMultReco++;
            }
        } // FwdTracks loop

        if ( fwdMult > 0){
            // compute the average efficiency vs multiplicity
            hMcPrimaryMult->Fill(mcMult);
            hMcPrimaryMultInAcc->Fill(fwdMult);
            hMatchedMcPrimaryMult->Fill(fwdMult, (float)(fwdMultReco/(float)fwdMult));
        }
    } // end of event loop

    // Matched over McPrimary
    format_eff_plot(hMatchedMcEta, hMcPrimaryEta, "Eta_Matched_McPrimary", "Efficiency vs #eta; #eta; Efficiency", kRed);
    format_eff_plot(hMatchedMcPt, hMcPrimaryPt, "Pt_Matched_McPrimary", "Efficiency vs p_{T}; p_{T} (GeV/c); Efficiency", kRed);
    format_eff_plot(hMatchedMcPhi, hMcPrimaryPhi, "Phi_Matched_McPrimary", "Efficiency vs #phi; #phi (rad); Efficiency", kRed);

    // Matched McPrimary over McPrimary 
    format_eff_plot(hMatchedMcPrimaryEta, hMcPrimaryEtaInAcc, "Eta_MatchedPrimary_McPrimary", "Efficiency vs #eta; #eta; Efficiency", kRed);
    format_eff_plot(hMatchedMcPrimaryPt, hMcPrimaryPtInAcc, "Pt_MatchedPrimary_McPrimary", "Efficiency vs p_{T}; p_{T} (GeV/c); Efficiency", kRed);
    format_eff_plot(hMatchedMcPrimaryPhi, hMcPrimaryPhiInAcc, "Phi_MatchedPrimary_McPrimary", "Efficiency vs #phi; #phi (rad); Efficiency", kRed);
    format_eff_plot(hMatchedMcPrimaryMult, hMcPrimaryMultInAcc, "Mult_MatchedPrimary_McPrimary", "Efficiency vs Multiplicity; Multiplicity; Efficiency", kRed);

    // Make the QATruth histograms
    hQATruth->SetTitle("QA Truth; QA Truth; Count");
    hQATruth->SetLineColor(kRed);
    hQATruth->SetFillColorAlpha(kRed, 0.3);
    hQATruth->SetLineWidth(2);
    hQATruthGoodChi2->SetTitle("QA Truth Good Chi2; QA Truth; Count");
    hQATruthGoodChi2->SetLineColor(kRed);
    hQATruthGoodChi2->SetFillColorAlpha(kRed, 0.3);
    hQATruthGoodChi2->SetLineWidth(2);

    // Create cumulative histograms
    hQATruth->Scale(1.0 / hQATruth->Integral()); // Normalize to 1
    hQATruthGoodChi2->Scale(1.0 / hQATruthGoodChi2->Integral()); 
    TH1 * hQaTruthC = hQATruth->GetCumulative(false);
    TH1 * hQaTruthGoodChi2C = hQATruthGoodChi2->GetCumulative(false);
    
    // Write histograms to the output file
    fOutput->Write();
    fOutput->Close();
}
