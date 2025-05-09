
void fwd_sim_ana1(){

	// Open the ROOT file
    	TFile *file = TFile::Open("SimpleTree_mudst.root", "READ");
    	if (!file || file->IsZombie()) {
        	std::cerr << "Error opening file!" << std::endl;
        	return 1;
    	}

    	// Get the TTree
    	TTree *tree = (TTree *)file->Get("data");
    	if (!tree) {
        	std::cerr << "Error: TTree 'tree' not found in file!" << std::endl;
       		return 1;
    	}

    	// Variables to hold the data
    	int mcpart_num;	
	float mcpart_px[1000];
	float mcpart_py[1000];
	float mcpart_pz[1000];
	int mcpart_charge[1000];
	float mcpart_Vtx_x[1000];
	float mcpart_Vtx_y[1000];
	float mcpart_Vtx_z[1000];
	float mcpart_VtxEnd_x[1000];
	float mcpart_VtxEnd_y[1000];
	float mcpart_VtxEnd_z[1000];

	int Trk_ntrks;
	int Trk_vtxindex[500];
	int Trk_nfitpoints[500];
	float Trk_px[500];
	float Trk_py[500];
	float Trk_pz[500];
	
	// Set branch addresses
	tree->SetBranchAddress("mcpart_num", &mcpart_num);
	tree->SetBranchAddress("mcpart_px", mcpart_px);
	tree->SetBranchAddress("mcpart_py", mcpart_py);
	tree->SetBranchAddress("mcpart_pz", mcpart_pz);
	tree->SetBranchAddress("mcpart_charge", mcpart_charge);
	tree->SetBranchAddress("mcpart_Vtx_x", mcpart_Vtx_x);
	tree->SetBranchAddress("mcpart_Vtx_y", mcpart_Vtx_y);
	tree->SetBranchAddress("mcpart_Vtx_z", mcpart_Vtx_z);
	tree->SetBranchAddress("mcpart_VtxEnd_x", mcpart_VtxEnd_x);
        tree->SetBranchAddress("mcpart_VtxEnd_y", mcpart_VtxEnd_y);
        tree->SetBranchAddress("mcpart_VtxEnd_z", mcpart_VtxEnd_z);

	tree->SetBranchAddress("Trk_ntrks",&Trk_ntrks);
	tree->SetBranchAddress("Trk_vtxindex",Trk_vtxindex);
	tree->SetBranchAddress("Trk_nfitpoints",Trk_nfitpoints);
	tree->SetBranchAddress("Trk_px",Trk_px);
	tree->SetBranchAddress("Trk_py",Trk_py);
	tree->SetBranchAddress("Trk_pz",Trk_pz);

	// Define Histograms
	TH1 *h1 = new TH1D("h1","Generated charged particles w/ beamline cut",100,-10,10);
	h1->GetXaxis()->SetTitle("Particle #eta");h1->GetXaxis()->CenterTitle();
	h1->SetLineColor(kBlue);h1->SetLineWidth(2);

	TH1 *h2 = new TH1D("h2","Generated charged particles w/ beamline cut and 2.5 < #eta < 4.0",100,0,4);
        h2->GetXaxis()->SetTitle("Particle P_{T} [GeV/c]");h2->GetXaxis()->CenterTitle();
        h2->SetLineColor(kBlue);h2->SetLineWidth(2);

	TH2 *h3 = new TH2D("h3","Generated charged particles w/ beamline cut and 2.5 < #eta < 4.0",200,-0.2,0.2,200,-0.2,0.2);
	h3->GetXaxis()->SetTitle("Particle P_{x} / P_{z}");h3->GetXaxis()->CenterTitle();
	h3->GetYaxis()->SetTitle("Particle P_{y} / P_{z}");h3->GetYaxis()->CenterTitle();
	h3->GetYaxis()->SetTitleOffset(0.95);

	// Generated charged particles w/ beamline cut in good acceptance region
	TH1 *h4 = new TH1D("h4","",100,0,4);
        h4->GetXaxis()->SetTitle("P_{T} [GeV/c]");h4->GetXaxis()->CenterTitle();
        h4->SetLineColor(kBlue);h4->SetLineWidth(2);

	// Transverse momentum of reconstructed primary forward tracks
	TH1 *g1 = new TH1D("g1","Transverse momentum of reconstructed primary forward tracks",100,0,4);
        g1->GetXaxis()->SetTitle("Track P_{T} [GeV/c]");g1->GetXaxis()->CenterTitle();
        g1->SetLineColor(kGreen);g1->SetLineWidth(2);

	// Transverse momentum of reconstructed primary forward tracks w/ nfitpoint > 5 cut
	TH1 *g2 = new TH1D("g2","",100,0,4);
        g2->GetXaxis()->SetTitle("Track P_{T} [GeV/c]");g2->GetXaxis()->CenterTitle();
        g2->SetLineColor(kRed);g2->SetLineWidth(2);

	TH2 *g3 = new TH2D("g3","Reconstructed primary forward tracks w/ >=6 fit points",200,-0.2,0.2,200,-0.2,0.2);
        g3->GetXaxis()->SetTitle("Track P_{x} / P_{z}");g3->GetXaxis()->CenterTitle();
        g3->GetYaxis()->SetTitle("Track P_{y} / P_{z}");g3->GetYaxis()->CenterTitle();
	g3->GetYaxis()->SetTitleOffset(0.95);

	//Reconstructed primary forward tracks w/ nfitpoint > 5 cut in good acceptance region
	TH1 *g4 = new TH1D("g4","",100,0,4);
        g4->GetXaxis()->SetTitle("P_{T} [GeV/c]");g4->GetXaxis()->CenterTitle();
        g4->SetLineColor(kRed);g4->SetLineWidth(2);

	// Loop over all entries in the tree
    	Long64_t nEntries = tree->GetEntries();
    	for (Long64_t iEvent = 0; iEvent < nEntries; ++iEvent) {
        	
		if(iEvent%10000==0) cout<<"Analysed "<<iEvent<<" events!"<<endl;
		tree->GetEntry(iEvent);

		// Loop over generated particles
		for(int iParticle = 0; iParticle < mcpart_num; iParticle++){

			TVector3 mom_vec(mcpart_px[iParticle],mcpart_py[iParticle],mcpart_pz[iParticle]);
			TVector3 start_vtx(mcpart_Vtx_x[iParticle],mcpart_Vtx_y[iParticle],mcpart_Vtx_z[iParticle]);
			TVector3 stop_vtx(mcpart_VtxEnd_x[iParticle],mcpart_VtxEnd_y[iParticle],mcpart_VtxEnd_z[iParticle]);

			// Only use charged particles
			if( mcpart_charge[iParticle]==0 ) continue;

			// Only use particles created in vertex area
			if( !(fabs(start_vtx.Z()) < 100 && start_vtx.Perp() < 1) ) continue;
	
			// Only use particles that leave vertex area
			if( fabs(stop_vtx.Z()) < 100 && stop_vtx.Perp() < 1 ) continue;

			h1->Fill(mom_vec.Eta());
		
			// Only keep particles with eta values that go into forward upgrade
			if( !(mom_vec.Eta() > 2.5 && mom_vec.Eta() < 4.0) ) continue;
			
			h2->Fill(mom_vec.Pt());
			h3->Fill( (mom_vec.Px()/mom_vec.Pz()) , (mom_vec.Py()/mom_vec.Pz()) );

			// Select a small good acceptance region
			if( !( fabs(mom_vec.Px()/mom_vec.Pz()) < 0.05 && 
			       (mom_vec.Py()/mom_vec.Pz()) > 0.075 &&
			       (mom_vec.Py()/mom_vec.Pz()) < 0.125 
			     )) continue;

			h4->Fill(mom_vec.Pt());			

		} //End loop over generated particles

		// Loop over reconstructed tracks
		for(int iTrk = 0; iTrk < Trk_ntrks; iTrk++){
			
			TVector3 trk_mom(Trk_px[iTrk],Trk_py[iTrk],Trk_pz[iTrk]);
			
			// Only keep primary tracks
			if( Trk_vtxindex[iTrk] == 255 ) continue;

			g1->Fill(trk_mom.Pt());

			// Remove tracks with < 6 hit points
			if ( Trk_nfitpoints[iTrk] < 6 ) continue;

			g2->Fill(trk_mom.Pt());
			g3->Fill( (trk_mom.Px()/trk_mom.Pz()) , (trk_mom.Py()/trk_mom.Pz()) );

			// Select a small good acceptance region
			if( !( fabs(trk_mom.Px()/trk_mom.Pz()) < 0.05 &&
                               (trk_mom.Py()/trk_mom.Pz()) > 0.075 &&
                               (trk_mom.Py()/trk_mom.Pz()) < 0.125
                             )) continue;		
	
			g4->Fill(trk_mom.Pt());			
	
		} //End loop over reconstructed tracks


	} //End loop over tree entries

	// Set Style
	gStyle->SetOptStat(0);

	// Make plots
	TCanvas *c1 = new TCanvas("c1");
	h1->Draw();

	TCanvas *c2 = new TCanvas("c2");
	h2->Draw();

	TCanvas *c3 = new TCanvas("c3");

	auto h2a = (TH1D*) h2->Clone("h2a");
	h2a->SetTitle("");
	h2a->GetXaxis()->SetTitle("P_{T} [GeV/c]");
	h2a->Draw();	

	g1->Draw("same");
	g2->Draw("same");

	TLegend *leg1 = new TLegend(0.25,0.6,0.75,0.85);
	leg1->SetBorderSize(0);
	leg1->SetTextSize(0.025);
	leg1->AddEntry(h2a,"Generated charged particles w/ beamline cut and 2.5 < #eta < 4.0");
	leg1->AddEntry(g1,"Reconstructed primary forward tracks");
	leg1->AddEntry(g2,"Reconstructed primary forward tracks w/ >=6 fit points");
	leg1->Draw();

	// Print some histogram statistics
	cout<<endl;
	cout<<"Number of generated paricles w/ beamline cut and 2.5 < #eta < 4.0 = "<<h2a->Integral()<<endl;
	cout<<"Number of reconstructed primary forward tracks = "<<g1->Integral()<<endl;
	cout<<"Number of reconstructed primary forward tracks w/ >=6 fit points = "<<g2->Integral()<<endl;
	cout<<endl;

	TCanvas *c4 = new TCanvas("c4");
	h3->Draw("colz");

	TCanvas *c5 = new TCanvas("c5");
        g3->Draw("colz");

	TCanvas *c6 = new TCanvas("c6");
	h4->Draw();
	g4->Draw("same");
	
	TLegend *leg2 = new TLegend(0.2,0.7,0.75,0.85);
        leg2->SetBorderSize(0);
        leg2->SetTextSize(0.0225);
        leg2->AddEntry(h4,"Generated charged particles w/ beamline cut in good acceptance region");
        leg2->AddEntry(g4,"Reconstructed primary forward tracks w/ >=6 fit points in good acceptance region");
        leg2->Draw();

	// Print plots to file
	c1->Print("fwd_sim_ana1.pdf[");
	c1->Print("fwd_sim_ana1.pdf");
	c2->Print("fwd_sim_ana1.pdf");
	c3->Print("fwd_sim_ana1.pdf");
	c4->Print("fwd_sim_ana1.pdf");
	c5->Print("fwd_sim_ana1.pdf");
	c6->Print("fwd_sim_ana1.pdf");
	c6->Print("fwd_sim_ana1.pdf]");
}

