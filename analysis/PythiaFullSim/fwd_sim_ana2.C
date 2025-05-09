void fwd_sim_ana2(){

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
	int Trk_ntrks;
        int Trk_vtxindex[500];
        int Trk_nfitpoints[500];
        float Trk_chi2[500];
	float Trk_px[500];
        float Trk_py[500];
        float Trk_pz[500];
	float Trk_proj_ecal_x[500];
	float Trk_proj_ecal_y[500];
  	float Trk_proj_ecal_z[500];
  	float Trk_proj_hcal_x[500];
  	float Trk_proj_hcal_y[500];
  	float Trk_proj_hcal_z[500];

	int Cal_nhits;
  	int Cal_detid[5000];
  	int Cal_hitid[5000];
  	float Cal_hit_energy[5000];
  	float Cal_hit_posx[5000];
  	float Cal_hit_posy[5000];
  	float Cal_hit_posz[5000];

        // Set branch addresses	
	tree->SetBranchAddress("Trk_ntrks",&Trk_ntrks);
        tree->SetBranchAddress("Trk_vtxindex",Trk_vtxindex);
        tree->SetBranchAddress("Trk_nfitpoints",Trk_nfitpoints);
	tree->SetBranchAddress("Trk_chi2",Trk_chi2);
        tree->SetBranchAddress("Trk_px",Trk_px);
        tree->SetBranchAddress("Trk_py",Trk_py);
        tree->SetBranchAddress("Trk_pz",Trk_pz);
	tree->SetBranchAddress("Trk_proj_ecal_x",Trk_proj_ecal_x);
	tree->SetBranchAddress("Trk_proj_ecal_y",Trk_proj_ecal_y);
	tree->SetBranchAddress("Trk_proj_ecal_z",Trk_proj_ecal_z);
	tree->SetBranchAddress("Trk_proj_hcal_x",Trk_proj_hcal_x);
        tree->SetBranchAddress("Trk_proj_hcal_y",Trk_proj_hcal_y);
        tree->SetBranchAddress("Trk_proj_hcal_z",Trk_proj_hcal_z);

	tree->SetBranchAddress("Cal_nhits",&Cal_nhits);
	tree->SetBranchAddress("Cal_detid",Cal_detid);
	tree->SetBranchAddress("Cal_hitid",Cal_hitid);
	tree->SetBranchAddress("Cal_hit_energy",Cal_hit_energy);
	tree->SetBranchAddress("Cal_hit_posx",Cal_hit_posx);	
	tree->SetBranchAddress("Cal_hit_posy",Cal_hit_posy);
	tree->SetBranchAddress("Cal_hit_posz",Cal_hit_posz);	

        // Define Histograms

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

	// Loop over all entries in the tree
        Long64_t nEntries = tree->GetEntries();
	cout << "Analysing " << nEntries << " events!" << endl << endl;        

	for (Long64_t iEvent = 0; iEvent < nEntries; ++iEvent) {

                if(iEvent%10000==0) cout<<"Analysed "<<iEvent<<" events!"<<endl;
                tree->GetEntry(iEvent);

		// Vector to hold indices of good tracks
		std::vector<int> good_trk_index;

		// Loop over reconstructed tracks and find set of 'good' tracks
		for(int iTrk = 0; iTrk < Trk_ntrks; iTrk++){
			if( Trk_vtxindex[iTrk] != 255 && Trk_nfitpoints[iTrk] > 5 ){
				good_trk_index.push_back(iTrk);
			}
		}	

		// Loop over FCS hits
		for (int ihit = 0; ihit < Cal_nhits; ++ihit){
	
			// FCS ECal
			if( Cal_detid[ihit] == 0 || Cal_detid[ihit] == 1 ){

				// Get tower information
				auto tower_energy = Cal_hit_energy[ihit];
				auto tower_det = Cal_detid[ihit];
				auto tower_row = int(Cal_hitid[ihit]/22);
				auto tower_column = Cal_hitid[ihit]%22;
				auto tower_x = Cal_hit_posx[ihit];
				auto tower_y = Cal_hit_posy[ihit];

				h1a->Fill(tower_energy);
				h1b->Fill(tower_energy);

				// Require hit to be close to a track projection
				// Full width of ECal towers ~5.57cm
				bool matched_to_trk = false;
				for(int iTrk : good_trk_index){
					if( std::abs(Trk_proj_ecal_x[iTrk]-tower_x) < 3.0 &&
					    std::abs(Trk_proj_ecal_y[iTrk]-tower_y) < 3.0){
						matched_to_trk = true;
						break;
					}
				}

				if(matched_to_trk) h1c->Fill(tower_energy);

				// Require small relative energy deposit in adjacent towers
				bool keep_tower = true;
				for (int iAdj = 0; iAdj < Cal_nhits; ++iAdj){

					// Don't use the same tower
					if(iAdj == ihit) continue;					

					auto adj_energy = Cal_hit_energy[iAdj];
					auto adj_det = Cal_detid[iAdj];
					auto adj_row = int(Cal_hitid[iAdj]/22);
					auto adj_column = Cal_hitid[iAdj]%22;

					if( adj_det == tower_det &&
					    abs(adj_row - tower_row) <=1 &&
					    abs(adj_column - tower_column) <=1 &&
					    (adj_energy / tower_energy) > 0.1 ){
						keep_tower = false;
						break;
					}
				}

				if(matched_to_trk && keep_tower) h1d->Fill(tower_energy);
				if(keep_tower) h1e->Fill(tower_energy);
				
			} // End if FCS ECal
			
			// FCS HCal
			else if( Cal_detid[ihit] == 2 || Cal_detid[ihit] == 3 ){
		
				// Get tower information
                                auto tower_energy = Cal_hit_energy[ihit];
                                auto tower_det = Cal_detid[ihit];
                                auto tower_row = int(Cal_hitid[ihit]/13);
                                auto tower_column = Cal_hitid[ihit]%13;
                                auto tower_x = Cal_hit_posx[ihit];
                                auto tower_y = Cal_hit_posy[ihit];	

				h2a->Fill(tower_energy);
                                h2b->Fill(tower_energy);

				// Require hit to be close to a track projection
				// Full width of HCal towers ~9.99cm
				bool matched_to_trk = false;
                                for(int iTrk : good_trk_index){
                                        if( std::abs(Trk_proj_hcal_x[iTrk]-tower_x) < 5.0 &&
                                            std::abs(Trk_proj_hcal_y[iTrk]-tower_y) < 5.0){
                                                matched_to_trk = true;
                                                break;
                                        }
                                }

				if(matched_to_trk) h2c->Fill(tower_energy);
				
				// Require small relative energy deposit in adjacent towers
                                bool keep_tower = true;
                                for (int iAdj = 0; iAdj < Cal_nhits; ++iAdj){

                                        // Don't use the same tower
                                        if(iAdj == ihit) continue;

                                        auto adj_energy = Cal_hit_energy[iAdj];
                                        auto adj_det = Cal_detid[iAdj];
                                        auto adj_row = int(Cal_hitid[iAdj]/13);
                                        auto adj_column = Cal_hitid[iAdj]%13;

                                        if( adj_det == tower_det &&
                                            abs(adj_row - tower_row) <=1 &&
                                            abs(adj_column - tower_column) <=1 &&
                                            (adj_energy / tower_energy) > 0.1 ){
                                                keep_tower = false;
                                                break;
                                        }
                                }

                                if(matched_to_trk && keep_tower) h2d->Fill(tower_energy);
				if(keep_tower) h2e->Fill(tower_energy);

			} // End if FCS HCal
	
		} // End loop over FCS hits

	} // End loop over events

	// Set Style
        gStyle->SetOptStat(0);

        // Make plots
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

	// Print plots to file
        c1a->Print("fwd_sim_ana2.pdf[");
        c1a->Print("fwd_sim_ana2.pdf");
        c1b->Print("fwd_sim_ana2.pdf");
	c1c->Print("fwd_sim_ana2.pdf");
	c1d->Print("fwd_sim_ana2.pdf");
	c2a->Print("fwd_sim_ana2.pdf");
        c2b->Print("fwd_sim_ana2.pdf");
	c2c->Print("fwd_sim_ana2.pdf");
	c2d->Print("fwd_sim_ana2.pdf");
        c2d->Print("fwd_sim_ana2.pdf]");	

}

