void CP_Plot_Class::Plot_Dimuon_Phi(){
				// ------- Create a canvas:
				TCanvas* p_c = new TCanvas("p_c","p_c",800,800);
				p_c->cd();

				// -------------- Top panel
				TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
				pad1->Draw(); pad1->cd();
				pad1->SetBottomMargin(0.00);
				gPad->RedrawAxis();
				// ------- Make atlerations and draw:
				TGaxis::SetMaxDigits(3);
				gPad->SetLogy();
				std::stringstream titlestring1;
				char title1[128];
				titlestring1 << "Dimuon #phi: " << this->histograms->mass_lower << "M" << this->histograms->mass_upper;
				titlestring1.getline(title1,128);
				this->histograms->p_mc_dimuon_phi_stack->SetNameTitle("dimuon_phi", title1);
				//this->histograms->p_mc_dimuon_phi_stack->SetMaximum(this->histograms->p_hist_dimuon_phi_Data->GetMaximum()*1.40);
				//this->histograms->p_mc_dimuon_phi_stack->SetMinimum(0);
				this->histograms->p_mc_dimuon_phi_stack->SetMaximum(5.0E8);
				this->histograms->p_mc_dimuon_phi_stack->SetMinimum(9E-2);
				this->histograms->p_hist_dimuon_phi_Data->SetLineColor(kBlack);
				this->histograms->p_hist_dimuon_phi_Data->SetMarkerStyle(20);
				this->histograms->p_hist_dimuon_phi_Data->SetMarkerSize(0.54);
				this->histograms->p_hist_dimuon_phi_Data->SetMarkerColor(kBlack);
				this->histograms->p_mc_dimuon_phi_stack->Draw("");
				this->histograms->p_hist_dimuon_phi_Data->Draw("same; pe");
				// ------- switch off x axis label for top plot
				this->histograms->p_mc_dimuon_phi_stack->GetYaxis()->SetTitle("No. Third Chain Dimuon events");
				this->histograms->p_mc_dimuon_phi_stack->GetYaxis()->SetTitleOffset(1.024);
				this->histograms->p_mc_dimuon_phi_stack->GetYaxis()->SetTitleSize(0.045);
				this->histograms->p_mc_dimuon_phi_stack->GetYaxis()->SetLabelSize(0.04);
				this->histograms->p_mc_dimuon_phi_stack->GetXaxis()->SetLabelSize(0);
				this->histograms->p_hist_dimuon_phi_Data->GetXaxis()->SetLabelSize(0);
				// ------- Set additional figure items:
				this->Set_Legend(0.49, 0.62);
				this->Set_Atlas_Internal_Text(0.16, 0.8, "Work in progress", kBlack);
				// ------- Return to parent canvas before making second pad:
				p_c->cd();
				// -------------- Bottom panel
				TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.3);
				pad2->SetTopMargin(0.00);
				pad2->Draw(); pad2->cd();
				pad2->SetBottomMargin(0.3);
				pad2->SetGrid(kTRUE);
				gStyle->SetOptStat(kFALSE);
				// ------- Make atlerations and draw:
				this->histograms->p_data_mc_dimuon_phi_ratio->SetMaximum(1.18);
				this->histograms->p_data_mc_dimuon_phi_ratio->SetMinimum(0.82);
				this->histograms->p_data_mc_dimuon_phi_ratio->Draw();
				this->histograms->p_data_mc_dimuon_phi_ratio->GetYaxis()->SetLabelSize(0.092);
				this->histograms->p_data_mc_dimuon_phi_ratio->GetXaxis()->SetLabelSize(0.10);
				this->histograms->p_data_mc_dimuon_phi_ratio->GetXaxis()->SetTitle("#phi_{#mu#mu}");
				this->histograms->p_data_mc_dimuon_phi_ratio->GetYaxis()->SetTitle("Data / MC");
				this->histograms->p_data_mc_dimuon_phi_ratio->GetXaxis()->SetTitleSize(0.097);
				this->histograms->p_data_mc_dimuon_phi_ratio->GetYaxis()->SetTitleSize(0.097);
				this->histograms->p_data_mc_dimuon_phi_ratio->GetXaxis()->SetTitleOffset(1.12);
				this->histograms->p_data_mc_dimuon_phi_ratio->GetYaxis()->SetTitleOffset(0.49);
				// ------- Fit the ratio:
				this->Fit_MC_Data_Ratio(0.18, 0.33, this->histograms->p_data_mc_dimuon_phi_ratio);
				// ------- Save the output:
				std::stringstream filestring;
				char out[128];
				filestring << "/data/armitage/ms3c_ntuples/" << this->extension << "analysis/thirdChain_" << this->histograms->file_i << "_" << this->histograms->mass_lower << "M" << this->histograms->mass_upper << "/" << this->subselection << "/dimuon_phi.eps";
				filestring.getline(out,128);
				p_c->Print(out);
				delete p_c;
				//delete pad1;
				//delete pad2;
} 
