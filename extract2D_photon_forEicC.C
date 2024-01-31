int extract2D_photon_forEicC(){

   TString dataname="test_pi0.root";

   TFile* simufile=new TFile(dataname.Data());

   TH2D* hist01_plot=(TH2D*)simufile->Get("pi0_pho1_rap_corr");
   TH2D* hist02_plot=(TH2D*)simufile->Get("pi0_rap_corr");

   TString dataname_forward="test_omega_pi0_v2.root";

   TFile* simufile_forward =new TFile(dataname_forward.Data());
   TH1D* hist03_plot=(TH1D*)simufile_forward->Get("Eff_HistFinal_p3");
   hist03_plot->SetName("hist03_plot");

   gStyle->SetOptStat(0);
   gStyle->SetOptDate(0);

   TF1 *ZDCLine = new TF1("ZDC","4.9",0,12);
   TF1 *FDTLine = new TF1("FDT",  "6",0,12);
   TF1 *EDTLine = new TF1("EDT","3.5",0,12);

   TBox *ZDCbox = new TBox(4.9,4.9,10,10);
   ZDCbox->SetFillColorAlpha(kRed, 0.35);
   ZDCbox->SetLineColor(kRed);

   TBox *FDTbox = new TBox(4.9,4.9,6,6);
   FDTbox->SetFillColorAlpha(kBlack, 0.35);
   FDTbox->SetLineColor(kBlack);

   TBox *EDTbox = new TBox(3,3,4.9,4.9);
   EDTbox->SetFillColorAlpha(kBlack, 0.35);
   EDTbox->SetLineColor(kBlack);

   TBox *CENbox = new TBox(0.0,0.0,3.,3.);
   CENbox->SetFillColorAlpha(kBlue, 0.35);
   CENbox->SetLineColor(kBlue);

/////////////////
   TBox *CENEDTbox = new TBox(0.0,3,3.,4.9);
   CENEDTbox->SetFillColorAlpha(kGreen+2, 0.35);
   CENEDTbox->SetLineColor(kGreen+2);

   TBox *EDTCENbox = new TBox(3,0.0,4.9,3.);
   EDTCENbox->SetFillColorAlpha(kGreen+2, 0.35);
   EDTCENbox->SetLineColor(kGreen+2);   


   TBox *CENZDCbox = new TBox(0.0,4.9,3.,10);
   CENZDCbox->SetFillColorAlpha(kYellow+2, 0.35);
   CENZDCbox->SetLineColor(kYellow+2);

   TBox *ZDCCENbox = new TBox(4.9,0.0,10,3.);
   ZDCCENbox->SetFillColorAlpha(kYellow+2, 0.35);
   ZDCCENbox->SetLineColor(kYellow+2);

   TCanvas* myC=new TCanvas("myC", "",10,10, 2700, 1000);

   TH2F* plothisto1 = new TH2F("ph1","ph1",1,0,2,1,0,5);
   TH2F* plothisto2 = new TH2F("ph2","ph2",1,0,2,1,0,5);
   TH2F* plothisto3 = new TH2F("ph3","ph3",1,0,2,1,0,5);

   TPad *pad1 = new TPad("pad1", "The pad 50% of the height",0.00,0.0,0.35,1.0,21);
   TPad *pad2 = new TPad("pad2", "The pad 50% of the height",0.35,0.0,0.70,1.0,21);
   TPad *pad3 = new TPad("pad3", "The pad 50% of the height",0.70,0.0,1.00,1.0,21);

   hist01_plot->GetXaxis()->SetRangeUser(0,9.99);
   hist01_plot->GetYaxis()->SetRangeUser(0.01,9.99);
   hist01_plot->GetZaxis()->SetRangeUser(1,150);
   hist02_plot->GetXaxis()->SetRangeUser(0.01,10);
   hist02_plot->GetYaxis()->SetRangeUser(0,9.99);
   hist02_plot->GetZaxis()->SetRangeUser(1,150);

   hist03_plot->GetXaxis()->SetRangeUser(0,6.5);

   pad1->SetFillColor(0);
   pad1->SetGrid(0,0);
   pad1->Draw();

   pad2->SetFillColor(0);
   pad2->SetGrid(0,0);
   pad2->Draw();

   pad3->SetFillColor(0);
   pad3->SetGrid(0,0);
   pad3->Draw();

   pad1->SetBottomMargin(0.2);
   pad1->SetTopMargin(0.1);
   pad1->SetLeftMargin(0.15);
   pad1->SetRightMargin(0.0);

   pad2->SetBottomMargin(0.2);
   pad2->SetTopMargin(0.1);
   pad2->SetLeftMargin(0.0);
   pad2->SetRightMargin(0.15);

   pad3->SetBottomMargin(0.2);
   pad3->SetTopMargin(0.1);
   pad3->SetLeftMargin(0.1);
   pad3->SetRightMargin(0.15);

   plothisto1->SetTitle("");
   hist01_plot->GetYaxis()->SetTitle("#eta (#gamma from #pi^{0} decay)");
   hist02_plot->GetYaxis()->SetTitle("#eta (#gamma from #pi^{0} decay)");

   hist02_plot->GetXaxis()->SetTitle("#eta (#gamma from #pi^{0} decay)");
   hist01_plot->GetXaxis()->SetTitle("#eta (#gamma from #omega decay)");

   hist01_plot->GetXaxis()->SetTitleOffset(1.5);
   hist02_plot->GetXaxis()->SetTitleOffset(1.5);

   hist03_plot->GetXaxis()->SetTitle("#eta");
   hist03_plot->GetYaxis()->SetTitle("Ratio [%]");
   hist03_plot->GetXaxis()->SetTitleOffset(1.5);

   pad1->cd();
   pad1->SetTickx(0);
   pad1->SetLogz();

   hist01_plot->Draw("COLZ");

   ZDCbox->Draw("Same");
   EDTbox->Draw("Same");
   CENbox->Draw("Same");

   CENEDTbox->Draw("Same");
   EDTCENbox->Draw("Same");

   CENZDCbox->Draw("Same");
   ZDCCENbox->Draw("Same");

   TLegend *legend=new TLegend(0.77,0.68,0.97,0.88);
   legend->AddEntry(ZDCbox,"ZDC","f");
   legend->AddEntry(EDTbox,"EDT ECal","f");
   legend->AddEntry(CENbox,"Central Detector","f");
   legend->AddEntry(CENEDTbox,"EDT ECal + Central","f");
   legend->AddEntry(CENZDCbox,"ZDC + Central","f");
   legend->SetFillColor(0);
   legend->SetBorderSize(1);
   legend->Draw();

   TLatex T1;
   T1.SetTextAlign(11);
   T1.SetTextSize(0.047);
   T1.SetTextColor(kBlack);
   T1.DrawLatex(0.75,9,"3.5#times20 GeV");
   T1.DrawLatex(0.75,8.1,"0 < Q^{2} < 1 GeV^{2}");

   pad2->cd();
   pad2->SetTickx(0);
   pad2->SetLogz();

   hist02_plot->Draw("COLZ");

   ZDCbox->Draw("Same");
   EDTbox->Draw("Same");
   CENbox->Draw("Same");

   CENEDTbox->Draw("Same");
   EDTCENbox->Draw("Same");

   CENZDCbox->Draw("Same");
   ZDCCENbox->Draw("Same");

   pad3->cd();
   pad3->SetTickx(0);

   hist03_plot->Draw("Hist");

   TString outputfile;//=filename;
   myC->Update();
   myC->Print("test.pdf");
   myC->Close();

   return 1;
}
