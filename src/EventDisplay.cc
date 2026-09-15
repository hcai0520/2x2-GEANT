void EventDisplay() {
    // Open file and get trees
    TFile *f = new TFile("MCP_01.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Failed to open OutPut.root" << std::endl;
        return;
    }

    // Get data tree and config tree
    TTree *t_event = (TTree*)f->Get("PixelEDep");
    TTree *t_step = (TTree*)f->Get("ntuple");

    //TTree *config = (TTree*)f->Get("config");
    //if (!t || !config) {
    //    std::cerr << "Missing ntuple or config tree!" << std::endl;
    //    return;
    //}

    // Read detector configuration
    int sizeY= 70;
    int sizeZ= 60;
    double pixelSizeY=.4;
    double pixelSizeZ=.4;
    //config->SetBranchAddress("sizeY", &sizeY);
    //config->SetBranchAddress("sizeZ", &sizeZ);
    //config->SetBranchAddress("pixelSizeY", &pixelSizeY);
    //config->SetBranchAddress("pixelSizeZ", &pixelSizeZ);
    //config->GetEntry(0); // read the only entry

    // Compute histogram bounds and binning
    double yMin = -sizeY ;
    double yMax =  sizeY ;
    //double zMin = -sizeZ ;
    //double zMax =  sizeZ ;
    int yBins = sizeY / pixelSizeY;
    int zBins = sizeZ / pixelSizeZ;
    double zMin[4] = { -65.0,5.0, -65.0, 5.0 };
    double zMax[4] =  {-5.0, 65.0, -5.0, 65.0 };
    // Create canvas with 2x2 pads
    TCanvas *c = new TCanvas("c", "Event Display per Module", 1000, 1000);
    c->Divide(2, 2);

    
    for (int module = 0; module < 4; ++module) {
        c->cd(module + 1);
        
        // Create unique histogram name per module
        TString hname = Form("h%d", module);
        TString htitle = Form("Module %d;Z [cm];Y [cm];Edep [MeV]", module);
        TH2D *h = new TH2D(hname, htitle,
                           zBins, zMin[module], zMax[module],
                           yBins, yMin, yMax);

        TString drawCmd = Form("y:z>>%s", hname.Data());
        //TString cut = Form("VolName==%d", module);
        TString cut = Form("(VolName==%d) && (EnergyDep > 0.22)", module);
        t_event->Draw(drawCmd, "EnergyDep*(" + cut + ")", "COLZ");
        

        gPad->SetLogz();

        TString graphCut = Form("VolName==%d ", module);
        t_step->Draw("y:z", graphCut, "goff");  // silent draw

        int n = t_step->GetSelectedRows();
        if (n > 0) {
        TGraph *g = new TGraph(n, t_step->GetV2(), t_step->GetV1());  // x=z, y=y
        g->SetLineColor(kRed);
        g->SetLineWidth(1);
        g->SetMarkerStyle(20);
        g->SetMarkerColor(kBlack);
        g->SetMarkerSize(0.001);
        //g->Draw("PL SAME");
    }
    }

    c->Update();

}