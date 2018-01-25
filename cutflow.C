void cutflow.C(bool signal = false) {

  TChain chain_("Delphes");
  
  chain.Add("/Users/adrianthompson/physics/zprime/samples/ttbar/tt_b.root");
  chain.Add("/Users/adrianthompson/physics/zprime/samples/ttbar/tt_g.root");
  chain.Add("/Users/adrianthompson/physics/zprime/samples/ttbar/tt_j.root");
  chain.Add("/Users/adrianthompson/physics/zprime/samples/ttbar/tt_k.root");
  chain.Add("/Users/adrianthompson/physics/zprime/samples/ttbar/tt_l.root");
  chain.Add("/Users/adrianthompson/physics/zprime/samples/ttbar/tt_c.root");

  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain_);  
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMPT = treeReader->UseBranch("MissingET");
  
  TProfile *cutflow = new TProfile("cutflow","(2b & 2OS muons) & max(SBM)>170, MET<100, DHTLT<0;cut number;percent passing",4,-0.5,3.5);
  TEfficiency *efficiency1 = new TEfficiency("eff_maxSBMvsMmm","#epsilon(max(M(b#mu)) | preselection);M(#mu^{+}#mu^{-}) [GeV];efficiency",200,0.,1000.);
    

  // Begin event loop.
  for (Int_t entry = 0; entry < numberOfEntries; ++entry) {
    treeReader->ReadEntry(entry);
  
    
  
  
  } // End event loop.










}
