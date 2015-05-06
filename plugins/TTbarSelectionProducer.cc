

#include "RecoBTag/PerformanceMeasurements/interface/TTbarSelectionProducer.h"

using namespace std;
using namespace edm;



TTbarSelectionProducer::TTbarSelectionProducer(const edm::ParameterSet& iConfig) :
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerColl"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonColl"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronColl"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetColl"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metColl")))  
{
   //register your products
   isData_            = iConfig.getParameter<bool > ("isData");
   verbose_           = iConfig.getParameter<int > ("verbose");

   //Configuration for electrons
   electron_cut_pt_   = iConfig.getParameter<double>        ("electron_cut_pt");
   electron_cut_eta_  = iConfig.getParameter<double>        ("electron_cut_eta");
   electron_cut_iso_  = iConfig.getParameter<double>        ("electron_cut_iso");

   //Configuration for muons
   muon_cut_pt_   = iConfig.getParameter<double>        ("muon_cut_pt");
   muon_cut_eta_  = iConfig.getParameter<double>        ("muon_cut_eta");
   muon_cut_iso_  = iConfig.getParameter<double>        ("muon_cut_iso");

   //Configuration for jets
   jet_cut_pt_   = iConfig.getParameter<double>        ("jet_cut_pt");
   jet_cut_eta_  = iConfig.getParameter<double>        ("jet_cut_eta");

   //Configuration for met
   met_cut_   = iConfig.getParameter<double>        ("met_cut");

   produces<int>();
   produces<vector<double>>();

   // some histograms
   hcheck_cutflow        = fs->make<TH1F>("hcheck_cutflow","Selection level", 8, -0.5, 7.5);
   hcheck_m_ee           = fs->make<TH1F>("hcheck_m_ee","M_{e e}",200,0.,1000);
   hcheck_m_emu          = fs->make<TH1F>("hcheck_m_emu","M_{e #mu}",200,0.,1000);
   hcheck_m_mumu         = fs->make<TH1F>("hcheck_m_mumu","M_{#mu #mu}",200,0.,1000);
   hcheck_met_ee         = fs->make<TH1F>("hcheck_met_ee","MET (ee channel)", 100,0., 500);
   hcheck_met_emu        = fs->make<TH1F>("hcheck_met_emu","MET (e #mu channel)", 100,0., 500);
   hcheck_met_mumu       = fs->make<TH1F>("hcheck_met_mumu","MET (#mu #mu channel)", 100,0., 500);
}


TTbarSelectionProducer::~TTbarSelectionProducer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
TTbarSelectionProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   if(verbose_>5) std::cout << "in TTbarSelectionProducer::produce " << std::endl;

   //check trigger
   edm::Handle<edm::TriggerResults> triggerBits;
   iEvent.getByToken(triggerBits_, triggerBits);
   
   
   //------------------------------------------
   //get beam spot
   //------------------------------------------
   const reco::BeamSpot* bs = 0;
   edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
   iEvent.getByLabel("offlineBeamSpot", recoBeamSpotHandle);
   bs = recoBeamSpotHandle.product();
   
   //------------------------------------------
   //Selection of muons
   //------------------------------------------

   std::vector< TLorentzVector  > p4Muon;
   std::vector< int  > chargeMuon;
   edm::Handle<pat::MuonCollection> muHa;
   iEvent.getByToken(muonToken_, muHa);
   for (const pat::Muon &mu : *muHa) 
     { 
     
       bool passKin( mu.pt() > muon_cut_pt_  && fabs(mu.eta()) < muon_cut_eta_ );
       if(!passKin) continue;

       bool passID ( mu.isPFMuon() 
		     && mu.isGlobalMuon() 
		     && mu.isTrackerMuon()
		     && mu.normChi2() < 10 
		     && mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 
		     && mu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0 
		     && mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0 
		     && mu.numberOfMatchedStations() > 1 
		     && mu.innerTrack()->dxy(*bs) < 0.02
		     );
       if(!passID) continue;

       double nhIso   = mu.neutralHadronIso();
       double puchIso = mu.puChargedHadronIso();
       double chIso   = mu.chargedHadronIso() ;
       double gIso    = mu.photonIso() ;
       double relIso  = (TMath::Max(Float_t(nhIso+gIso-0.5*puchIso),Float_t(0.))+chIso)/mu.pt();
       bool passIso( relIso < muon_cut_iso_ );
       if(!passIso) continue;
	       
       TLorentzVector themuon;
       themuon.SetPtEtaPhiM(mu.pt(), mu.eta(),  mu.phi(), 0.);
       p4Muon.push_back(themuon);
       chargeMuon.push_back(mu.charge());
     }
   if(verbose_>5) std::cout << "\t Selected  " << chargeMuon.size() << " muons" << std::endl;


   //------------------------------------------
   //Selection of electrons
   //------------------------------------------


   std::vector< TLorentzVector  > p4Elec;
   std::vector< int  > chargeElec;
   edm::Handle<pat::ElectronCollection> elHa;
   iEvent.getByToken(electronToken_, elHa);
   for (const pat::Electron &el : *elHa) 
     {    
       double theta = 2*atan(exp(-1*el.superCluster()->eta()));
       double ET_SC = el.superCluster()->energy()*sin(theta);
       bool passKin( el.pt() > electron_cut_pt_  && fabs(el.eta()) < electron_cut_eta_ && ET_SC>15);
       if(!passKin) continue;

       bool passID( el.ecalDrivenSeed()
		    && el.gsfTrack().isNonnull()
		    && fabs(el.gsfTrack()->dxy(*bs)) < 0.04
		    && ET_SC > 15 
		    && el.passConversionVeto() 
		    && el.gsfTrack()-> hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <1);

       const std::vector< std::pair<std::string,float> > patids  = el.electronIDs();
       for (unsigned int i=0;i<patids.size();i++)
	 {
	   if(patids[i].first != "mvaTrigV0") continue;
	   passID &= (patids[i].second > 0.5 );
	 }
       if(!passID) continue;

       double nhIso   = el.neutralHadronIso();
       double puchIso = el.puChargedHadronIso();
       double chIso   = el.chargedHadronIso() ;
       double gIso    = el.photonIso() ;
       double relIso  = (TMath::Max(Float_t(nhIso+gIso-0.5*puchIso),Float_t(0.))+chIso)/el.pt();
       bool passIso( relIso < electron_cut_iso_ );
       if(!passIso) continue;
	     
       TLorentzVector theelectron;
       theelectron.SetPtEtaPhiM(el.pt(), el.eta(),  el.phi(), 0.);
       
       p4Elec.push_back(theelectron);
       chargeElec.push_back(el.charge());
     }
   if(verbose_>5) std::cout << "\t Selected " << chargeElec.size() << " electrons" << std::endl;
   
   
   //------------------------------------------
   //Selection of jet
   //------------------------------------------
   edm::Handle<pat::JetCollection> jetHa;
   iEvent.getByToken(jetToken_, jetHa);
   std::vector< TLorentzVector  > p4Jet;
   for (const pat::Jet &j : *jetHa) 
     {
       bool passKin( j.pt() > jet_cut_pt_  && fabs(j.eta()) < jet_cut_eta_ );
       if(!passKin) continue;
	   
       // check overlap with electron and muon
       TLorentzVector thejet;
       thejet.SetPtEtaPhiM(j.pt(), j.eta(),  j.phi(), 0.);
       double deltaRmu = 10000;
       double deltaRel = 10000;
       for(unsigned int imu=0; imu< p4Muon.size(); imu++)
	 {
	   double deltaR = thejet.DeltaR(p4Muon[imu]);
	   if(deltaR < deltaRmu) deltaRmu = deltaR;
	 }
       for(unsigned int iel=0; iel< p4Elec.size(); iel++)
	 {
	   double deltaR = thejet.DeltaR(p4Elec[iel]);
	   if(deltaR < deltaRel) deltaRel = deltaR;
	 }
       bool hasOverlap(deltaRmu<0.5 || deltaRel<0.5);
       if(!hasOverlap) continue;

       
       p4Jet.push_back(thejet);
     }
   if(verbose_>5) std::cout << "\t Selected "<< p4Jet.size() << " jets" << std::endl;
   

   edm::Handle<pat::METCollection> metHa;
   iEvent.getByToken(metToken_, metHa);

   //std::cout << "get met : done " << std::endl;
   int channel = -1;
   int ind_cutflow=0;
   std::vector< TLorentzVector  > TheLeptons;
   double Mll_2 = -1;
   double themet = -1.;

   bool passSel = false;
   int nLeps(p4Elec.size()+p4Muon.size());
   int nSelJets(p4Jet.size());
   if (nLeps >=1) ind_cutflow++;
   if( nLeps >=2){
     ind_cutflow++;
     int idxLept1 = -1;
     int idxLept2 = -1;
     //std::cout << " Lepton : " << p4Elec.size()+p4Muon.size() << std::endl;
     GetLeptonPair(p4Elec, p4Muon, chargeElec, chargeMuon, idxLept1, idxLept2, channel);

     if(channel >=0){

       ind_cutflow++;
       double Mll = -1;
       if(channel == 0) Mll = (p4Elec[idxLept1]+p4Elec[idxLept2]).M();
       if(channel == 1) Mll = (p4Muon[idxLept1]+p4Muon[idxLept2]).M();
       if(channel == 2) Mll = (p4Elec[idxLept1]+p4Muon[idxLept2]).M();

       if(channel == 0)  {
         TheLeptons.push_back(p4Elec[idxLept1]);
         TheLeptons.push_back(p4Elec[idxLept2]);
       }
       else if (channel == 1) {
         TheLeptons.push_back(p4Muon[idxLept1]);
         TheLeptons.push_back(p4Muon[idxLept2]);
       }
      else if (channel == 2) {
	TheLeptons.push_back(p4Elec[idxLept1]);
	TheLeptons.push_back(p4Muon[idxLept2]);
      }
       
       Mll_2 = (TheLeptons[0]+TheLeptons[1]).M();
       if (fabs(Mll_2-Mll)>0.01) std::cout << " Mll_2 " << Mll_2 << " Mll " << Mll << std::endl;
          
       const pat::MET    *met = 0;
       met = &(metHa->front());
       themet = met->pt();
       
       if( Mll > 20 && ( channel ==2 || ( channel <=1 && fabs(Mll-91)>15)) ) {
	   ind_cutflow++;
	   if (nSelJets >=2) {
	     ind_cutflow++;
	     if (themet >  met_cut_  ||  channel ==2) {
	       passSel = true;
	       ind_cutflow++;
	       if (channel==0) {
		 hcheck_m_ee->Fill(Mll)   ;
		 hcheck_met_ee->Fill(themet) ;
	       }
	       else if (channel==1) {
		 hcheck_m_mumu->Fill(Mll)   ;
		 hcheck_met_mumu->Fill(themet) ;
	       }
	       else if (channel==2) {
                 hcheck_m_emu->Fill(Mll)   ;
                 hcheck_met_emu->Fill(themet) ;
		 if(verbose_>5)
		   {
		     cout << "\t lepton 1 " << TheLeptons[0].Pt() << " " << TheLeptons[0].Eta() << endl
			  << "\t lepton 2 " << TheLeptons[1].Pt() << " " << TheLeptons[1].Eta() << endl;
		     for(unsigned int ij=0; ij< p4Jet.size(); ij++) {
		       cout << "\t jet " << ij << "   " << p4Jet[ij].Pt() << " " << p4Jet[ij].Eta() << endl;
		     }
		   }
	       }
	     }
	   }
	 }
     }
   }
   

   if(!passSel) channel = -1;
   hcheck_cutflow->Fill(0)        ;
   if (ind_cutflow>0) {
     for (int ii=1; ii<=ind_cutflow; ii++) {
       hcheck_cutflow->Fill(ii)        ;
     }
   }
   
   std::auto_ptr<int > pOut( new int (channel) );
   iEvent.put(pOut);
   
   vector<double> thelep_and_met;
   if (channel>=0) {
     thelep_and_met.push_back(TheLeptons[0].Pt());
     thelep_and_met.push_back(TheLeptons[0].Eta());
     thelep_and_met.push_back(TheLeptons[0].Phi());
     thelep_and_met.push_back(TheLeptons[1].Pt());
     thelep_and_met.push_back(TheLeptons[1].Eta());
     thelep_and_met.push_back(TheLeptons[1].Phi());
     thelep_and_met.push_back(themet);
     thelep_and_met.push_back(Mll_2);
   }
   std::auto_ptr<std::vector<double>> pOut2 (new std::vector<double> (thelep_and_met) );
   iEvent.put(pOut2);
   /*
     std::auto_ptr<std::vector<TLorentzVector>> pOut2 (new std::vector<TLorentzVector> (TheLeptons) );
     iEvent.put(pOut2);

     std::auto_ptr<double > pOut3( new double (themet) );
     iEvent.put(pOut3);
   */
   
   //std::cout << "print output : done " << std::endl;
}

// ------------ method called once each job just before starting event loop  ------------
void
TTbarSelectionProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
TTbarSelectionProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void
TTbarSelectionProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
TTbarSelectionProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
TTbarSelectionProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
TTbarSelectionProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TTbarSelectionProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



void
TTbarSelectionProducer::GetLeptonPair(
                          std::vector<TLorentzVector> elec_in, std::vector<TLorentzVector> muon_in,
                          std::vector<int> elec_charge, std::vector<int> muon_charge,
			  int &idxLept1, int &idxLept2, int &thechannel){


  float sum_pT_ee = 0.;
  bool pass_elec = false;
  int ie1 = -1;
  int ie2 = -1;
  if (elec_in.size () >= 2) {
    for (unsigned int i = 0; i < elec_in.size (); i++) {
      for (unsigned int j = i + 1; j < elec_in.size (); j++) {
	if (pass_elec)
	  continue;
	if ( elec_charge[i] != elec_charge[j] ){
	  pass_elec = true;
	  sum_pT_ee = elec_in[i].Pt () + elec_in[j].Pt ();
	  ie1 = i;
	  ie2 = j;
	}
      }
    }
  }

  float sum_pT_mumu = 0.;
  bool pass_muon = false;
  int imu1 = -1;
  int imu2 = -1;
  if (muon_in.size () >= 2) {
    for (unsigned int i = 0; i < muon_in.size (); i++) {
      for (unsigned int j = i + 1; j < muon_in.size (); j++) {
	if (pass_muon)
	  continue;
	if ( muon_charge[i] != muon_charge[j] ){
	  pass_muon = true;
	  sum_pT_mumu = muon_in[i].Pt () + muon_in[j].Pt ();
	  imu1 = i;
	  imu2 = j;
	}
      }
    }
  }


  float sum_pT_emu_start = 0.;
  float sum_pT_emu = 0.;
  int je1 = -1;
  int jmu2 = -1;
  if (muon_in.size () >= 1 && elec_in.size () >= 1) {
    for (unsigned int i = 0; i < muon_in.size (); i++) {
      for (unsigned int j = 0; j < elec_in.size (); j++) {
	if ( (muon_charge[i] != elec_charge[j]) ){
	  sum_pT_emu = muon_in[i].Pt () + elec_in[j].Pt ();
	  if (sum_pT_emu > sum_pT_emu_start) {
	    sum_pT_emu_start = sum_pT_emu;
	    je1 = j;
	    jmu2 = i;
	  }
	}
      }
    }
  }


  float sum[3] = { sum_pT_ee, sum_pT_mumu, sum_pT_emu };
  int sortedIndices[3];
  TMath::Sort (3, sum, sortedIndices);
  if (sortedIndices[0] == 0 && sum_pT_ee != 0.) {
    idxLept1 = ie1;
    idxLept2 = ie2;
    thechannel = 0;
  }
  else if (sortedIndices[0] == 1 && sum_pT_mumu != 0.) {
    idxLept1 = imu1;
    idxLept2 = imu2;
    thechannel = 1;
  }
  else if (sortedIndices[0] == 2 && sum_pT_emu != 0.) {
    idxLept1 = je1;
    idxLept2 = jmu2;
    thechannel = 2;
  }






}




//define this as a plug-in
DEFINE_FWK_MODULE(TTbarSelectionProducer);
