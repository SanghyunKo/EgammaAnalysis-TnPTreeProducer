// -*- C++ -*-
//
// Package:    EventCounter/SimpleEventCounter
// Class:      SimpleEventCounter
// 
/**\class SimpleEventCounter SimpleEventCounter.cc EventCounter/SimpleEventCounter/plugins/SimpleEventCounter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Francesco Micheli
//         Created:  Wed, 09 May 2018 13:01:25 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class SimpleEventCounter : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit SimpleEventCounter(const edm::ParameterSet&);
  ~SimpleEventCounter();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  edm::Service<TFileService> fs_;

  const bool isMC_;
  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle>> genptcToken_;
  // to keep track of the sum of weights
  TH1F* h_sumW;
  TH1F* h_sum4E;
  TH1F* h_sum4M;
  TH1F* h_sum2E2M;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SimpleEventCounter::SimpleEventCounter(const edm::ParameterSet& iConfig)
: isMC_(iConfig.getParameter<bool>("isMC")),
  generatorToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
  genptcToken_(consumes<edm::View<reco::GenParticle>>(edm::InputTag("prunedGenParticles")))
{
  //now do what ever initialization is needed
  usesResource("TFileService");
}


SimpleEventCounter::~SimpleEventCounter()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SimpleEventCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;



#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif

  double aWeight = 1.0;

  if (isMC_) {
    edm::Handle<GenEventInfoProduct> genInfo;
    iEvent.getByToken(generatorToken_, genInfo);
    double mcweight = genInfo->weight();

    aWeight = mcweight/std::abs(mcweight);

    std::vector<reco::GenParticleRef> promptMuons;
    std::vector<reco::GenParticleRef> promptElectrons;

    edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
    iEvent.getByToken(genptcToken_, genptcHandle);

    for (unsigned int idx = 0; idx < genptcHandle->size(); ++idx) {
      const auto& genptc = genptcHandle->refAt(idx);

      if ( ( std::abs(genptc->pdgId())==13 ) && genptc->fromHardProcessFinalState() )
        promptMuons.push_back(genptc.castTo<reco::GenParticleRef>());

      if ( ( std::abs(genptc->pdgId())==11 ) && genptc->fromHardProcessFinalState() )
        promptElectrons.push_back(genptc.castTo<reco::GenParticleRef>());
    }

    if (promptMuons.size()==4)
      h_sum4M->Fill(0.5,aWeight);

    if (promptElectrons.size()==4)
      h_sum4E->Fill(0.5,aWeight);

    if (promptElectrons.size()==2 && promptMuons.size()==2)
      h_sum2E2M->Fill(0.5,aWeight);
  }

  // To keep track of the sum of weights
  h_sumW->Fill(0.5,aWeight);
}


// ------------ method called once each job just before starting event loop  ------------
void 
SimpleEventCounter::beginJob()
{
  // to keep track of the sum of weights
  h_sumW = fs_->make<TH1F>("h_sumW", "h_sumW", 1,  0., 1.);
  h_sumW->Sumw2();

  h_sum4E = fs_->make<TH1F>("h_sum4E", "h_sum4E", 1,  0., 1.);
  h_sum4E->Sumw2();

  h_sum4M = fs_->make<TH1F>("h_sum4M", "h_sum4M", 1,  0., 1.);
  h_sum4M->Sumw2();

  h_sum2E2M = fs_->make<TH1F>("h_sum2E2M", "h_sum2E2M", 1,  0., 1.);
  h_sum2E2M->Sumw2();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SimpleEventCounter::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SimpleEventCounter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleEventCounter);

