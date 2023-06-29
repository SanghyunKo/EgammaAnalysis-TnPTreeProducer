#include "MiniAODTriggerCandProducer.h"
#include "EgammaAnalysis/TnPTreeProducer/plugins/WriteValueMap.h"

#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

template <>
bool MiniAODTriggerCandProducer<reco::GsfElectron, trigger::TriggerObject>::onlineOfflineMatchingRECO(reco::GsfElectronRef ref, 
												      const std::vector<trigger::TriggerObject>* triggerObjects,
												      const trigger::Keys* keys, float dRmin) {

  for (const auto & key : *keys) {
    float dR = deltaR2(ref->superCluster()->position().eta(), ref->superCluster()->position().phi(),
		       (*triggerObjects)[key].eta(), (*triggerObjects)[key].phi());
    //std::cout << "dr" << dR << std::endl;
    if (dR < dRmin*dRmin) {
      return true;
    }
  }

  return false;
}

template <>
MiniAODTriggerCandProducer<reco::GsfElectron, trigger::TriggerObject>::MiniAODTriggerCandProducer(const edm::ParameterSet& iConfig ):
  filterNames_(iConfig.getParameter<std::vector<std::string> >("filterNames")),
  inputs_(consumes<TRefVector>(iConfig.getParameter<edm::InputTag>("inputs"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  //triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  //triggerObjects_(consumes<U>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerEvent_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("objects"))),
  dRMatch_(iConfig.getParameter<double>("dR")),
  isAND_(iConfig.getParameter<bool>("isAND")) {
  
  produces<TRefVector>();
}

template <>
void MiniAODTriggerCandProducer<reco::GsfElectron, trigger::TriggerObject>::produce(edm::Event &iEvent, const edm::EventSetup &eventSetup) {
  
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<trigger::TriggerEvent> trEv;
  iEvent.getByToken(triggerEvent_, trEv);
  const trigger::TriggerObjectCollection& triggerObjects(trEv->getObjects());
  
  edm::Handle<TRefVector> inputs;

  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(inputs_, inputs);

  // Create the output collection
  std::unique_ptr<TRefVector> outColRef(new TRefVector);
  
  if (!triggerBits.isValid()) {
    LogDebug("") << "TriggerResults product not found - returning result=false!";
    return;
  }

  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerBits);
  if (triggerNamesID_ != triggerNames.parameterSetID()) {
    triggerNamesID_ = triggerNames.parameterSetID();
    init(*triggerBits, triggerNames);
  } 

  for (size_t i=0; i<inputs->size(); i++) {
    bool saveObj = true;
    TRef ref = (*inputs)[i];

    //    std::cout << "REF:" << ref->eta() << " " << ref->phi() << " " << ref->et() << std::endl;
    if (filterNames_.size() > 0) {
      unsigned int moduleFilterIndex = 9999;
      if (filterNames_[0] != "") {
	//std::cout << "DIVERSO DA 0" << std::endl;
	for (int i=0;i<trEv->sizeFilters(); i++) {
	  if (trEv->filterLabel(i) == filterNames_[0]) {
	    moduleFilterIndex = i;
	    //std::cout << "myIndex:" << i << std::endl;
	    break;
	  }
	}
	//unsigned int moduleFilterIndex = trEv->filterIndex(filterNames_[0]);
	//std::cout << "Filter: " << filterNames_[0] << std::endl;//" " << moduleFilterIndex << " " << trEv->sizeFilters() << std::endl;
	
	if (moduleFilterIndex+1 > trEv->sizeFilters()) 
	  saveObj = false;
	else {
	  const trigger::Keys &keys = trEv->filterKeys( moduleFilterIndex );
	  saveObj = onlineOfflineMatchingRECO(ref, &triggerObjects, &keys, dRMatch_);
	}
      }
      for (size_t f=1; f<filterNames_.size(); f++) {
      	if (filterNames_[f] != "") {
      	  unsigned int moduleFilterIndex = trEv->filterIndex(filterNames_[f]);
	  for (int i=0;i<trEv->sizeFilters(); i++) {
	    if (trEv->filterLabel(i) == filterNames_[f]) {
	      moduleFilterIndex = i;
	      break;
	    }//if (trEv->filterLabel(i) == filterNames_[f])
	  }//for (int i=0;i<trEv->sizeFilters(); i++)

      	  if (moduleFilterIndex+1 > trEv->sizeFilters()) 
      	    saveObj = false;
      	  else {
      	    const trigger::Keys &keys = trEv->filterKeys( moduleFilterIndex );
      	    if (isAND_)
      	      saveObj = (saveObj && onlineOfflineMatchingRECO(ref, &triggerObjects, &keys, dRMatch_));
      	    else		   
      	      saveObj = (saveObj || onlineOfflineMatchingRECO(ref, &triggerObjects, &keys, dRMatch_));
      	  } 
      	}
      }
    }

    //std::cout << saveObj << std::endl;
    if (saveObj)
      outColRef->push_back(ref);
  }	  

  iEvent.put(std::move(outColRef));
}

template<>
bool MiniAODTriggerCandProducer<pat::Electron, pat::TriggerObjectStandAlone>::onlineOfflineMatching(pat::ElectronRef ref, 
												    const std::vector<pat::TriggerObjectStandAlone>* triggerObjects, 
												    std::string filterLabel, float dRmin,const edm::Handle<edm::TriggerResults> & triggerBits,const edm::TriggerNames &triggerNames, edm::Event &iEvent) {
  
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
    //obj.unpackPathNames(triggerNames); 
    obj.unpackPathNames(triggerNames);
    obj.unpackFilterLabels(iEvent, *triggerBits);
    if (obj.hasFilterLabel(filterLabel)) {
      float dR = deltaR(ref->superCluster()->position(), obj.p4());
      if (dR < dRmin)
	return true;
    }
  }

  return false;
}

template <>
bool MiniAODTriggerCandProducer<pat::Photon, pat::TriggerObjectStandAlone>::onlineOfflineMatching(pat::PhotonRef ref,
												  const std::vector<pat::TriggerObjectStandAlone>* triggerObjects, 
												  std::string filterLabel, float dRmin,const edm::Handle<edm::TriggerResults> & triggerBits,const edm::TriggerNames &triggerNames, edm::Event &iEvent) {
  
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
    //obj.unpackPathNames(triggerNames); 

    obj.unpackPathNames(triggerNames);
    obj.unpackFilterLabels(iEvent, *triggerBits);
    if (obj.hasFilterLabel(filterLabel)) {
      float dR = deltaR(ref->superCluster()->position(), obj.p4());
      if (dR < dRmin)
	return true;
    }
  }

  return false;
}

template <>
bool MiniAODTriggerCandProducer<reco::RecoEcalCandidate, pat::TriggerObjectStandAlone>::onlineOfflineMatching(edm::Ref<std::vector<reco::RecoEcalCandidate>> ref,
													      const std::vector<pat::TriggerObjectStandAlone>* triggerObjects, 
													      std::string filterLabel, float dRmin,const edm::Handle<edm::TriggerResults> & triggerBits,const edm::TriggerNames &triggerNames, edm::Event &iEvent) {

  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
    //obj.unpackPathNames(triggerNames); 

    obj.unpackPathNames(triggerNames);
    obj.unpackFilterLabels(iEvent, *triggerBits);
    if (obj.hasFilterLabel(filterLabel)) {
      float dR = deltaR(ref->superCluster()->position(), obj.p4());
      if (dR < dRmin)
	return true;
    }
  }

  return false;
}

template <>
MiniAODTriggerCandProducer<pat::Electron, pat::TriggerObjectStandAlone>::MiniAODTriggerCandProducer(const edm::ParameterSet& iConfig ) :
  filterNames_(iConfig.getParameter<std::vector<std::string> >("filterNames")),
  inputs_(consumes<edm::RefVector<pat::ElectronCollection>>(iConfig.getParameter<edm::InputTag>("inputs"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  jets_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
  tightJetIdToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("tightJetId"))),
  slimmedEleToken_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("slimmedElectrons"))),
  dRMatch_(iConfig.getParameter<double>("dR")),
  isAND_(iConfig.getParameter<bool>("isAND")),
  hasJet_(iConfig.getParameter<bool>("hasJet")) {

  if(iConfig.existsAs<edm::InputTag>("triggerEvent"))
    triggerEvent_ = mayConsume<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("triggerEvent"));

  produces<edm::RefVector<pat::ElectronCollection>>();
  produces<edm::ValueMap<float>>("jetPt");
  produces<edm::ValueMap<float>>("jetEta");
  produces<edm::ValueMap<float>>("jetPhi");
}

template <>
void MiniAODTriggerCandProducer<pat::Electron, pat::TriggerObjectStandAlone>::produce(edm::Event &iEvent, const edm::EventSetup &eventSetup) {
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;

  edm::Handle<edm::RefVector<pat::ElectronCollection>> inputs;
  edm::Handle<edm::View<pat::Electron>> slimmedElectrons;

  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  iEvent.getByToken(inputs_, inputs);
  iEvent.getByToken(slimmedEleToken_, slimmedElectrons);

  // Create the output collection
  std::unique_ptr<edm::RefVector<pat::ElectronCollection>> outColRef(new edm::RefVector<pat::ElectronCollection>);
  std::vector<float> outJetPt, outJetEta, outJetPhi;
  outJetPt.reserve(slimmedElectrons->size());
  outJetEta.reserve(slimmedElectrons->size());
  outJetPhi.reserve(slimmedElectrons->size());
  std::vector<std::pair<pat::ElectronRef,float>> mapJetPt, mapJetEta, mapJetPhi;

  if (!triggerBits.isValid()) {
    LogDebug("") << "TriggerResults product not found - returning result=false!";
    return;
  }

  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerBits);
  if (triggerNamesID_ != triggerNames.parameterSetID()) {
    triggerNamesID_ = triggerNames.parameterSetID();
    init(*triggerBits, triggerNames);
  }

  std::vector<pat::JetRef> matchedJets;

  if (hasJet_) {
    edm::Handle<edm::View<pat::Jet>> jetHandle;
    iEvent.getByToken(jets_, jetHandle);

    edm::Handle<edm::ValueMap<int>> tightJetIdHandle;
    iEvent.getByToken(tightJetIdToken_, tightJetIdHandle);

    for (size_t i = 0; i < jetHandle->size(); i++) {
      pat::JetRef matchedJet;
      const auto& aJet = jetHandle->refAt(i);

      if ( (*tightJetIdHandle)[aJet]==0 )
        continue;

      if (filterNames_.size() > 0) {
        for (size_t f = 0; f < filterNames_.size(); f++) {
          for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
            obj.unpackPathNames(triggerNames);
            obj.unpackFilterLabels(iEvent, *triggerBits);

            if (obj.hasFilterLabel(filterNames_[f])) {
              float dR = deltaR(aJet->p4(), obj.p4());
              if (dR < dRMatch_) {
                matchedJet = aJet.castTo<pat::JetRef>();
                break;
              }
            }
          }

          if (matchedJet.isNonnull())
            break;
        }
      }

      if (matchedJet.isNonnull())
        matchedJets.push_back(matchedJet);
    }

    auto sortByPt = [] (const pat::JetRef& aJet, const pat::JetRef& bJet) { return aJet->pt() > bJet->pt(); };

    if ( matchedJets.size() > 1 )
      std::partial_sort(matchedJets.begin(),matchedJets.begin()+1,matchedJets.end(),sortByPt);
  }

  for (size_t i = 0 ; i < inputs->size(); i++) {
    bool saveObj = true;
    pat::ElectronRef ref = (*inputs)[i];

    if (filterNames_.size() > 0) {
      saveObj = onlineOfflineMatching(ref, triggerObjects.product(), filterNames_[0], dRMatch_,triggerBits,triggerNames,iEvent);

      for (size_t f = 1; f < filterNames_.size(); f++) {
        if (isAND_)
          saveObj = (saveObj && onlineOfflineMatching(ref, triggerObjects.product(), filterNames_[f], dRMatch_,triggerBits,triggerNames,iEvent));
        else
          saveObj = (saveObj || onlineOfflineMatching(ref, triggerObjects.product(), filterNames_[f], dRMatch_,triggerBits,triggerNames,iEvent));
      }
    }

    if (saveObj) {
      if (!matchedJets.empty()) {
        const pat::JetRef& aJet = matchedJets.front();
        mapJetPt.push_back(std::make_pair<pat::ElectronRef,float>(pat::ElectronRef(ref),aJet->pt()));
        mapJetEta.push_back(std::make_pair<pat::ElectronRef,float>(pat::ElectronRef(ref),aJet->eta()));
        mapJetPhi.push_back(std::make_pair<pat::ElectronRef,float>(pat::ElectronRef(ref),aJet->phi()));
      }

      outColRef->push_back(ref);
    }
  }

  for (size_t i = 0; i < slimmedElectrons->size(); i++) {
    pat::ElectronRef ref = slimmedElectrons->refAt(i).castTo<pat::ElectronRef>();

    float aJetPt = 0.;
    float aJetEta = 0.;
    float aJetPhi = 0.;

    for (size_t iref = 0; iref < mapJetPt.size(); iref++) {
      if ( ref==mapJetPt.at(iref).first ) {
        aJetPt = mapJetPt.at(iref).second;
        aJetEta = mapJetEta.at(iref).second;
        aJetPhi = mapJetPhi.at(iref).second;
        break;
      }
    }

    outJetPt.emplace_back(aJetPt);
    outJetEta.emplace_back(aJetEta);
    outJetPhi.emplace_back(aJetPhi);
  }

  writeValueMap(iEvent, slimmedElectrons, outJetPt, "jetPt");
  writeValueMap(iEvent, slimmedElectrons, outJetEta, "jetEta");
  writeValueMap(iEvent, slimmedElectrons, outJetPhi, "jetPhi");
  iEvent.put(std::move(outColRef));
}

typedef MiniAODTriggerCandProducer<reco::GsfElectron, trigger::TriggerObject> GsfElectronTriggerCandProducer;
DEFINE_FWK_MODULE(GsfElectronTriggerCandProducer);

typedef MiniAODTriggerCandProducer<pat::Electron, pat::TriggerObjectStandAlone> PatElectronTriggerCandProducer;
DEFINE_FWK_MODULE(PatElectronTriggerCandProducer);

typedef MiniAODTriggerCandProducer<pat::Photon, pat::TriggerObjectStandAlone> PatPhotonTriggerCandProducer;
DEFINE_FWK_MODULE(PatPhotonTriggerCandProducer);

typedef MiniAODTriggerCandProducer<reco::RecoEcalCandidate, pat::TriggerObjectStandAlone> RecoEcalCandidateTriggerCandProducer;
DEFINE_FWK_MODULE(RecoEcalCandidateTriggerCandProducer);
