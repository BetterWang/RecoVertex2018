// -*- C++ -*-
//
// Package:    V0Producer
//
// Class:      CascadeProducer
// 
/**\class CascadeProducer CascadeProducer.cc RecoVertex/V0Producer/src/CascadeProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Brian Drell
//         Created:  Fri May 18 22:57:40 CEST 2007
// $Id: CascadeProducer.cc,v 1.12 2010/02/20 21:02:02 wmtan Exp $
//
//


// system include files
#include <memory>

#include "RecoVertex/CascadeProducer/interface/CascadeProducer.h"

// Constructor
CascadeProducer::CascadeProducer(const edm::ParameterSet& iConfig) :
  theParams(iConfig) {

  // Trying this with Candidates instead of the simple reco::Vertex
  produces< reco::VertexCompositeCandidateCollection >("Xi");
  produces< reco::VertexCompositeCandidateCollection >("Omega");
  produces< reco::VertexCompositeCandidateCollection >("LambdaC");

}

// (Empty) Destructor
CascadeProducer::~CascadeProducer() {
}


//
// Methods
//

// Producer Method
void CascadeProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;

   // Create CascadeFitter object which reconstructs the vertices and creates
   //  (and contains) collections of Kshorts, Lambda0s
   CascadeFitter theVees(theParams, iEvent, iSetup);

   // Create unique_ptr for each collection to be stored in the Event
   std::unique_ptr< reco::VertexCompositeCandidateCollection >
     xiCandidates( new reco::VertexCompositeCandidateCollection );
   xiCandidates->reserve( theVees.getXis().size() );

   std::unique_ptr< reco::VertexCompositeCandidateCollection >
     omegaCandidates( new reco::VertexCompositeCandidateCollection );
   omegaCandidates->reserve( theVees.getOmegas().size() );

   std::unique_ptr< reco::VertexCompositeCandidateCollection >
     lambdaCCandidates( new reco::VertexCompositeCandidateCollection );
   lambdaCCandidates->reserve( theVees.getLambdaC().size() );

   std::copy( theVees.getXis().begin(),
              theVees.getXis().end(),
              std::back_inserter(*xiCandidates) );
   std::copy( theVees.getOmegas().begin(),
              theVees.getOmegas().end(),
              std::back_inserter(*omegaCandidates) );
   std::copy( theVees.getLambdaC().begin(),
              theVees.getLambdaC().end(),
              std::back_inserter(*lambdaCCandidates) );

   // Write the collections to the Event
   iEvent.put( std::move(xiCandidates), std::string("Xi") );
   iEvent.put( std::move(omegaCandidates), std::string("Omega") );
   iEvent.put( std::move(lambdaCCandidates), std::string("LambdaC") );
}


void CascadeProducer::beginJob() {
}


void CascadeProducer::endJob() {
}

//define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(CascadeProducer);
