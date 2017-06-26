// -*- C++ -*-
//
// Package:    V0Producer
//
// Class:      V0Producer
// 
/**\class V0Producer V0Producer.cc MyProducers/V0Producer/src/V0Producer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Brian Drell
//         Created:  Fri May 18 22:57:40 CEST 2007
// $Id: V0Producer.cc,v 1.12 2010/02/20 21:02:02 wmtan Exp $
//
//


// system include files
#include <memory>

#include "RecoVertex/V0Producer/interface/V0Producer.h"

// Constructor
V0Producer::V0Producer(const edm::ParameterSet& iConfig) :
 theVees(iConfig, consumesCollector())
//  theParams(iConfig) {
{
   // Registering V0 Collections
  //produces<reco::VertexCollection>("Kshort");
  //produces<reco::VertexCollection>("Lambda");
  //produces<reco::VertexCollection>("LambdaBar");

  // Trying this with Candidates instead of the simple reco::Vertex
  produces< reco::VertexCompositeCandidateCollection >("Kshort");
  produces< reco::VertexCompositeCandidateCollection >("Lambda");
  produces< reco::VertexCompositeCandidateCollection >("Xi");
  produces< reco::VertexCompositeCandidateCollection >("Omega");
  produces< reco::VertexCompositeCandidateCollection >("D0");
  produces< reco::VertexCompositeCandidateCollection >("DS");
  produces< reco::VertexCompositeCandidateCollection >("DPM");
  produces< reco::VertexCompositeCandidateCollection >("LambdaCToLamPi");
  produces< reco::VertexCompositeCandidateCollection >("LambdaCToKsP");
  //produces< reco::VertexCompositeCandidateCollection >("LambdaBar");

}

// (Empty) Destructor
V0Producer::~V0Producer() {
}


//
// Methods
//

// Producer Method
void V0Producer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;

   // Create V0Fitter object which reconstructs the vertices and creates
   //  (and contains) collections of Kshorts, Lambda0s
//   V0Fitter theVees(theParams, iEvent, iSetup);

   theVees.fitAll(iEvent, iSetup);

   // Create auto_ptr for each collection to be stored in the Event
   std::auto_ptr< reco::VertexCompositeCandidateCollection > 
     kShortCandidates( new reco::VertexCompositeCandidateCollection );
   kShortCandidates->reserve( theVees.getKshorts().size() ); 

   std::auto_ptr< reco::VertexCompositeCandidateCollection >
     lambdaCandidates( new reco::VertexCompositeCandidateCollection );
   lambdaCandidates->reserve( theVees.getLambdas().size() );

   std::auto_ptr< reco::VertexCompositeCandidateCollection >
     xiCandidates( new reco::VertexCompositeCandidateCollection );
   xiCandidates->reserve( theVees.getXis().size() );

   std::auto_ptr< reco::VertexCompositeCandidateCollection >
     omegaCandidates( new reco::VertexCompositeCandidateCollection );
   omegaCandidates->reserve( theVees.getOmegas().size() );

   std::auto_ptr< reco::VertexCompositeCandidateCollection >
     d0Candidates( new reco::VertexCompositeCandidateCollection );
   d0Candidates->reserve( theVees.getD0().size() );

   std::auto_ptr< reco::VertexCompositeCandidateCollection >
     dsCandidates( new reco::VertexCompositeCandidateCollection );
   dsCandidates->reserve( theVees.getDS().size() );

   std::auto_ptr< reco::VertexCompositeCandidateCollection >
     dpmCandidates( new reco::VertexCompositeCandidateCollection );
   dpmCandidates->reserve( theVees.getDPM().size() );

   std::auto_ptr< reco::VertexCompositeCandidateCollection >
     lambdaCCandidates1( new reco::VertexCompositeCandidateCollection );
   lambdaCCandidates1->reserve( theVees.getLambdaCToLamPi().size() );

   std::auto_ptr< reco::VertexCompositeCandidateCollection >
     lambdaCCandidates2( new reco::VertexCompositeCandidateCollection );
   lambdaCCandidates2->reserve( theVees.getLambdaCToKsP().size() );

   std::copy( theVees.getKshorts().begin(),
	      theVees.getKshorts().end(),
	      std::back_inserter(*kShortCandidates) );
   std::copy( theVees.getLambdas().begin(),
	      theVees.getLambdas().end(),
	      std::back_inserter(*lambdaCandidates) );
   std::copy( theVees.getXis().begin(),
              theVees.getXis().end(),
              std::back_inserter(*xiCandidates) );
   std::copy( theVees.getOmegas().begin(),
              theVees.getOmegas().end(),
              std::back_inserter(*omegaCandidates) );
   std::copy( theVees.getD0().begin(),
              theVees.getD0().end(),
              std::back_inserter(*d0Candidates) );
   std::copy( theVees.getDS().begin(),
              theVees.getDS().end(),
              std::back_inserter(*dsCandidates) );
   std::copy( theVees.getDPM().begin(),
              theVees.getDPM().end(),
              std::back_inserter(*dpmCandidates) );
   std::copy( theVees.getLambdaCToLamPi().begin(),
              theVees.getLambdaCToLamPi().end(),
              std::back_inserter(*lambdaCCandidates1) );
   std::copy( theVees.getLambdaCToKsP().begin(),
              theVees.getLambdaCToKsP().end(),
              std::back_inserter(*lambdaCCandidates2) );

   // Write the collections to the Event
   iEvent.put( kShortCandidates, std::string("Kshort") );
   iEvent.put( lambdaCandidates, std::string("Lambda") );
   iEvent.put( xiCandidates, std::string("Xi") );
   iEvent.put( omegaCandidates, std::string("Omega") );
   iEvent.put( d0Candidates, std::string("D0") );
   iEvent.put( dsCandidates, std::string("DS") );
   iEvent.put( dpmCandidates, std::string("DPM") );
   iEvent.put( lambdaCCandidates1, std::string("LambdaCToLamPi") );
   iEvent.put( lambdaCCandidates2, std::string("LambdaCToKsP") );

   theVees.resetAll();
}


//void V0Producer::beginJob() {
void V0Producer::beginJob() {
}


void V0Producer::endJob() {
}

//define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(V0Producer);
//DEFINE_FWK_MODULE(V0finder);
