// -*- C++ -*-
//
// Package:    D0Producer
// Class:      D0Producer
// 
/**\class D0Producer D0Producer.h RecoVertex/D0Producer/interface/D0Producer.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wei Li 
//
//

#ifndef RECOVERTEX__D0_PRODUCER_H
#define RECOVERTEX__D0_PRODUCER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"

#include "RecoVertex/D0Producer/interface/D0Fitter.h"

class D0Producer : public edm::EDProducer {
public:
  explicit D0Producer(const edm::ParameterSet&);
  ~D0Producer();

private:
  //virtual void beginJob() ;
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  D0Fitter theVees; 
//  edm::ParameterSet theParams;
};

#endif
