// -*- C++ -*-
//
// Package:    V0Producer
// Class:      CascadeProducer
// 
/**\class CascadeProducer CascadeProducer.h RecoVertex/V0Producer/interface/CascadeProducer.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Brian Drell
//         Created:  Fri May 18 22:57:40 CEST 2007
// $Id: CascadeProducer.h,v 1.8 2009/12/18 20:45:08 wmtan Exp $
//
//

#ifndef RECOVERTEX__CASCADE_PRODUCER_H
#define RECOVERTEX__CASCADE_PRODUCER_H

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

#include "RecoVertex/CascadeProducer/interface/CascadeFitter.h"

class CascadeProducer : public edm::EDProducer {
public:
  explicit CascadeProducer(const edm::ParameterSet&);
  ~CascadeProducer();

private:
  //virtual void beginJob() ;
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  edm::ParameterSet theParams;
      
};

#endif
