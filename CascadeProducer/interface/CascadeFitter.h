// -*- C++ -*-
//
// Package:    V0Producer
// Class:      CascadeFitter
// 
/**\class CascadeFitter CascadeFitter.h RecoVertex/V0Producer/interface/CascadeFitter.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Brian Drell
//         Created:  Fri May 18 22:57:40 CEST 2007
// $Id: CascadeFitter.h,v 1.24 2010/08/05 22:06:39 wmtan Exp $
//
//

#ifndef RECOVERTEX__CASCADE_FITTER_H
#define RECOVERTEX__CASCADE_FITTER_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/VolumeBasedEngine/interface/VolumeBasedMagneticField.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/Math/interface/angle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
//#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
//#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"

#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <string>
#include <fstream>

class CascadeFitter {
 public:
  CascadeFitter(const edm::ParameterSet& theParams,
	   const edm::Event& iEvent, const edm::EventSetup& iSetup);
  ~CascadeFitter();

  // Switching to L. Lista's reco::Candidate infrastructure for Cascade storage
  const reco::VertexCompositeCandidateCollection& getXis() const;
  const reco::VertexCompositeCandidateCollection& getOmegas() const;
  const reco::VertexCompositeCandidateCollection& getLambdaC() const;

 private:
  // STL vector of VertexCompositeCandidate that will be filled with VertexCompositeCandidates by fitAll()
  reco::VertexCompositeCandidateCollection theXis;
  reco::VertexCompositeCandidateCollection theOmegas;
  reco::VertexCompositeCandidateCollection theLambdaCs;

  // Tracker geometry for discerning hit positions
  const TrackerGeometry* trackerGeom;

  const MagneticField* magField;

  edm::InputTag recoAlg;
  edm::InputTag vtxAlg;
  edm::InputTag v0CollName;
  bool doXis;
  bool doOmegas;
  bool doLambdaCs;

  // Cuts
  double tkChi2Cut;
  int    tkNhitsCut;
  double tkPtCut;
  double xiChi2Cut;
  double xiRVtxCut;
  double xiRVtxSigCut;
  double xiLVtxCut;
  double xiLVtxSigCut;
  double xiCollinCut;
  double lambdaCMassCut;
  double xiMassCut;
  double omegaMassCut;
  double dauTransImpactSigCut;
  double dauLongImpactSigCut;

  std::vector<reco::TrackBase::TrackQuality> qualities;

  edm::InputTag vtxFitter;

  // Helper method that does the actual fitting using the KalmanVertexFitter
  void fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  double findCascadeMassError(const GlobalPoint &vtxPos, std::vector<reco::TransientTrack> dauTracks);
};

#endif
