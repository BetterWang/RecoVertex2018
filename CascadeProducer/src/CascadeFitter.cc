// -*- C++ -*-
//
// Package:    V0Producer
// Class:      CascadeFitter
// 
/**\class CascadeFitter CascadeFitter.cc RecoVertex/V0Producer/src/CascadeFitter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Brian Drell
//         Created:  Fri May 18 22:57:40 CEST 2007
// $Id: CascadeFitter.cc,v 1.56 2011/12/23 08:13:39 innocent Exp $
//
//

#include "RecoVertex/CascadeProducer/interface/CascadeFitter.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include <typeinfo>
#include <memory>
#include <vector>

const ParticleMass piMass = 0.13957018;
const double piMassSquared = piMass*piMass;
const ParticleMass protonMass = 0.938272013;
const double protonMassSquared = protonMass*protonMass;
const ParticleMass kaonMass = 0.493677;
const double kaonMassSquared = kaonMass*kaonMass;
const ParticleMass kShortMass = 0.497614;
const ParticleMass lambdaMass = 1.115683;
const ParticleMass d0Mass = 1.86484;
const ParticleMass lambdaCMass = 2.28646;
const ParticleMass xiMass = 1.32171;
const ParticleMass omegaMass = 1.67245;
float piMass_sigma = piMass*1.e-6;
float protonMass_sigma = protonMass*1.e-6;
float kaonMass_sigma = kaonMass*1.e-6;
float kShortMass_sigma = kShortMass*1.e-6;
float lambdaMass_sigma = lambdaMass*1.e-6;
float d0Mass_sigma = d0Mass*1.e-6;
float lambdaCMass_sigma = lambdaCMass*1.e-6;
float xiMass_sigma = xiMass*1.e-6;
float omegaMass_sigma = omegaMass*1.e-6;

// Constructor and (empty) destructor
CascadeFitter::CascadeFitter(const edm::ParameterSet& theParameters,
		   const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using std::string;

  // Get the track reco algorithm from the ParameterSet
  recoAlg = theParameters.getParameter<edm::InputTag>("trackRecoAlgorithm");
  vtxAlg  = theParameters.getParameter<edm::InputTag>("vertexRecoAlgorithm");
  v0CollName = theParameters.getParameter<edm::InputTag>("v0CollName");

  // ------> Initialize parameters from PSet. ALL TRACKED, so no defaults.
  // First set bits to do various things:
  //  -decide whether to use the KVF track smoother, and whether to store those
  //     tracks in the reco::Vertex
  //  -whether to reconstruct Xi 
  doXis = theParameters.getParameter<bool>(string("selectXis"));
  //  -whether to reconstruct Omega
  doOmegas = theParameters.getParameter<bool>(string("selectOmegas"));
  //  -whether to reconstruct LambdaC 
  doLambdaCs = theParameters.getParameter<bool>(string("selectLambdaCs"));

  // Second, initialize post-fit cuts
  tkChi2Cut = theParameters.getParameter<double>(string("tkChi2Cut"));
  tkNhitsCut = theParameters.getParameter<int>(string("tkNhitsCut"));
  tkPtCut = theParameters.getParameter<double>(string("tkPtCut"));
  // Sames cuts for Xi
  xiChi2Cut = theParameters.getParameter<double>(string("xiVtxChi2Cut"));
  xiRVtxCut = theParameters.getParameter<double>(string("xiRVtxCut"));
  xiRVtxSigCut = theParameters.getParameter<double>(string("xiVtxSignificance2DCut"));
  xiLVtxCut = theParameters.getParameter<double>(string("xiLVtxCut"));
  xiLVtxSigCut = theParameters.getParameter<double>(string("xiVtxSignificance3DCut"));
  xiCollinCut = theParameters.getParameter<double>(string("xiCollinearityCut"));
  // other common cuts
  lambdaCMassCut = theParameters.getParameter<double>(string("lambdaCMassCut"));
  xiMassCut = theParameters.getParameter<double>(string("xiMassCut"));
  omegaMassCut = theParameters.getParameter<double>(string("omegaMassCut"));
  dauTransImpactSigCut = theParameters.getParameter<double>(string("dauTransImpactSigCut"));
  dauLongImpactSigCut = theParameters.getParameter<double>(string("dauLongImpactSigCut"));
  std::vector<std::string> qual = theParameters.getParameter<std::vector<std::string> >("trackQualities");
  for (unsigned int ndx = 0; ndx < qual.size(); ndx++) {
    qualities.push_back(reco::TrackBase::qualityByName(qual[ndx]));
  }

  fitAll(iEvent, iSetup);
}

CascadeFitter::~CascadeFitter() {
}

// Method containing the algorithm for vertex reconstruction
void CascadeFitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using std::vector;
  using std::cout;
  using std::endl;
  using namespace reco;
  using namespace edm;

  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

  // Create std::vectors for Tracks and TrackRefs (required for
  //  passing to the KalmanVertexFitter)
  std::vector<TrackRef> theTrackRefs;
  std::vector<TransientTrack> theTransTracks;

  // Handles for tracks, B-field, and tracker geometry
  Handle<reco::TrackCollection> theTrackHandle;
  Handle<reco::VertexCollection> theVertexHandle;
  Handle<reco::BeamSpot> theBeamSpotHandle;
  ESHandle<MagneticField> bFieldHandle;
  ESHandle<TrackerGeometry> trackerGeomHandle;
  ESHandle<GlobalTrackingGeometry> globTkGeomHandle;
  //cout << "Check 0" << endl;

  // Get the tracks, vertices from the event, and get the B-field record
  //  from the EventSetup
  iEvent.getByLabel(recoAlg, theTrackHandle);
  iEvent.getByLabel(vtxAlg,  theVertexHandle);
  iEvent.getByLabel(std::string("offlineBeamSpot"), theBeamSpotHandle);
  if( !theTrackHandle->size() ) return;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
  iSetup.get<TrackerDigiGeometryRecord>().get(trackerGeomHandle);
  iSetup.get<GlobalTrackingGeometryRecord>().get(globTkGeomHandle);

  trackerGeom = trackerGeomHandle.product();
  magField = bFieldHandle.product();

  edm::Handle<reco::VertexCompositeCandidateCollection> theLambdas;
  iEvent.getByLabel(v0CollName,theLambdas);
  if(!theLambdas.isValid()) return;

  if( !doXis && !doOmegas && !doLambdaCs ) return;

  bool isVtxPV = 0;
  double xVtx=-99999.0;
  double yVtx=-99999.0;
  double zVtx=-99999.0;
  double xVtxError=-999.0;
  double yVtxError=-999.0;
  double zVtxError=-999.0;
  const reco::VertexCollection vtxCollection = *(theVertexHandle.product());
  reco::VertexCollection::const_iterator vtxPrimary = vtxCollection.begin();
  if(vtxCollection.size()>0 && !vtxPrimary->isFake() && vtxPrimary->tracksSize()>=2)
  {
    isVtxPV = 1;
    xVtx = vtxPrimary->x();
    yVtx = vtxPrimary->y();
    zVtx = vtxPrimary->z();
    xVtxError = vtxPrimary->xError();
    yVtxError = vtxPrimary->yError();
    zVtxError = vtxPrimary->zError();
  }
  else {
    isVtxPV = 0;
    xVtx = theBeamSpotHandle->position().x();
    yVtx = theBeamSpotHandle->position().y();
    zVtx = 0.0;
    xVtxError = theBeamSpotHandle->BeamWidthX();
    yVtxError = theBeamSpotHandle->BeamWidthY();
    zVtxError = 0.0;
  }

  // Fill vectors of TransientTracks and TrackRefs after applying preselection cuts.
  for(unsigned int indx = 0; indx < theTrackHandle->size(); indx++) {
    TrackRef tmpRef( theTrackHandle, indx );
    bool quality_ok = true;
    if (qualities.size()!=0) {
      quality_ok = false;
      for (unsigned int ndx_ = 0; ndx_ < qualities.size(); ndx_++) {
	if (tmpRef->quality(qualities[ndx_])){
	  quality_ok = true;
	  break;          
	}
      }
    }
    if( !quality_ok ) continue;

    if( tmpRef->normalizedChi2() < tkChi2Cut &&
        tmpRef->numberOfValidHits() >= tkNhitsCut &&
        tmpRef->pt() > tkPtCut ) {
      TransientTrack tmpTk( *tmpRef, &(*bFieldHandle), globTkGeomHandle );

      math::XYZPoint bestvtx(xVtx,yVtx,zVtx);
      double dzvtx = tmpRef->dz(bestvtx);
      double dxyvtx = tmpRef->dxy(bestvtx);      
      double dzerror = sqrt(tmpRef->dzError()*tmpRef->dzError()+zVtxError*zVtxError);
      double dxyerror = sqrt(tmpRef->d0Error()*tmpRef->d0Error()+xVtxError*yVtxError);

      double dauTransImpactSig = dzvtx/dzerror;
      double dauLongImpactSig = dxyvtx/dxyerror;

      if( fabs(dauTransImpactSig) > dauTransImpactSigCut && fabs(dauLongImpactSig) > dauLongImpactSigCut ) {
        theTrackRefs.push_back( tmpRef );
        theTransTracks.push_back( tmpTk );
      }
    }
  }

    for(unsigned it=0; it<theLambdas->size(); ++it){
        
      const reco::VertexCompositeCandidate & theLambda = (*theLambdas)[it];

      // check lambda mass to be within +- 20MeV...
      float massWindow = 0.020;
      if(theLambda.mass() > lambdaMass + massWindow || theLambda.mass() < lambdaMass - massWindow) continue;
      // after finding values also cut on xi vertex probability > 0.0001 and batPion Dxyz/sigma > 0.5

      int lamPDG = theLambda.pdgId();
      bool lamIsParticle = (lamPDG > 0);

      //get daughter tracks from V0 candidate
      vector<RecoChargedCandidate> v0daughters;
      vector<TrackRef> theDaughterTracks;

      //check momentum and add pion first, proton second
      if (theLambda.daughter(0)->momentum().mag2() < theLambda.daughter(1)->momentum().mag2()){
         v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
                                  (theLambda.daughter(0))) );
         v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
                                  (theLambda.daughter(1))) );
       } else {
         v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
                                  (theLambda.daughter(1))) );
         v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
                                  (theLambda.daughter(0))) );
       }
       for(unsigned int j = 0; j < v0daughters.size(); ++j) {
         theDaughterTracks.push_back(v0daughters[j].track());
       }

       vector<TransientTrack> tracksForKalmanFit;
       for (unsigned int ndx = 0; ndx < theDaughterTracks.size(); ndx++) {
           tracksForKalmanFit.push_back(TransientTrack(theDaughterTracks[ndx], &(*bFieldHandle)));
       }

       vector<double> lamDauMasses;
       if (lamIsParticle && tracksForKalmanFit[0].charge() > 0) {
           lamDauMasses.push_back(protonMass);
           lamDauMasses.push_back(piMass);
       } else if (!lamIsParticle && tracksForKalmanFit[0].charge() < 0) {
           lamDauMasses.push_back(protonMass);
           lamDauMasses.push_back(piMass);
       } else {
           lamDauMasses.push_back(piMass);
           lamDauMasses.push_back(protonMass);
       }

      // Loop over tracks and vertex good charged track pairs
      for(unsigned int trdx = 0; trdx < theTrackRefs.size(); trdx++) {

         if ( theTrackRefs[trdx].isNull() ) continue;

         math::XYZPoint bestvtx(xVtx,yVtx,zVtx);
         double dzvtx = theTrackRefs[trdx]->dz(bestvtx);
         double dxyvtx = theTrackRefs[trdx]->dxy(bestvtx);
         if(fabs(dzvtx)>50 || fabs(dxyvtx)>50) continue;

         bool match = false;
         // check if the pion is used in any v0 candidate in the collection and flag it
         for(unsigned int j = 0; j < theDaughterTracks.size(); ++j) {
           if (theTrackRefs[trdx]->charge() == theDaughterTracks[j]->charge() &&
               theTrackRefs[trdx]->momentum() == theDaughterTracks[j]->momentum() ) match = true;
           if (match) break;
         }
         if (match) continue; // Track is already used in making the V0

         // check if pion is in *any* good V0
         // Placeholder
       
         if(lamIsParticle && theTrackRefs[trdx]->charge()>0) continue;
         if(!lamIsParticle && theTrackRefs[trdx]->charge()<0) continue;

         TransientTrack pionTT(theDaughterTracks[0], &(*bFieldHandle) );
         TransientTrack protonTT(theDaughterTracks[1], &(*bFieldHandle) );
         TransientTrack batPionTT(theTrackRefs[trdx], &(*bFieldHandle) );

         if (!batPionTT.isValid()) continue;
         if (!pionTT.isValid()) continue;
         if (!protonTT.isValid()) continue;

         //Creating a KinematicParticleFactory
         KinematicParticleFactoryFromTransientTrack pFactory;

         //initial chi2 and ndf before kinematic fits.
         float chi = 0.;
         float ndf = 0.;
         vector<RefCountedKinematicParticle> lambdaParticles;
         lambdaParticles.push_back(pFactory.particle(pionTT,piMass,chi,ndf,piMass_sigma));
         lambdaParticles.push_back(pFactory.particle(protonTT,protonMass,chi,ndf,protonMass_sigma));

         KinematicParticleVertexFitter fitter;
         RefCountedKinematicTree lambdaVertexFitTree;
         lambdaVertexFitTree = fitter.fit(lambdaParticles);
         if (!lambdaVertexFitTree->isValid()) continue;
  
         lambdaVertexFitTree->movePointerToTheTop();

         RefCountedKinematicParticle lambda_vFit = lambdaVertexFitTree->currentParticle();
         RefCountedKinematicVertex lambda_vFit_vertex = lambdaVertexFitTree->currentDecayVertex();

         lambdaVertexFitTree->movePointerToTheFirstChild();
         RefCountedKinematicParticle lambdaPi1 = lambdaVertexFitTree->currentParticle();
         lambdaVertexFitTree->movePointerToTheNextChild();
         RefCountedKinematicParticle lambdaPi2 = lambdaVertexFitTree->currentParticle();

//         KinematicParameters lambdaPiKP = lambdaPi1->currentState().kinematicParameters();
//         KinematicParameters lambdaPKP  = lambdaPi2->currentState().kinematicParameters();

         // now apply Lambda mass constraint
         KinematicParticleFitter csFitterLambda;
         KinematicConstraint * lam_c = new MassKinematicConstraint(lambdaMass,lambdaMass_sigma);
         // add mass constraint to the lambda fit to do a constrained fit:  

         lambdaVertexFitTree->movePointerToTheTop();
         lambdaVertexFitTree = csFitterLambda.fit(lam_c,lambdaVertexFitTree);
         if (!lambdaVertexFitTree->isValid()) continue;
         lambdaVertexFitTree->movePointerToTheTop();
         RefCountedKinematicParticle lambda_vFit_withMC = lambdaVertexFitTree->currentParticle();

         if (!lambda_vFit_withMC->currentState().isValid()) continue; 

         if(doXis || doLambdaCs)
         {
           vector<RefCountedKinematicParticle> xiFitParticles;

           xiFitParticles.push_back(pFactory.particle(batPionTT,piMass,chi,ndf,piMass_sigma));
           xiFitParticles.push_back(lambda_vFit_withMC);

           //fit Xi 
           RefCountedKinematicTree xiFitTree = fitter.fit(xiFitParticles);
           if (!xiFitTree->isValid()) continue;

           xiFitTree->movePointerToTheTop();
           RefCountedKinematicParticle xiCand = xiFitTree->currentParticle();
           if (!xiCand->currentState().isValid()) continue;

           RefCountedKinematicVertex xiDecayVertex = xiFitTree->currentDecayVertex();
           if (!xiDecayVertex->vertexIsValid()) continue;

           if ( xiDecayVertex->chiSquared()<0 || xiDecayVertex->chiSquared()>1000 ) continue;

           float xiC2Prob =
              ChiSquaredProbability((double)(xiDecayVertex->chiSquared()),(double)(xiDecayVertex->degreesOfFreedom()));
           if (xiC2Prob < 0.0001) continue;

           if ( xiCand->currentState().mass() > 3 ) continue;

           // get children from final Xi fit
           xiFitTree->movePointerToTheFirstChild();
           RefCountedKinematicParticle batPionCand = xiFitTree->currentParticle();
           xiFitTree->movePointerToTheNextChild();
           RefCountedKinematicParticle lambdaCandMC = xiFitTree->currentParticle();
         
           if(!batPionCand->currentState().isValid() || !lambdaCandMC->currentState().isValid()) continue;

           // get batchlor pion and V0 parameters from Xi fit
           KinematicParameters batPionKP = batPionCand->currentState().kinematicParameters();
           KinematicParameters lambdaCandKP = lambdaCandMC->currentState().kinematicParameters();

           GlobalVector xiTotalP = GlobalVector (xiCand->currentState().globalMomentum().x(),
                                                 xiCand->currentState().globalMomentum().y(),
                                                 xiCand->currentState().globalMomentum().z());

           GlobalVector batPionTotalP = GlobalVector(batPionKP.momentum().x(),batPionKP.momentum().y(),batPionKP.momentum().z());
           GlobalVector lambdaTotalP = GlobalVector(lambdaCandKP.momentum().x(),lambdaCandKP.momentum().y(),lambdaCandKP.momentum().z());

           double batPionTotalE = sqrt( batPionTotalP.mag2() + piMassSquared );
           double lambdaTotalE = sqrt( lambdaTotalP.mag2() + lambdaMass*lambdaMass );
           double xiTotalE = batPionTotalE + lambdaTotalE;

           const Particle::LorentzVector xiP4(xiTotalP.x(), xiTotalP.y(), xiTotalP.z(), xiTotalE);

           Particle::Point xiVtx((*xiDecayVertex).position().x(), (*xiDecayVertex).position().y(), (*xiDecayVertex).position().z());
           std::vector<double> xiVtxEVec;
           xiVtxEVec.push_back( lambda_vFit_vertex->error().cxx() );
           xiVtxEVec.push_back( lambda_vFit_vertex->error().cyx() );
           xiVtxEVec.push_back( lambda_vFit_vertex->error().cyy() );
           xiVtxEVec.push_back( lambda_vFit_vertex->error().czx() );
           xiVtxEVec.push_back( lambda_vFit_vertex->error().czy() );
           xiVtxEVec.push_back( lambda_vFit_vertex->error().czz() );
           SMatrixSym3D xiVtxCovMatrix(xiVtxEVec.begin(), xiVtxEVec.end());
           const Vertex::CovarianceMatrix xiVtxCov(xiVtxCovMatrix);
           double xiVtxChi2(xiDecayVertex->chiSquared());
           double xiVtxNdof(xiDecayVertex->degreesOfFreedom());
           double xiNormalizedChi2 = xiVtxChi2/xiVtxNdof;

           double xiRVtxMag = 99999.0;
           double xiLVtxMag = 99999.0;
           double xiSigmaRvtxMag = 999.0;
           double xiSigmaLvtxMag = 999.0;
           double xiV0Angle = -100.0;

           GlobalVector xiV0LineOfFlight = GlobalVector (xiVtx.x() - xVtx,
                                                         xiVtx.y() - yVtx,
                                                         xiVtx.z() - zVtx);

           SMatrixSym3D xiTotalCov;
           if(isVtxPV) xiTotalCov = xiVtxCovMatrix + vtxPrimary->covariance();
           else xiTotalCov = xiVtxCovMatrix + theBeamSpotHandle->rotatedCovariance3D();

           SVector3 distanceVector3D(xiV0LineOfFlight.x(), xiV0LineOfFlight.y(), xiV0LineOfFlight.z());
           SVector3 distanceVector2D(xiV0LineOfFlight.x(), xiV0LineOfFlight.y(), 0.0);

           xiV0Angle = angle(xiV0LineOfFlight.x(), xiV0LineOfFlight.y(), xiV0LineOfFlight.z(),
                             xiTotalP.x(), xiTotalP.y(), xiTotalP.z());
           xiLVtxMag = xiV0LineOfFlight.mag();
           xiRVtxMag = xiV0LineOfFlight.perp();
           xiSigmaLvtxMag = sqrt(ROOT::Math::Similarity(xiTotalCov, distanceVector3D)) / xiLVtxMag;
           xiSigmaRvtxMag = sqrt(ROOT::Math::Similarity(xiTotalCov, distanceVector2D)) / xiRVtxMag;

           if( xiNormalizedChi2 > xiChi2Cut ||
               xiRVtxMag < xiRVtxCut ||
               xiRVtxMag / xiSigmaRvtxMag < xiRVtxSigCut ||
               xiLVtxMag < xiLVtxCut ||
               xiLVtxMag / xiSigmaLvtxMag < xiLVtxSigCut ||
               cos(xiV0Angle) < xiCollinCut
           ) continue;

           // Cuts finished, Create the VertexCompositeCandidate object for Xi that will be stored in the Event
           RecoChargedCandidate PionCand(theTrackRefs[trdx]->charge(), Particle::LorentzVector(batPionTotalP.x(),
                                                   batPionTotalP.y(), batPionTotalP.z(),
                                                   batPionTotalE), xiVtx);
           PionCand.setTrack(theTrackRefs[trdx]);

           AddFourMomenta addp4;
           // Store the daughter Candidates in the VertexCompositeCandidates 
           // if they pass mass cuts
           if(doXis)
           {
             VertexCompositeCandidate* theXi = 0;
             theXi = new VertexCompositeCandidate(0, xiP4, xiVtx, xiVtxCov, xiVtxChi2, xiVtxNdof);
             theXi->addDaughter(theLambda);
             theXi->addDaughter(PionCand);
             if(lamIsParticle) theXi->setPdgId(3312);
             else theXi->setPdgId(-3312);
             addp4.set( *theXi );
             if( theXi->mass() < xiMass + xiMassCut &&
                 theXi->mass() > xiMass - xiMassCut ) theXis.push_back( *theXi );
             if(theXi) delete theXi;
             theXi = 0;
           }
           if(doLambdaCs)
           {
             VertexCompositeCandidate* theLambdaC = 0;
             theLambdaC = new VertexCompositeCandidate(0, xiP4, xiVtx, xiVtxCov, xiVtxChi2, xiVtxNdof);
             theLambdaC->addDaughter(theLambda);
             theLambdaC->addDaughter(PionCand);
             if(lamIsParticle) theLambdaC->setPdgId(4122);
             else theLambdaC->setPdgId(-4122);
             addp4.set( *theLambdaC );
             if( theLambdaC->mass() < lambdaCMass + lambdaCMassCut &&
                 theLambdaC->mass() > lambdaCMass - lambdaCMassCut ) theLambdaCs.push_back( *theLambdaC );
             if(theLambdaC) delete theLambdaC;
             theLambdaC = 0;
           }
        }

         if(doOmegas)
         {
           vector<RefCountedKinematicParticle> omegaFitParticles;

           omegaFitParticles.push_back(pFactory.particle(batPionTT,kaonMass,chi,ndf,kaonMass_sigma));
           omegaFitParticles.push_back(lambda_vFit_withMC);

           //fit Omega 
           RefCountedKinematicTree omegaFitTree = fitter.fit(omegaFitParticles);
           if (!omegaFitTree->isValid()) continue;

           omegaFitTree->movePointerToTheTop();
           RefCountedKinematicParticle omegaCand = omegaFitTree->currentParticle();
           if (!omegaCand->currentState().isValid()) continue;

           RefCountedKinematicVertex omegaDecayVertex = omegaFitTree->currentDecayVertex();
           if (!omegaDecayVertex->vertexIsValid()) continue;

           if ( omegaDecayVertex->chiSquared()<0 || omegaDecayVertex->chiSquared()>1000 ) continue;

           float omegaC2Prob =
              ChiSquaredProbability((double)(omegaDecayVertex->chiSquared()),(double)(omegaDecayVertex->degreesOfFreedom()));
           if (omegaC2Prob < 0.0001) continue;

           if ( omegaCand->currentState().mass() > 3 ) continue;

           // get children from final Omega fit
           omegaFitTree->movePointerToTheFirstChild();
           RefCountedKinematicParticle batPionCand = omegaFitTree->currentParticle();
           omegaFitTree->movePointerToTheNextChild();
           RefCountedKinematicParticle lambdaCandMC = omegaFitTree->currentParticle();
         
           if(!batPionCand->currentState().isValid() || !lambdaCandMC->currentState().isValid()) continue;

           // get batchlor pion and V0 parameters from Xi fit
           KinematicParameters batPionKP = batPionCand->currentState().kinematicParameters();
           KinematicParameters lambdaCandKP = lambdaCandMC->currentState().kinematicParameters();

           GlobalVector omegaTotalP = GlobalVector (omegaCand->currentState().globalMomentum().x(),
                                                    omegaCand->currentState().globalMomentum().y(),
                                                    omegaCand->currentState().globalMomentum().z());

           GlobalVector batPionTotalP = GlobalVector(batPionKP.momentum().x(),batPionKP.momentum().y(),batPionKP.momentum().z());
           GlobalVector lambdaTotalP = GlobalVector(lambdaCandKP.momentum().x(),lambdaCandKP.momentum().y(),lambdaCandKP.momentum().z());


           double batPionTotalE = sqrt( batPionTotalP.mag2() + kaonMassSquared );
           double lambdaTotalE = sqrt( lambdaTotalP.mag2() + lambdaMass*lambdaMass );
           double omegaTotalE = batPionTotalE + lambdaTotalE;

           const Particle::LorentzVector omegaP4(omegaTotalP.x(), omegaTotalP.y(), omegaTotalP.z(), omegaTotalE);

           Particle::Point omegaVtx((*omegaDecayVertex).position().x(), (*omegaDecayVertex).position().y(), (*omegaDecayVertex).position().z());
           std::vector<double> omegaVtxEVec;
           omegaVtxEVec.push_back( lambda_vFit_vertex->error().cxx() );
           omegaVtxEVec.push_back( lambda_vFit_vertex->error().cyx() );
           omegaVtxEVec.push_back( lambda_vFit_vertex->error().cyy() );
           omegaVtxEVec.push_back( lambda_vFit_vertex->error().czx() );
           omegaVtxEVec.push_back( lambda_vFit_vertex->error().czy() );
           omegaVtxEVec.push_back( lambda_vFit_vertex->error().czz() );
           SMatrixSym3D omegaVtxCovMatrix(omegaVtxEVec.begin(), omegaVtxEVec.end());
           const Vertex::CovarianceMatrix omegaVtxCov(omegaVtxCovMatrix);
           double omegaVtxChi2(omegaDecayVertex->chiSquared());
           double omegaVtxNdof(omegaDecayVertex->degreesOfFreedom());
           double omegaNormalizedChi2 = omegaVtxChi2/omegaVtxNdof;

           double omegaRVtxMag = 99999.0;
           double omegaLVtxMag = 99999.0;
           double omegaSigmaRvtxMag = 999.0;
           double omegaSigmaLvtxMag = 999.0;
           double omegaV0Angle = -100.0;

           GlobalVector omegaV0LineOfFlight = GlobalVector (omegaVtx.x() - xVtx,
                                                            omegaVtx.y() - yVtx,
                                                            omegaVtx.z() - zVtx);

           SMatrixSym3D omegaTotalCov;
           if(isVtxPV) omegaTotalCov = omegaVtxCovMatrix + vtxPrimary->covariance();
           else omegaTotalCov = omegaVtxCovMatrix + theBeamSpotHandle->rotatedCovariance3D();

           SVector3 distanceVector3D(omegaV0LineOfFlight.x(), omegaV0LineOfFlight.y(), omegaV0LineOfFlight.z());
           SVector3 distanceVector2D(omegaV0LineOfFlight.x(), omegaV0LineOfFlight.y(), 0.0);

           omegaV0Angle = angle(omegaV0LineOfFlight.x(), omegaV0LineOfFlight.y(), omegaV0LineOfFlight.z(),
                             omegaTotalP.x(), omegaTotalP.y(), omegaTotalP.z());
           omegaLVtxMag = omegaV0LineOfFlight.mag();
           omegaRVtxMag = omegaV0LineOfFlight.perp();
           omegaSigmaLvtxMag = sqrt(ROOT::Math::Similarity(omegaTotalCov, distanceVector3D)) / omegaLVtxMag;
           omegaSigmaRvtxMag = sqrt(ROOT::Math::Similarity(omegaTotalCov, distanceVector2D)) / omegaRVtxMag;

           if( omegaNormalizedChi2 > xiChi2Cut ||
               omegaRVtxMag < xiRVtxCut ||
               omegaRVtxMag / omegaSigmaRvtxMag < xiRVtxSigCut ||
               omegaLVtxMag < xiLVtxCut ||
               omegaLVtxMag / omegaSigmaLvtxMag < xiLVtxSigCut ||
               cos(omegaV0Angle) < xiCollinCut
           ) continue;

           // Cuts finished, Create the VertexCompositeCandidate object for Omega that will be stored in the Event
           VertexCompositeCandidate* theOmega = 0;
           theOmega = new VertexCompositeCandidate(0, omegaP4, omegaVtx, omegaVtxCov, omegaVtxChi2, omegaVtxNdof);

           RecoChargedCandidate PionCand(theTrackRefs[trdx]->charge(), Particle::LorentzVector(batPionTotalP.x(),
                                                   batPionTotalP.y(), batPionTotalP.z(),
                                                   batPionTotalE), omegaVtx);
           PionCand.setTrack(theTrackRefs[trdx]);

           AddFourMomenta addp4;
           // Store the daughter Candidates in the VertexCompositeCandidates 
           // if they pass mass cuts
           theOmega->addDaughter(theLambda);
           theOmega->addDaughter(PionCand);
           if(lamIsParticle) theOmega->setPdgId(3334);
           else theOmega->setPdgId(-3334);
           addp4.set( *theOmega );
           if( theOmega->mass() < omegaMass + omegaMassCut &&
               theOmega->mass() > omegaMass - omegaMassCut ) theOmegas.push_back( *theOmega );
           if(theOmega) delete theOmega;
           theOmega = 0;
	 }
      }
    }
}

// Get methods
const reco::VertexCompositeCandidateCollection& CascadeFitter::getXis() const {
  return theXis;
}

const reco::VertexCompositeCandidateCollection& CascadeFitter::getOmegas() const {
  return theOmegas;
}

const reco::VertexCompositeCandidateCollection& CascadeFitter::getLambdaC() const {
  return theLambdaCs;
}

// Experimental
double CascadeFitter::findCascadeMassError(const GlobalPoint &vtxPos, std::vector<reco::TransientTrack> dauTracks) { 
  return -1.;
}

/*
double CascadeFitter::findCascadeMassError(const GlobalPoint &vtxPos, std::vector<reco::TransientTrack> dauTracks) {
  // Returns -99999. if trajectory states fail at vertex position

  // Load positive track trajectory at vertex into vector, then negative track
  std::vector<TrajectoryStateClosestToPoint> sortedTrajStatesAtVtx;
  for( unsigned int ndx = 0; ndx < dauTracks.size(); ndx++ ) {
    if( dauTracks[ndx].trajectoryStateClosestToPoint(vtxPos).isValid() ) {
      std::cout << "From TSCP: " 
		<< dauTracks[ndx].trajectoryStateClosestToPoint(vtxPos).perigeeParameters().transverseCurvature()
		<< "; From Track: " << dauTracks[ndx].track().qoverp() << std::endl;
    }
    if( sortedTrajStatesAtVtx.size() == 0 ) {
      if( dauTracks[ndx].charge() > 0 ) {
	sortedTrajStatesAtVtx.push_back( dauTracks[ndx].trajectoryStateClosestToPoint(vtxPos) );
      }
      else {
	sortedTrajStatesAtVtx.push_back( dauTracks[ndx].trajectoryStateClosestToPoint(vtxPos) );
      }
    }
  }
  std::vector<PerigeeTrajectoryParameters> param;
  std::vector<PerigeeTrajectoryError> paramError;
  std::vector<GlobalVector> momenta;

  for( unsigned int ndx2 = 0; ndx2 < sortedTrajStatesAtVtx.size(); ndx2++ ) {
    if( sortedTrajStatesAtVtx[ndx2].isValid() ) {
      param.push_back( sortedTrajStatesAtVtx[ndx2].perigeeParameters() );
      paramError.push_back( sortedTrajStatesAtVtx[ndx2].perigeeError() );
      momenta.push_back( sortedTrajStatesAtVtx[ndx2].momentum() );
    }
    else return -99999.;
  }
  return 0;
}
*/



