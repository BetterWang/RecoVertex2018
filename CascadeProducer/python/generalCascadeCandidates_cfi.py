import FWCore.ParameterSet.Config as cms

generalCascadeCandidates = cms.EDProducer("CascadeProducer",
                                     
    # InputTag that tells which TrackCollection to use for vertexing
    trackRecoAlgorithm = cms.InputTag('generalTracks'),
    vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices'),
    v0CollName = cms.InputTag('generalV0Candidates','Lambda'), 

    # These bools decide whether or not to reconstruct
    #  specific V0 particles
    selectLambdaCs = cms.bool(True),
    selectXis = cms.bool(True),
    selectOmegas = cms.bool(True),

    # Select tracks using TrackBase::TrackQuality.
    # Select ALL tracks by leaving this vstring empty, which
    #   is equivalent to using 'loose'
    #trackQualities = cms.vstring('highPurity', 'goodIterative'),
    trackQualities = cms.vstring('loose'),
                                     
    # The next parameters are cut values.
    # All distances are in cm, all energies in GeV, as usual.

    # --Track quality/compatibility cuts--
    #   Normalized track Chi2 <
    tkChi2Cut = cms.double(5.0),
    #   Number of valid hits on track >=
    tkNhitsCut = cms.int32(6),
    #   Number of valid hits on track >=
    tkPtCut = cms.double(0.0),
    #   Track impact parameter significance >
    dauTransImpactSigCut = cms.double(2.),
    dauLongImpactSigCut = cms.double(2.),

    # cuts for Xi 
    xiVtxChi2Cut = cms.double(7.0),
    xiCollinearityCut = cms.double(-2.0),
    xiRVtxCut = cms.double(0.0),
    xiLVtxCut = cms.double(0.0),
    xiVtxSignificance2DCut = cms.double(15.0),
    xiVtxSignificance3DCut = cms.double(0.0),

    # V0 mass window, Candidate mass must be within these values of
    # the PDG mass to be stored in the collection
    lambdaCMassCut = cms.double(0.3),
    xiMassCut = cms.double(0.10),
    omegaMassCut = cms.double(0.10)
)
