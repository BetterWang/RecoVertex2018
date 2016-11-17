import FWCore.ParameterSet.Config as cms

generalD0Candidates = cms.EDProducer("D0Producer",
                                     
    # InputTag that tells which TrackCollection to use for vertexing
    trackRecoAlgorithm = cms.InputTag('generalTracks'),
    vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices'),

    trackQualities = cms.vstring('loose'),
                                     
    tkChi2Cut = cms.double(5.0),
    tkNhitsCut = cms.int32(3),
    tkPtCut = cms.double(0.0),

    #   Track impact parameter significance >
    dauTransImpactSigCut = cms.double(1.),
    dauLongImpactSigCut = cms.double(1.),

    #   PCA distance between tracks <
    tkDCACut = cms.double(1.),
    vtxChi2Cut = cms.double(7.0),
    collinearityCut = cms.double(-2.0),
    rVtxCut = cms.double(0.0),
    lVtxCut = cms.double(0.0),
    vtxSignificance2DCut = cms.double(0.0),
    vtxSignificance3DCut = cms.double(0.0),
    d0MassCut = cms.double(0.2),
)
