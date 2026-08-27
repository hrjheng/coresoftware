#include "DFinderNtuple.h"

#include <TFile.h>
#include <TNamed.h>
#include <TTree.h>

#include <string>

bool DFinderNtuple::open(const std::string &fileName, const bool keepSameSign)
{
  m_file.reset(TFile::Open(fileName.c_str(), "RECREATE"));
  if (!m_file || m_file->IsZombie())
  {
    m_file.reset();
    return false;
  }

  m_file->SetCompressionLevel(4);
  m_file->cd();
  TNamed sameSignProvenance("DFinderConfig.keepSameSign",
                            keepSameSign ? "true" : "false");
  sameSignProvenance.Write();
  m_candidateTree = new TTree("DFinder", "DFinder keep-all candidates");
  m_eventTree = new TTree("DFinderEvent", "DFinder event bookkeeping");
  m_cutFlowTree = new TTree("DFinderCutFlow", "DFinder cumulative cut flow");
  initializeCandidateBranches();
  initializeEventBranches();
  initializeCutFlowBranches();
  return true;
}

void DFinderNtuple::fillCandidate(const DFinderCandidateRecord &candidate)
{
  m_candidate = candidate;
  m_candidateTree->Fill();
}

void DFinderNtuple::fillEvent(const DFinderEventRecord &event)
{
  m_event = event;
  m_eventTree->Fill();
}

void DFinderNtuple::close(const DFinderCutFlow &cutFlow)
{
  if (!m_file)
  {
    return;
  }

  m_cutFlow = cutFlow;
  m_cutFlowTree->Fill();
  m_file->cd();
  m_candidateTree->Write();
  m_eventTree->Write();
  m_cutFlowTree->Write();
  m_file->Close();
  m_file.reset();
}

void DFinderNtuple::initializeCandidateBranches()
{
  m_candidateTree->Branch("run", &m_candidate.run, "run/I");
  m_candidateTree->Branch("event", &m_candidate.event, "event/I");
  m_candidateTree->Branch("centrality", &m_candidate.centrality, "centrality/F");
  m_candidateTree->Branch("nTracks", &m_candidate.nTracks, "nTracks/I");
  m_candidateTree->Branch("nPV", &m_candidate.nPV, "nPV/I");
  m_candidateTree->Branch("PVx", &m_candidate.PVx, "PVx/F");
  m_candidateTree->Branch("PVy", &m_candidate.PVy, "PVy/F");
  m_candidateTree->Branch("PVz", &m_candidate.PVz, "PVz/F");
  m_candidateTree->Branch("PVchi2", &m_candidate.PVchi2, "PVchi2/F");
  m_candidateTree->Branch("PVndf", &m_candidate.PVndf, "PVndf/I");
  m_candidateTree->Branch("PVid", &m_candidate.PVid, "PVid/I");
  m_candidateTree->Branch("pvCrossing", &m_candidate.pvCrossing, "pvCrossing/I");
  m_candidateTree->Branch("daughterCrossing", &m_candidate.daughterCrossing, "daughterCrossing/I");
  m_candidateTree->Branch("pvMatched", &m_candidate.pvMatched, "pvMatched/I");
  m_candidateTree->Branch("nVertexSameCrossing", &m_candidate.nVertexSameCrossing, "nVertexSameCrossing/I");
  m_candidateTree->Branch("pvMatchedBy", &m_candidate.pvMatchedBy, "pvMatchedBy/I");
  m_candidateTree->Branch("pvNtracks", &m_candidate.pvNtracks, "pvNtracks/I");
  m_candidateTree->Branch("truncated", &m_candidate.truncated, "truncated/O");

  m_candidateTree->Branch("mass", &m_candidate.mass, "mass/F");
  m_candidateTree->Branch("unfitted_mass", &m_candidate.unfittedMass, "unfitted_mass/F");
  m_candidateTree->Branch("pt", &m_candidate.pt, "pt/F");
  m_candidateTree->Branch("unfitted_pt", &m_candidate.unfittedPt, "unfitted_pt/F");
  m_candidateTree->Branch("eta", &m_candidate.eta, "eta/F");
  m_candidateTree->Branch("phi", &m_candidate.phi, "phi/F");
  m_candidateTree->Branch("y", &m_candidate.rapidity, "y/F");
  m_candidateTree->Branch("px", &m_candidate.px, "px/F");
  m_candidateTree->Branch("py", &m_candidate.py, "py/F");
  m_candidateTree->Branch("pz", &m_candidate.pz, "pz/F");
  m_candidateTree->Branch("vtxX", &m_candidate.vtxX, "vtxX/F");
  m_candidateTree->Branch("vtxY", &m_candidate.vtxY, "vtxY/F");
  m_candidateTree->Branch("vtxZ", &m_candidate.vtxZ, "vtxZ/F");
  m_candidateTree->Branch("vtxXErr", &m_candidate.vtxXErr, "vtxXErr/F");
  m_candidateTree->Branch("vtxYErr", &m_candidate.vtxYErr, "vtxYErr/F");
  m_candidateTree->Branch("vtxZErr", &m_candidate.vtxZErr, "vtxZErr/F");
  m_candidateTree->Branch("vtxYXErr", &m_candidate.vtxYXErr, "vtxYXErr/F");
  m_candidateTree->Branch("vtxZXErr", &m_candidate.vtxZXErr, "vtxZXErr/F");
  m_candidateTree->Branch("vtxZYErr", &m_candidate.vtxZYErr, "vtxZYErr/F");
  m_candidateTree->Branch("vtxchi2", &m_candidate.vtxchi2, "vtxchi2/F");
  m_candidateTree->Branch("vtxdof", &m_candidate.vtxdof, "vtxdof/I");
  m_candidateTree->Branch("svpvDistance", &m_candidate.svpvDistance, "svpvDistance/F");
  m_candidateTree->Branch("svpvDisErr", &m_candidate.svpvDisErr, "svpvDisErr/F");
  m_candidateTree->Branch("svpvDistance_2D", &m_candidate.svpvDistance2D, "svpvDistance_2D/F");
  m_candidateTree->Branch("svpvDisErr_2D", &m_candidate.svpvDisErr2D, "svpvDisErr_2D/F");
  m_candidateTree->Branch("alpha", &m_candidate.alpha, "alpha/F");
  m_candidateTree->Branch("MaxDoca", &m_candidate.maxDoca, "MaxDoca/F");
  m_candidateTree->Branch("MaxDocaXY", &m_candidate.maxDocaXY, "MaxDocaXY/F");

  m_candidateTree->Branch("fdchi2", &m_candidate.fdchi2, "fdchi2/F");
  m_candidateTree->Branch("dira", &m_candidate.dira, "dira/F");
  m_candidateTree->Branch("dira_xy", &m_candidate.diraXY, "dira_xy/F");
  m_candidateTree->Branch("decayLength", &m_candidate.decayLength, "decayLength/F");
  m_candidateTree->Branch("decayLengthErr", &m_candidate.decayLengthErr, "decayLengthErr/F");
  m_candidateTree->Branch("decayLength_xy", &m_candidate.decayLengthXY, "decayLength_xy/F");
  m_candidateTree->Branch("decayLengthErr_xy", &m_candidate.decayLengthErrXY, "decayLengthErr_xy/F");
  m_candidateTree->Branch("decayTime", &m_candidate.decayTime, "decayTime/F");
  m_candidateTree->Branch("decayTimeErr", &m_candidate.decayTimeErr, "decayTimeErr/F");
  m_candidateTree->Branch("pseudoProperDecayTime", &m_candidate.pseudoProperDecayTime, "pseudoProperDecayTime/F");
  m_candidateTree->Branch("PV_dca", &m_candidate.PVdca, "PV_dca/F");
  m_candidateTree->Branch("PV_dca_stddev", &m_candidate.PVdcaStddev, "PV_dca_stddev/F");
  m_candidateTree->Branch("PV_dca_xy", &m_candidate.PVdcaXY, "PV_dca_xy/F");
  m_candidateTree->Branch("PV_dca_stddev_xy", &m_candidate.PVdcaStddevXY, "PV_dca_stddev_xy/F");
  m_candidateTree->Branch("massErr", &m_candidate.massErr, "massErr/F");
  m_candidateTree->Branch("ellipsoidVolume", &m_candidate.ellipsoidVolume, "ellipsoidVolume/F");

  for (size_t index = 0; index < m_candidate.daughters.size(); ++index)
  {
    const std::string prefix = "rftk" + std::to_string(index + 1);
    DFinderDaughterRecord &daughter = m_candidate.daughters[index];
    m_candidateTree->Branch((prefix + "_pt").c_str(), &daughter.pt, (prefix + "_pt/F").c_str());
    m_candidateTree->Branch((prefix + "_eta").c_str(), &daughter.eta, (prefix + "_eta/F").c_str());
    m_candidateTree->Branch((prefix + "_phi").c_str(), &daughter.phi, (prefix + "_phi/F").c_str());
    m_candidateTree->Branch((prefix + "_mass").c_str(), &daughter.mass, (prefix + "_mass/F").c_str());
    m_candidateTree->Branch((prefix + "_index").c_str(), &daughter.index, (prefix + "_index/I").c_str());
    m_candidateTree->Branch((prefix + "_parentIntermediateIndex").c_str(), &daughter.parentIntermediateIndex, (prefix + "_parentIntermediateIndex/I").c_str());
    m_candidateTree->Branch((prefix + "_MassHypo").c_str(), &daughter.massHypothesis, (prefix + "_MassHypo/F").c_str());
    m_candidateTree->Branch((prefix + "_charge").c_str(), &daughter.charge, (prefix + "_charge/I").c_str());
    m_candidateTree->Branch((prefix + "_chi2ndof").c_str(), &daughter.chi2ndof, (prefix + "_chi2ndof/F").c_str());
    m_candidateTree->Branch((prefix + "_nMvtx").c_str(), &daughter.nMvtx, (prefix + "_nMvtx/I").c_str());
    m_candidateTree->Branch((prefix + "_nIntt").c_str(), &daughter.nIntt, (prefix + "_nIntt/I").c_str());
    m_candidateTree->Branch((prefix + "_nTpc").c_str(), &daughter.nTpc, (prefix + "_nTpc/I").c_str());
    m_candidateTree->Branch((prefix + "_crossing").c_str(), &daughter.crossing, (prefix + "_crossing/I").c_str());
    m_candidateTree->Branch((prefix + "_dEdx").c_str(), &daughter.dEdx, (prefix + "_dEdx/F").c_str());
    m_candidateTree->Branch((prefix + "_ip3d").c_str(), &daughter.ip3d, (prefix + "_ip3d/F").c_str());
    m_candidateTree->Branch((prefix + "_ip3dSig").c_str(), &daughter.ip3dSig, (prefix + "_ip3dSig/F").c_str());
    m_candidateTree->Branch((prefix + "_ipxy").c_str(), &daughter.ipxy, (prefix + "_ipxy/F").c_str());
    m_candidateTree->Branch((prefix + "_ipxySig").c_str(), &daughter.ipxySig, (prefix + "_ipxySig/F").c_str());
  }

  m_candidateTree->Branch("nIntermediate", &m_candidate.nIntermediate, "nIntermediate/I");
  for (size_t index = 0; index < m_candidate.intermediates.size(); ++index)
  {
    const std::string prefix = "intermediate" + std::to_string(index + 1);
    DFinderIntermediateRecord &intermediate = m_candidate.intermediates[index];
    m_candidateTree->Branch((prefix + "_pdgId").c_str(), &intermediate.pdgId, (prefix + "_pdgId/I").c_str());
    m_candidateTree->Branch((prefix + "_parentIntermediateIndex").c_str(), &intermediate.parentIntermediateIndex, (prefix + "_parentIntermediateIndex/I").c_str());
    m_candidateTree->Branch((prefix + "_mass").c_str(), &intermediate.mass, (prefix + "_mass/F").c_str());
    m_candidateTree->Branch((prefix + "_massErr").c_str(), &intermediate.massErr, (prefix + "_massErr/F").c_str());
    m_candidateTree->Branch((prefix + "_pt").c_str(), &intermediate.pt, (prefix + "_pt/F").c_str());
    m_candidateTree->Branch((prefix + "_eta").c_str(), &intermediate.eta, (prefix + "_eta/F").c_str());
    m_candidateTree->Branch((prefix + "_phi").c_str(), &intermediate.phi, (prefix + "_phi/F").c_str());
    m_candidateTree->Branch((prefix + "_px").c_str(), &intermediate.px, (prefix + "_px/F").c_str());
    m_candidateTree->Branch((prefix + "_py").c_str(), &intermediate.py, (prefix + "_py/F").c_str());
    m_candidateTree->Branch((prefix + "_pz").c_str(), &intermediate.pz, (prefix + "_pz/F").c_str());
    m_candidateTree->Branch((prefix + "_vtxX").c_str(), &intermediate.vtxX, (prefix + "_vtxX/F").c_str());
    m_candidateTree->Branch((prefix + "_vtxY").c_str(), &intermediate.vtxY, (prefix + "_vtxY/F").c_str());
    m_candidateTree->Branch((prefix + "_vtxZ").c_str(), &intermediate.vtxZ, (prefix + "_vtxZ/F").c_str());
    m_candidateTree->Branch((prefix + "_vtxchi2").c_str(), &intermediate.vtxchi2, (prefix + "_vtxchi2/F").c_str());
    m_candidateTree->Branch((prefix + "_vtxdof").c_str(), &intermediate.vtxdof, (prefix + "_vtxdof/I").c_str());
  }

  m_candidateTree->Branch("isTruthMatched", &m_candidate.isTruthMatched, "isTruthMatched/O");
  m_candidateTree->Branch("truthPdgId", &m_candidate.truthPdgId, "truthPdgId/I");
  m_candidateTree->Branch("truthPt", &m_candidate.truthPt, "truthPt/F");
  m_candidateTree->Branch("truthY", &m_candidate.truthY, "truthY/F");
  m_candidateTree->Branch("truthDecayLength", &m_candidate.truthDecayLength, "truthDecayLength/F");
  m_candidateTree->Branch("isPrompt", &m_candidate.isPrompt, "isPrompt/O");
}

void DFinderNtuple::initializeEventBranches()
{
  m_eventTree->Branch("run", &m_event.run, "run/I");
  m_eventTree->Branch("event", &m_event.event, "event/I");
  m_eventTree->Branch("nTracks", &m_event.nTracks, "nTracks/I");
  m_eventTree->Branch("nSelectedTracks", &m_event.nSelectedTracks, "nSelectedTracks/I");
  m_eventTree->Branch("nPV", &m_event.nPV, "nPV/I");
  m_eventTree->Branch("combinations", &m_event.combinations, "combinations/l");
  m_eventTree->Branch("fitsAttempted", &m_event.fitsAttempted, "fitsAttempted/l");
  m_eventTree->Branch("fitsValid", &m_event.fitsValid, "fitsValid/l");
  m_eventTree->Branch("rowsWritten", &m_event.rowsWritten, "rowsWritten/l");
  m_eventTree->Branch("truncated", &m_event.truncated, "truncated/O");
}

void DFinderNtuple::initializeCutFlowBranches()
{
  m_cutFlowTree->Branch("events", &m_cutFlow.events, "events/l");
  m_cutFlowTree->Branch("tracksInput", &m_cutFlow.tracksInput, "tracksInput/l");
  m_cutFlowTree->Branch("tracksSelected", &m_cutFlow.tracksSelected, "tracksSelected/l");
  m_cutFlowTree->Branch("tracksCrossingUndefined", &m_cutFlow.tracksCrossingUndefined, "tracksCrossingUndefined/l");
  m_cutFlowTree->Branch("tracksCrossingZero", &m_cutFlow.tracksCrossingZero, "tracksCrossingZero/l");
  m_cutFlowTree->Branch("crossingGroups", &m_cutFlow.crossingGroups, "crossingGroups/l");
  m_cutFlowTree->Branch("combinations", &m_cutFlow.combinations, "combinations/l");
  m_cutFlowTree->Branch("osPairs", &m_cutFlow.osPairs, "osPairs/l");
  m_cutFlowTree->Branch("ssPairs", &m_cutFlow.ssPairs, "ssPairs/l");
  m_cutFlowTree->Branch("hypothesesFormed", &m_cutFlow.hypothesesFormed, "hypothesesFormed/l");
  m_cutFlowTree->Branch("massPass", &m_cutFlow.massPass, "massPass/l");
  m_cutFlowTree->Branch("kinematicPass", &m_cutFlow.kinematicPass, "kinematicPass/l");
  m_cutFlowTree->Branch("docaPass", &m_cutFlow.docaPass, "docaPass/l");
  m_cutFlowTree->Branch("propagationFailed", &m_cutFlow.propagationFailed, "propagationFailed/l");
  m_cutFlowTree->Branch("fitsAttempted", &m_cutFlow.fitsAttempted, "fitsAttempted/l");
  m_cutFlowTree->Branch("fitsValid", &m_cutFlow.fitsValid, "fitsValid/l");
  m_cutFlowTree->Branch("noMatchedVertex", &m_cutFlow.noMatchedVertex, "noMatchedVertex/l");
  m_cutFlowTree->Branch("ambiguousVertex", &m_cutFlow.ambiguousVertex, "ambiguousVertex/l");
  m_cutFlowTree->Branch("rowsWritten", &m_cutFlow.rowsWritten, "rowsWritten/l");
  m_cutFlowTree->Branch("truncatedEvents", &m_cutFlow.truncatedEvents, "truncatedEvents/l");
}
