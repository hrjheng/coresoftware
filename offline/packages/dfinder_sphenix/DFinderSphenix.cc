#include "DFinderSphenix.h"
#include "DFinderFun4All.h"

#include <centrality/CentralityInfo.h>
#include <ffaobjects/EventHeader.h>
#include <ffamodules/CDBInterface.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>
#include <phool/getClass.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxPHG4ParticleMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackreco/ActsPropagator.h>

#include <kfparticle_sphenix/KFParticle_truthAndDetTools.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <Acts/Definitions/Units.hpp>

#include <TEntryList.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

DFinderSphenix::DFinderSphenix(const std::string &name,
                               const DFinderConfig &config)
  : SubsysReco(name)
  , m_config(config)
{
}

SubsysReco *makeDFinderSphenix(const std::string &name,
                               const DFinderConfig &config)
{
  return new DFinderSphenix(name, config);
}

int DFinderSphenix::Init(PHCompositeNode * /*topNode*/)
{
  std::string parseError;
  if (!m_combinatorics.parseDecayDescriptor(m_config.decayDescriptor, m_decay, parseError))
  {
    std::cerr << Name() << ": cannot parse decay descriptor '"
              << m_config.decayDescriptor << "': " << parseError << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (m_config.massWindowLo > m_config.massWindowHi ||
      m_config.maxCandidatesPerEvent == 0)
  {
    std::cerr << Name() << ": invalid mass window or zero candidate limit" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (!m_config.intermediateMassWindows.empty() &&
      m_config.intermediateMassWindows.size() != m_decay.nIntermediateStates)
  {
    std::cerr << Name() << ": configured "
              << m_config.intermediateMassWindows.size()
              << " intermediate mass windows for "
              << m_decay.nIntermediateStates << " intermediate states" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  for (const std::pair<float, float> &window : m_config.intermediateMassWindows)
  {
    if (window.first > window.second)
    {
      std::cerr << Name() << ": invalid intermediate mass window" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }

  m_pidPermutations = m_combinatorics.makePidPermutations(m_decay);
  if (m_pidPermutations.empty())
  {
    std::cerr << Name() << ": no distinct PID hypotheses were constructed" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (!m_ntuple.open(m_config.outputFileName, m_config.keepSameSign))
  {
    std::cerr << Name() << ": failed to create output file "
              << m_config.outputFileName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_trk_map_node_name = m_config.trackMapNodeName;
  m_vtx_map_node_name = m_config.vertexMapNodeName;
  m_use_mbd_vertex = m_config.useMbdVertex;
  m_dont_use_global_vertex = !m_config.useGlobalVertex;
  return Fun4AllReturnCodes::EVENT_OK;
}

int DFinderSphenix::InitRun(PHCompositeNode *topNode)
{
  const float minimumMass = m_decay.motherMass + m_config.massWindowLo;
  const float maximumMass = m_decay.motherMass + m_config.massWindowHi;
  std::cout << Name() << ": mass window: "
            << minimumMass << " - " << maximumMass << " GeV" << std::endl;
  if (m_config.maxPairDocaXY >= m_config.maxPairDoca)
  {
    std::cerr << Name() << ": maxPairDocaXY = " << m_config.maxPairDocaXY
              << " cm cannot reject independently of maxPairDoca = "
              << m_config.maxPairDoca << " cm" << std::endl;
  }
  if (m_config.keepSameSign)
  {
    std::cout << Name()
              << ": same-sign pairs enabled; output volume may approximately double. "
              << "Check truncatedEvents before using the reference sample."
              << std::endl;
  }

  if (!findNode::getClass<SvtxTrackMap>(topNode, m_config.trackMapNodeName))
  {
    std::cerr << Name() << ": required track node "
              << m_config.trackMapNodeName << " is missing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (m_config.vertexAssociation == DFinderConfig::VertexAssociation::MatchedCrossing &&
      m_config.useMbdVertex)
  {
    std::cerr << Name()
              << ": MatchedCrossing vertex association requires a vertex source "
              << "with a beam crossing; MbdVertex does not provide one"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (m_config.useMbdVertex)
  {
    if (!findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap"))
    {
      std::cerr << Name() << ": required MbdVertexMap node is missing" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }
  else
  {
    m_svtxVertexMap = findNode::getClass<SvtxVertexMap>(topNode, m_config.vertexMapNodeName);
    if (!m_svtxVertexMap)
    {
      std::cerr << Name() << ": required vertex node "
                << m_config.vertexMapNodeName << " is missing" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }

  if (m_config.useGlobalVertex &&
      !findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap"))
  {
    std::cerr << Name() << ": required GlobalVertexMap node is missing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (m_config.doPropagateToSV)
  {
    if (!m_svtxVertexMap)
    {
      m_svtxVertexMap = findNode::getClass<SvtxVertexMap>(topNode, m_config.vertexMapNodeName);
    }
    m_actsGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
    if (!m_svtxVertexMap || !m_actsGeometry)
    {
      std::cerr << Name()
                << ": ActsGeometry and an SvtxVertexMap are required for SV propagation"
                << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }

  std::string fieldFileName = m_config.magneticFieldFile;
  if (!std::filesystem::exists(fieldFileName))
  {
    fieldFileName = CDBInterface::instance()->getUrl(fieldFileName);
  }

  TFile fieldFile(fieldFileName.c_str(), "READ");
  TTree *fieldMap = nullptr;
  fieldFile.GetObject("fieldmap", fieldMap);
  if (fieldFile.IsZombie() || !fieldMap)
  {
    std::cerr << Name() << ": cannot read fieldmap from " << fieldFileName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  constexpr unsigned int maxFieldScanRadiusCm = 2;
  constexpr unsigned int maxFieldScanZCm = 10;
  constexpr double quarterTurn = M_PI / 2.;
  float fieldZTesla = 0.F;
  for (unsigned int z = 0; z <= maxFieldScanZCm && fieldZTesla == 0.F; ++z)
  {
    for (unsigned int radius = 0; radius <= maxFieldScanRadiusCm && fieldZTesla == 0.F; ++radius)
    {
      const unsigned int nPhi = radius == 0 ? 1 : 4;
      for (unsigned int phiIndex = 0; phiIndex < nPhi && fieldZTesla == 0.F; ++phiIndex)
      {
        const double x = radius * std::cos(phiIndex * quarterTurn);
        const double y = radius * std::sin(phiIndex * quarterTurn);
        const std::string selection =
            "x == " + std::to_string(x) +
            " && y == " + std::to_string(y) +
            " && z == " + std::to_string(z);
        fieldMap->Draw(">>dfinderFieldEntries", selection.c_str(), "entrylist");
        auto *entries = dynamic_cast<TEntryList *>(gDirectory->Get("dfinderFieldEntries"));
        if (!entries)
        {
          continue;
        }
        const Long64_t entry = entries->GetEntry(0);
        if (entry < 0)
        {
          continue;
        }
        fieldMap->GetEntry(entry);
        TLeaf *fieldLeaf = fieldMap->GetLeaf("bz");
        if (fieldLeaf)
        {
          fieldZTesla = fieldLeaf->GetValue();
        }
      }
    }
  }

  if (fieldZTesla == 0.F)
  {
    std::cerr << Name() << ": no nonzero Bz sample found near the origin" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  constexpr double teslaToKilogauss = 10.;
  KFParticle::SetField(fieldZTesla * teslaToKilogauss);
  if (Verbosity() > 0)
  {
    std::cout << Name() << ": using Bz = " << fieldZTesla
              << " T from " << fieldFileName << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int DFinderSphenix::process_event(PHCompositeNode *topNode)
{
  auto *trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_config.trackMapNodeName);
  if (!trackMap)
  {
    std::cerr << Name() << ": required track node disappeared" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  std::vector<KFParticle> primaryVertices =
      makeAllPrimaryVertices(topNode, m_config.vertexMapNodeName);

  constexpr short int undefinedCrossing = std::numeric_limits<short int>::max();
  constexpr short int crossingZero = 0;
  constexpr int pvMatchRuleUnavailable = -1;
  constexpr int pvMatchedUnique = 0;
  constexpr int pvMatchedByMembership = 1;
  constexpr int pvMatchedByNearestZ = 2;
  std::vector<SvtxVertex *> primaryVertexSources(primaryVertices.size(), nullptr);
  std::vector<short int> primaryVertexCrossings(primaryVertices.size(), undefinedCrossing);
  if (!m_config.useMbdVertex && !m_config.useGlobalVertex)
  {
    for (size_t index = 0; index < primaryVertices.size(); ++index)
    {
      SvtxVertex *vertex = m_svtxVertexMap->get(primaryVertices[index].Id());
      primaryVertexSources[index] = vertex;
      if (vertex)
      {
        primaryVertexCrossings[index] = vertex->get_beam_crossing();
      }
    }
  }
  else if (!m_config.useMbdVertex)
  {
    auto *globalVertexMap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    size_t primaryVertexIndex = 0;
    for (auto &globalVertexEntry : *globalVertexMap)
    {
      GlobalVertex *globalVertex = globalVertexEntry.second;
      auto vertexIterator = globalVertex->find_vertexes(GlobalVertex::SVTX);
      if (vertexIterator == globalVertex->end_vertexes())
      {
        continue;
      }
      for (const Vertex *vertexReference : vertexIterator->second)
      {
        if (primaryVertexIndex >= primaryVertices.size())
        {
          break;
        }
        SvtxVertex *vertex = m_svtxVertexMap->get(vertexReference->get_id());
        primaryVertexSources[primaryVertexIndex] = vertex;
        if (vertex)
        {
          primaryVertexCrossings[primaryVertexIndex] = vertex->get_beam_crossing();
        }
        ++primaryVertexIndex;
      }
    }
  }

  int run = -1;
  int event = static_cast<int>(m_cutFlow.events);
  if (auto *eventHeader = findNode::getClass<EventHeader>(topNode, "EventHeader"))
  {
    run = eventHeader->get_RunNumber();
    event = eventHeader->get_EvtSequence();
  }

  float centrality = -99.F;
  if (auto *centralityInfo = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo"))
  {
    if (centralityInfo->has_centile(CentralityInfo::PROP::mbd_NS))
    {
      centrality = centralityInfo->get_centile(CentralityInfo::PROP::mbd_NS);
    }
    else if (centralityInfo->has_centile(CentralityInfo::PROP::bimp))
    {
      centrality = centralityInfo->get_centile(CentralityInfo::PROP::bimp);
    }
  }

  const unsigned long long startCombinations = m_cutFlow.combinations;
  const unsigned long long startFitsAttempted = m_cutFlow.fitsAttempted;
  const unsigned long long startFitsValid = m_cutFlow.fitsValid;
  const unsigned long long startRowsWritten = m_cutFlow.rowsWritten;
  ++m_cutFlow.events;
  m_cutFlow.tracksInput += trackMap->size();

  auto *truthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  auto *recoTruthMap = findNode::getClass<SvtxPHG4ParticleMap>(topNode, "SvtxPHG4ParticleMap");
  KFParticle_truthAndDetTools truthTools;

  std::vector<TrackRecord> selectedTracks;
  const bool groupByCrossing =
      m_config.vertexAssociation == DFinderConfig::VertexAssociation::MatchedCrossing &&
      m_config.requireSameDaughterCrossing;
  selectedTracks.reserve(trackMap->size());
  for (auto &trackEntry : *trackMap)
  {
    SvtxTrack *track = trackEntry.second;
    if (!track)
    {
      continue;
    }
    const short int crossing = track->get_crossing();
    if (crossing == undefinedCrossing)
    {
      ++m_cutFlow.tracksCrossingUndefined;
      if (groupByCrossing)
      {
        continue;
      }
    }
    else if (crossing == crossingZero)
    {
      ++m_cutFlow.tracksCrossingZero;
    }
    if (m_config.requireCrossingZero && track->get_crossing() != crossingZero)
    {
      continue;
    }

    int nMvtx = 0;
    int nIntt = 0;
    int nTpc = 0;
    for (auto stateIterator = track->begin_states();
         stateIterator != track->end_states();
         ++stateIterator)
    {
      SvtxTrackState *state = stateIterator->second;
      if (!state || state->get_pathlength() == 0)
      {
        continue;
      }

      const uint8_t detectorId = TrkrDefs::getTrkrId(state->get_cluskey());
      if (detectorId == TrkrDefs::mvtxId)
      {
        ++nMvtx;
      }
      else if (detectorId == TrkrDefs::inttId)
      {
        ++nIntt;
      }
      else if (detectorId == TrkrDefs::tpcId)
      {
        ++nTpc;
      }
    }

    const float chi2ndof =
        track->get_ndf() > 0
            ? track->get_chisq() / static_cast<float>(track->get_ndf())
            : std::numeric_limits<float>::infinity();
    if (!std::isfinite(track->get_pt()) ||
        !std::isfinite(track->get_eta()) ||
        !std::isfinite(chi2ndof) ||
        track->get_pt() <= m_config.trackMinPt ||
        std::abs(track->get_eta()) >= m_config.trackMaxEta ||
        chi2ndof >= m_config.trackMaxChi2ndof ||
        nMvtx < m_config.minMvtxStates ||
        nIntt < m_config.minInttStates ||
        nTpc < m_config.minTpcStates)
    {
      continue;
    }

    m_dst_track = track;
    KFParticle particle = makeParticle(topNode);
    particle.SetId(trackEntry.first);

    TrackRecord record;
    record.track = track;
    record.particle = particle;
    record.nMvtx = nMvtx;
    record.nIntt = nIntt;
    record.nTpc = nTpc;
    record.dEdx = get_dEdx(topNode, particle);
    if (truthInfo && recoTruthMap)
    {
      record.truth = truthTools.getTruthTrack(track, topNode);
    }
    selectedTracks.push_back(record);
  }
  m_cutFlow.tracksSelected += selectedTracks.size();

  const size_t nProngs = m_decay.prongs.size();
  std::map<short int, std::vector<size_t>> tracksByCrossing;
  for (size_t index = 0; index < selectedTracks.size(); ++index)
  {
    const short int crossing = selectedTracks[index].track->get_crossing();
    if (crossing != undefinedCrossing)
    {
      tracksByCrossing[crossing].push_back(index);
    }
  }

  for (const auto &crossingEntry : tracksByCrossing)
  {
    if (crossingEntry.second.size() >= DFinderConfig::minimumDaughters)
    {
      ++m_cutFlow.crossingGroups;
    }
  }

  std::vector<std::vector<size_t>> trackGroups;
  if (groupByCrossing)
  {
    for (const auto &crossingEntry : tracksByCrossing)
    {
      if (crossingEntry.second.size() >= nProngs)
      {
        trackGroups.push_back(crossingEntry.second);
      }
    }
  }
  else
  {
    trackGroups.emplace_back();
    trackGroups.back().reserve(selectedTracks.size());
    for (size_t index = 0; index < selectedTracks.size(); ++index)
    {
      trackGroups.back().push_back(index);
    }
  }

  std::vector<DFinderCandidateRecord> outputRecords;
  outputRecords.reserve(std::min<size_t>(m_config.maxCandidatesPerEvent, 4096));
  bool truncated = false;

  constexpr size_t daughtersPerPair = 2;
  const bool canFitCandidates =
      m_config.vertexAssociation == DFinderConfig::VertexAssociation::MatchedCrossing ||
      !primaryVertices.empty();
  for (const std::vector<size_t> &trackGroup : trackGroups)
  {
    if (trackGroup.size() < nProngs || !canFitCandidates)
    {
      continue;
    }

    std::vector<size_t> combinationPositions(nProngs);
    std::vector<size_t> trackIndices(nProngs);
    for (size_t index = 0; index < nProngs; ++index)
    {
      combinationPositions[index] = index;
      trackIndices[index] = trackGroup[index];
    }

    bool haveCombination = true;
    while (haveCombination && !truncated)
    {
      for (size_t index = 0; index < nProngs; ++index)
      {
        trackIndices[index] = trackGroup[combinationPositions[index]];
      }
      ++m_cutFlow.combinations;
      const bool isTrackPair = nProngs == daughtersPerPair;
      int qProduct = 0;
      if (isTrackPair)
      {
        qProduct = selectedTracks[trackIndices[0]].track->get_charge() *
                   selectedTracks[trackIndices[1]].track->get_charge();
      }
      const bool sameSignPair = isTrackPair && qProduct > 0;
      const bool combinationAccepted =
          !isTrackPair || (qProduct != 0 && (!sameSignPair || m_config.keepSameSign));
      if (isTrackPair && combinationAccepted)
      {
        if (sameSignPair)
        {
          ++m_cutFlow.ssPairs;
        }
        else
        {
          ++m_cutFlow.osPairs;
        }
      }

      std::set<std::vector<int>> sameSignMassAssignments;
      for (const DFinderHypothesis &hypothesis : m_pidPermutations)
      {
        if (!combinationAccepted)
        {
          break;
        }
        if (sameSignPair)
        {
          std::vector<int> massAssignment;
          massAssignment.reserve(nProngs);
          for (const DFinderProng &prong : hypothesis.prongs)
          {
            massAssignment.push_back(std::abs(prong.pdgId));
          }
          if (!sameSignMassAssignments.insert(std::move(massAssignment)).second)
          {
            continue;
          }
        }
        else
        {
          bool chargeMatches = true;
          for (size_t index = 0; index < nProngs; ++index)
          {
            if (selectedTracks[trackIndices[index]].track->get_charge() != hypothesis.prongs[index].charge)
            {
              chargeMatches = false;
              break;
            }
          }
          if (!chargeMatches)
          {
            continue;
          }
        }
        ++m_cutFlow.hypothesesFormed;

        double unfittedEnergy = 0.;
        double unfittedPx = 0.;
        double unfittedPy = 0.;
        double unfittedPz = 0.;
        std::vector<KFParticle> daughters;
        std::vector<SvtxTrack *> daughterTracks;
        daughters.reserve(nProngs);
        daughterTracks.reserve(nProngs);
        for (size_t index = 0; index < nProngs; ++index)
        {
          const TrackRecord &trackRecord = selectedTracks[trackIndices[index]];
          const double px = trackRecord.track->get_px();
          const double py = trackRecord.track->get_py();
          const double pz = trackRecord.track->get_pz();
          const double momentumSquared = px * px + py * py + pz * pz;
          unfittedEnergy += std::sqrt(momentumSquared +
                                      hypothesis.prongs[index].mass * hypothesis.prongs[index].mass);
          unfittedPx += px;
          unfittedPy += py;
          unfittedPz += pz;
          daughters.push_back(trackRecord.particle);
          daughterTracks.push_back(trackRecord.track);
        }

        const double unfittedMomentumSquared =
            unfittedPx * unfittedPx + unfittedPy * unfittedPy + unfittedPz * unfittedPz;
        const double unfittedMassSquared =
            unfittedEnergy * unfittedEnergy - unfittedMomentumSquared;
        if (unfittedMassSquared < 0.)
        {
          continue;
        }
        const float unfittedMass = std::sqrt(unfittedMassSquared);
        const float minimumMass = m_decay.motherMass + m_config.massWindowLo;
        const float maximumMass = m_decay.motherMass + m_config.massWindowHi;
        if (unfittedMass < minimumMass || unfittedMass > maximumMass)
        {
          continue;
        }
        ++m_cutFlow.massPass;

        const float unfittedPt = std::hypot(unfittedPx, unfittedPy);
        float unfittedRapidity = std::numeric_limits<float>::infinity();
        if (unfittedEnergy > std::abs(unfittedPz))
        {
          unfittedRapidity =
              0.5F * std::log((unfittedEnergy + unfittedPz) /
                              (unfittedEnergy - unfittedPz));
        }
        if (unfittedPt < m_config.motherMinPt ||
            !std::isfinite(unfittedRapidity) ||
            std::abs(unfittedRapidity) > m_config.motherMaxRapidity)
        {
          continue;
        }
        ++m_cutFlow.kinematicPass;

        float maximumDoca = 0.F;
        float maximumDocaXY = 0.F;
        size_t firstDocaTrack = 0;
        size_t secondDocaTrack = 1;
        for (size_t first = 0; first < nProngs; ++first)
        {
          for (size_t second = first + 1; second < nProngs; ++second)
          {
            const float doca =
                std::abs(daughters[first].GetDistanceFromParticle(daughters[second]));
            const float docaXY =
                std::abs(daughters[first].GetDistanceFromParticleXY(daughters[second]));
            if (doca > maximumDoca)
            {
              maximumDoca = doca;
              firstDocaTrack = first;
              secondDocaTrack = second;
            }
            maximumDocaXY = std::max(maximumDocaXY, docaXY);
          }
        }
        if (!std::isfinite(maximumDoca) ||
            !std::isfinite(maximumDocaXY) ||
            maximumDoca > m_config.maxPairDoca ||
            maximumDocaXY > m_config.maxPairDocaXY)
        {
          continue;
        }
        ++m_cutFlow.docaPass;

        if (m_config.doPropagateToSV)
        {
          std::vector<KFParticle> propagated;
          if (!propagateToSecondaryVertex(
                  daughters, daughterTracks, firstDocaTrack, secondDocaTrack, propagated))
          {
            ++m_cutFlow.propagationFailed;
            continue;
          }
          daughters = std::move(propagated);
        }

        ++m_cutFlow.fitsAttempted;
        DFinderFitResult fitResult =
            m_fitter.fit(daughters,
                         m_decay,
                         hypothesis,
                         m_config.constrainIntermediateMasses,
                         m_config.intermediateMassWindows);
        if (!fitResult.valid)
        {
          continue;
        }
        ++m_cutFlow.fitsValid;

        std::vector<PHG4Particle *> truthByProng(nProngs, nullptr);
        for (size_t index = 0; index < nProngs; ++index)
        {
          truthByProng[hypothesis.prongs[index].prongIndex] =
              selectedTracks[trackIndices[index]].truth;
        }
        PHG4Particle *truthMother =
            truthInfo
                ? matchTruthDecay(
                      m_decay.mother,
                      hypothesis.chargeConjugated,
                      truthByProng,
                      truthInfo)
                : nullptr;
        const bool isTruthMatched = truthMother != nullptr;

        short int daughterCrossing = daughterTracks.front()->get_crossing();
        for (const SvtxTrack *daughterTrack : daughterTracks)
        {
          if (daughterTrack->get_crossing() != daughterCrossing)
          {
            daughterCrossing = undefinedCrossing;
            break;
          }
        }

        std::vector<size_t> matchingVertexIndices;
        if (daughterCrossing != undefinedCrossing)
        {
          for (size_t index = 0; index < primaryVertexCrossings.size(); ++index)
          {
            if (primaryVertexCrossings[index] == daughterCrossing)
            {
              matchingVertexIndices.push_back(index);
            }
          }
        }

        std::vector<size_t> outputVertexIndices;
        int selectedPvMatchedBy = pvMatchRuleUnavailable;
        if (m_config.vertexAssociation == DFinderConfig::VertexAssociation::AllVertices)
        {
          outputVertexIndices.reserve(primaryVertices.size());
          for (size_t index = 0; index < primaryVertices.size(); ++index)
          {
            outputVertexIndices.push_back(index);
          }
        }
        else if (matchingVertexIndices.empty())
        {
          ++m_cutFlow.noMatchedVertex;
          if (!m_config.keepUnmatchedCandidates)
          {
            continue;
          }
          outputVertexIndices.push_back(primaryVertices.size());
        }
        else
        {
          size_t selectedVertexIndex = matchingVertexIndices.front();
          if (matchingVertexIndices.size() == 1)
          {
            selectedPvMatchedBy = pvMatchedUnique;
          }
          else
          {
            ++m_cutFlow.ambiguousVertex;
            bool matchedByMembership = false;
            for (const size_t vertexIndex : matchingVertexIndices)
            {
              SvtxVertex *vertex = primaryVertexSources[vertexIndex];
              bool containsAllDaughters = vertex != nullptr;
              for (const SvtxTrack *daughterTrack : daughterTracks)
              {
                if (!vertex ||
                    vertex->find_track(daughterTrack->get_id()) == vertex->end_tracks())
                {
                  containsAllDaughters = false;
                  break;
                }
              }
              if (containsAllDaughters)
              {
                selectedVertexIndex = vertexIndex;
                selectedPvMatchedBy = pvMatchedByMembership;
                matchedByMembership = true;
                break;
              }
            }
            if (!matchedByMembership)
            {
              selectedPvMatchedBy = pvMatchedByNearestZ;
              float minimumZDistance = std::numeric_limits<float>::max();
              for (const size_t vertexIndex : matchingVertexIndices)
              {
                const float zDistance =
                    std::abs(primaryVertices[vertexIndex].GetZ() - fitResult.mother.GetZ());
                if (zDistance < minimumZDistance)
                {
                  minimumZDistance = zDistance;
                  selectedVertexIndex = vertexIndex;
                }
              }
            }
          }
          outputVertexIndices.push_back(selectedVertexIndex);
        }

        for (const size_t primaryVertexIndex : outputVertexIndices)
        {
          if (outputRecords.size() >= m_config.maxCandidatesPerEvent)
          {
            truncated = true;
            break;
          }

          const bool hasPrimaryVertex = primaryVertexIndex < primaryVertices.size();
          const KFParticle *primaryVertex =
              hasPrimaryVertex ? &primaryVertices[primaryVertexIndex] : nullptr;
          SvtxVertex *primaryVertexSource =
              hasPrimaryVertex ? primaryVertexSources[primaryVertexIndex] : nullptr;

          DFinderCandidateRecord candidate;
          candidate.run = run;
          candidate.event = event;
          candidate.centrality = centrality;
          candidate.nTracks = trackMap->size();
          candidate.nPV = primaryVertices.size();
          candidate.daughterCrossing = daughterCrossing;
          candidate.nVertexSameCrossing = matchingVertexIndices.size();
          candidate.pvMatched = hasPrimaryVertex &&
                                primaryVertexCrossings[primaryVertexIndex] == daughterCrossing;
          candidate.pvMatchedBy =
              m_config.vertexAssociation == DFinderConfig::VertexAssociation::MatchedCrossing
                  ? selectedPvMatchedBy
                  : pvMatchRuleUnavailable;
          if (hasPrimaryVertex)
          {
            candidate.PVx = primaryVertex->GetX();
            candidate.PVy = primaryVertex->GetY();
            candidate.PVz = primaryVertex->GetZ();
            candidate.PVchi2 = primaryVertex->GetChi2();
            candidate.PVndf = primaryVertex->GetNDF();
            candidate.PVid = primaryVertex->Id();
            candidate.pvCrossing = primaryVertexCrossings[primaryVertexIndex];
            candidate.pvNtracks = primaryVertexSource ? primaryVertexSource->size_tracks() : 0;
          }
          candidate.unfittedMass = unfittedMass;
          candidate.unfittedPt = unfittedPt;
          candidate.maxDoca = maximumDoca;
          candidate.maxDocaXY = maximumDocaXY;

          KFParticle &mother = fitResult.mother;
          candidate.mass = mother.GetMass();
          candidate.pt = mother.GetPt();
          candidate.eta = mother.GetEta();
          candidate.phi = mother.GetPhi();
          candidate.rapidity = mother.GetRapidity();
          candidate.px = mother.GetPx();
          candidate.py = mother.GetPy();
          candidate.pz = mother.GetPz();
          candidate.vtxX = mother.GetX();
          candidate.vtxY = mother.GetY();
          candidate.vtxZ = mother.GetZ();
          candidate.vtxXErr = mother.GetCovariance(0, 0);
          candidate.vtxYErr = mother.GetCovariance(1, 1);
          candidate.vtxZErr = mother.GetCovariance(2, 2);
          candidate.vtxYXErr = mother.GetCovariance(1, 0);
          candidate.vtxZXErr = mother.GetCovariance(2, 0);
          candidate.vtxZYErr = mother.GetCovariance(2, 1);
          candidate.vtxchi2 = mother.GetChi2();
          candidate.vtxdof = mother.GetNDF();
          candidate.massErr = mother.GetErrMass();
          candidate.ellipsoidVolume = calculateEllipsoidVolume(mother);
          candidate.nIntermediate =
              static_cast<int>(fitResult.intermediateParticles.size());
          for (size_t index = 0; index < fitResult.intermediateParticles.size(); ++index)
          {
            const KFParticle &intermediate = fitResult.intermediateParticles[index];
            DFinderIntermediateRecord &intermediateRecord = candidate.intermediates[index];
            intermediateRecord.pdgId = intermediate.GetPDG();
            intermediateRecord.parentIntermediateIndex =
                fitResult.intermediateParentIndices[index];
            intermediateRecord.mass = intermediate.GetMass();
            intermediateRecord.massErr = intermediate.GetErrMass();
            intermediateRecord.pt = intermediate.GetPt();
            intermediateRecord.eta = intermediate.GetEta();
            intermediateRecord.phi = intermediate.GetPhi();
            intermediateRecord.px = intermediate.GetPx();
            intermediateRecord.py = intermediate.GetPy();
            intermediateRecord.pz = intermediate.GetPz();
            intermediateRecord.vtxX = intermediate.GetX();
            intermediateRecord.vtxY = intermediate.GetY();
            intermediateRecord.vtxZ = intermediate.GetZ();
            intermediateRecord.vtxchi2 = intermediate.GetChi2();
            intermediateRecord.vtxdof = intermediate.GetNDF();
          }

          if (hasPrimaryVertex)
          {
            const double displacement[3]{
                mother.GetX() - primaryVertex->GetX(),
                mother.GetY() - primaryVertex->GetY(),
                mother.GetZ() - primaryVertex->GetZ()};
            const double distanceSquared =
                displacement[0] * displacement[0] +
                displacement[1] * displacement[1] +
                displacement[2] * displacement[2];
            candidate.svpvDistance = std::sqrt(distanceSquared);
            double distanceVariance = 0.;
            for (int row = 0; row < 3; ++row)
            {
              for (int column = 0; column < 3; ++column)
              {
                distanceVariance +=
                    displacement[row] * displacement[column] *
                    (mother.GetCovariance(row, column) +
                     primaryVertex->GetCovariance(row, column));
              }
            }
            if (distanceSquared > 0. && distanceVariance >= 0.)
            {
              candidate.svpvDisErr =
                  std::sqrt(distanceVariance / distanceSquared);
            }

            const double distanceSquaredXY =
                displacement[0] * displacement[0] +
                displacement[1] * displacement[1];
            candidate.svpvDistance2D = std::sqrt(distanceSquaredXY);
            double distanceVarianceXY = 0.;
            for (int row = 0; row < 2; ++row)
            {
              for (int column = 0; column < 2; ++column)
              {
                distanceVarianceXY +=
                    displacement[row] * displacement[column] *
                    (mother.GetCovariance(row, column) +
                     primaryVertex->GetCovariance(row, column));
              }
            }
            if (distanceSquaredXY > 0. && distanceVarianceXY >= 0.)
            {
              candidate.svpvDisErr2D =
                  std::sqrt(distanceVarianceXY / distanceSquaredXY);
            }

            candidate.dira = eventDIRA(mother, *primaryVertex);
            candidate.diraXY = eventDIRA(mother, *primaryVertex, false);
            candidate.alpha = std::acos(std::clamp(candidate.dira, -1.F, 1.F));
            candidate.fdchi2 = flightDistanceChi2(mother, *primaryVertex);
            candidate.PVdca = mother.GetDistanceFromVertex(*primaryVertex);
            candidate.PVdcaStddev = mother.GetDeviationFromVertex(*primaryVertex);
            candidate.PVdcaXY = std::abs(mother.GetDistanceFromVertexXY(*primaryVertex));
            candidate.PVdcaStddevXY = mother.GetDeviationFromVertexXY(*primaryVertex);
            candidate.pseudoProperDecayTime =
                mother.GetPseudoProperDecayTime(*primaryVertex, candidate.mass);

            KFParticle constrainedMother = mother;
            constrainedMother.SetProductionVertex(*primaryVertex);
            constrainedMother.TransportToDecayVertex();
            constrainedMother.GetDecayLength(
                candidate.decayLength, candidate.decayLengthErr);
            constrainedMother.GetDecayLengthXY(
                candidate.decayLengthXY, candidate.decayLengthErrXY);
            constrainedMother.GetLifeTime(
                candidate.decayTime, candidate.decayTimeErr);
            constexpr float speedOfLightCmPerPs = 2.99792458e-2F;
            candidate.decayTime /= speedOfLightCmPerPs;
            candidate.decayTimeErr /= speedOfLightCmPerPs;
          }
          else
          {
            const float nan = std::numeric_limits<float>::quiet_NaN();
            candidate.PVx = nan;
            candidate.PVy = nan;
            candidate.PVz = nan;
            candidate.PVchi2 = nan;
            candidate.svpvDistance = nan;
            candidate.svpvDisErr = nan;
            candidate.svpvDistance2D = nan;
            candidate.svpvDisErr2D = nan;
            candidate.alpha = nan;
            candidate.fdchi2 = nan;
            candidate.dira = nan;
            candidate.diraXY = nan;
            candidate.decayLength = nan;
            candidate.decayLengthErr = nan;
            candidate.decayLengthXY = nan;
            candidate.decayLengthErrXY = nan;
            candidate.decayTime = nan;
            candidate.decayTimeErr = nan;
            candidate.pseudoProperDecayTime = nan;
            candidate.PVdca = nan;
            candidate.PVdcaStddev = nan;
            candidate.PVdcaXY = nan;
            candidate.PVdcaStddevXY = nan;
          }

          for (size_t index = 0; index < nProngs; ++index)
          {
            const TrackRecord &trackRecord = selectedTracks[trackIndices[index]];
            const KFParticle &daughter = fitResult.refittedDaughters[index];
            DFinderDaughterRecord &daughterRecord = candidate.daughters[index];
            daughterRecord.pt = daughter.GetPt();
            daughterRecord.eta = daughter.GetEta();
            daughterRecord.phi = daughter.GetPhi();
            daughterRecord.mass = daughter.GetMass();
            daughterRecord.index = trackRecord.track->get_id();
            daughterRecord.parentIntermediateIndex =
                hypothesis.prongs[index].parentIntermediateIndex;
            daughterRecord.massHypothesis = hypothesis.prongs[index].mass;
            daughterRecord.charge = trackRecord.track->get_charge();
            daughterRecord.chi2ndof =
                trackRecord.track->get_chisq() /
                static_cast<float>(trackRecord.track->get_ndf());
            daughterRecord.nMvtx = trackRecord.nMvtx;
            daughterRecord.nIntt = trackRecord.nIntt;
            daughterRecord.nTpc = trackRecord.nTpc;
            daughterRecord.crossing = trackRecord.track->get_crossing();
            daughterRecord.dEdx = trackRecord.dEdx;
            if (hasPrimaryVertex)
            {
              daughterRecord.ip3d = daughter.GetDistanceFromVertex(*primaryVertex);
              daughterRecord.ip3dSig = daughter.GetDeviationFromVertex(*primaryVertex);
              daughterRecord.ipxy =
                  std::abs(daughter.GetDistanceFromVertexXY(*primaryVertex));
              daughterRecord.ipxySig =
                  daughter.GetDeviationFromVertexXY(*primaryVertex);
            }
            else
            {
              const float nan = std::numeric_limits<float>::quiet_NaN();
              daughterRecord.ip3d = nan;
              daughterRecord.ip3dSig = nan;
              daughterRecord.ipxy = nan;
              daughterRecord.ipxySig = nan;
            }
          }

          candidate.isTruthMatched = isTruthMatched;
          if (truthMother)
          {
            candidate.truthPdgId = truthMother->get_pid();
            candidate.truthPt =
                std::hypot(truthMother->get_px(), truthMother->get_py());
            if (truthMother->get_e() > std::abs(truthMother->get_pz()))
            {
              candidate.truthY =
                  0.5F * std::log(
                             (truthMother->get_e() + truthMother->get_pz()) /
                             (truthMother->get_e() - truthMother->get_pz()));
            }
            PHG4VtxPoint *productionVertex =
                truthInfo->GetVtx(truthMother->get_vtx_id());
            PHG4Particle *truthDecayProduct =
                matchTruthDecay(
                    m_decay.mother.daughters.front(),
                    hypothesis.chargeConjugated,
                    truthByProng,
                    truthInfo);
            PHG4VtxPoint *decayVertex =
                truthDecayProduct
                    ? truthInfo->GetVtx(truthDecayProduct->get_vtx_id())
                    : nullptr;
            if (productionVertex && decayVertex)
            {
              candidate.truthDecayLength =
                  std::sqrt(
                      std::pow(decayVertex->get_x() - productionVertex->get_x(), 2) +
                      std::pow(decayVertex->get_y() - productionVertex->get_y(), 2) +
                      std::pow(decayVertex->get_z() - productionVertex->get_z(), 2));
            }
            candidate.isPrompt = truthInfo->is_primary(truthMother);
          }

          outputRecords.push_back(candidate);
        }
        if (truncated)
        {
          break;
        }
      }

      if (!truncated)
      {
        haveCombination = false;
        for (size_t reverseIndex = nProngs; reverseIndex > 0; --reverseIndex)
        {
          const size_t index = reverseIndex - 1;
          const size_t maximumPosition = trackGroup.size() - nProngs + index;
          if (combinationPositions[index] < maximumPosition)
          {
            ++combinationPositions[index];
            for (size_t following = index + 1; following < nProngs; ++following)
            {
              combinationPositions[following] = combinationPositions[following - 1] + 1;
            }
            haveCombination = true;
            break;
          }
        }
      }
    }
  }

  if (truncated)
  {
    ++m_cutFlow.truncatedEvents;
    std::cerr << Name() << ": candidate output truncated at "
              << m_config.maxCandidatesPerEvent << " rows in run "
              << run << ", event " << event << std::endl;
  }
  for (DFinderCandidateRecord &candidate : outputRecords)
  {
    candidate.truncated = truncated;
    m_ntuple.fillCandidate(candidate);
    ++m_cutFlow.rowsWritten;
  }

  DFinderEventRecord eventRecord;
  eventRecord.run = run;
  eventRecord.event = event;
  eventRecord.nTracks = trackMap->size();
  eventRecord.nSelectedTracks = selectedTracks.size();
  eventRecord.nPV = primaryVertices.size();
  eventRecord.combinations = m_cutFlow.combinations - startCombinations;
  eventRecord.fitsAttempted = m_cutFlow.fitsAttempted - startFitsAttempted;
  eventRecord.fitsValid = m_cutFlow.fitsValid - startFitsValid;
  eventRecord.rowsWritten = m_cutFlow.rowsWritten - startRowsWritten;
  eventRecord.truncated = truncated;
  m_ntuple.fillEvent(eventRecord);

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": run " << run << ", event " << event
              << ", tracks " << trackMap->size()
              << ", selected " << selectedTracks.size()
              << ", PVs " << primaryVertices.size()
              << ", rows " << eventRecord.rowsWritten << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int DFinderSphenix::End(PHCompositeNode * /*topNode*/)
{
  m_ntuple.close(m_cutFlow);
  const unsigned long long acceptedPairs = m_cutFlow.osPairs + m_cutFlow.ssPairs;
  const double sameSignFraction =
      acceptedPairs > 0
          ? static_cast<double>(m_cutFlow.ssPairs) / static_cast<double>(acceptedPairs)
          : 0.;
  std::cout << Name() << " cumulative cut flow\n"
            << "  events: " << m_cutFlow.events << '\n'
            << "  tracks input: " << m_cutFlow.tracksInput << '\n'
            << "  tracks selected: " << m_cutFlow.tracksSelected << '\n'
            << "  tracks with undefined crossing (tracks; dropped in grouped mode): " << m_cutFlow.tracksCrossingUndefined << '\n'
            << "  tracks with crossing zero (tracks; kept): " << m_cutFlow.tracksCrossingZero << '\n'
            << "  crossing groups (groups with at least 2 selected tracks): " << m_cutFlow.crossingGroups << '\n'
            << "  unordered track combinations considered: " << m_cutFlow.combinations << '\n'
            << "  opposite-sign track pairs: " << m_cutFlow.osPairs << '\n'
            << "  same-sign track pairs: " << m_cutFlow.ssPairs << '\n'
            << "  pair x mass-hypothesis entries formed: " << m_cutFlow.hypothesesFormed << '\n'
            << "  pair x mass-hypothesis entries passing unfitted mass: " << m_cutFlow.massPass << '\n'
            << "  pair x mass-hypothesis entries passing unfitted pT and rapidity: " << m_cutFlow.kinematicPass << '\n'
            << "  pair x mass-hypothesis entries passing maximum pair DOCA: " << m_cutFlow.docaPass << '\n'
            << "  pair x mass-hypothesis propagation failures: " << m_cutFlow.propagationFailed << '\n'
            << "  pair x mass-hypothesis fits attempted: " << m_cutFlow.fitsAttempted << '\n'
            << "  pair x mass-hypothesis valid fits: " << m_cutFlow.fitsValid << '\n'
            << "  candidates without a matched production vertex (candidates): " << m_cutFlow.noMatchedVertex << '\n'
            << "  candidates with ambiguous same-crossing vertices (candidates): " << m_cutFlow.ambiguousVertex << '\n'
            << "  candidate rows written: " << m_cutFlow.rowsWritten << '\n'
            << "  truncated events: " << m_cutFlow.truncatedEvents << '\n'
            << "  same-sign fraction of accepted track pairs: "
            << sameSignFraction << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

bool DFinderSphenix::propagateToSecondaryVertex(
    const std::vector<KFParticle> &daughters,
    const std::vector<SvtxTrack *> &tracks,
    size_t firstDocaTrack,
    size_t secondDocaTrack,
    std::vector<KFParticle> &propagated) const
{
  propagated.clear();
  if (!m_actsGeometry || !m_svtxVertexMap ||
      daughters.size() != tracks.size() ||
      firstDocaTrack >= daughters.size() ||
      secondDocaTrack >= daughters.size())
  {
    return false;
  }

  float pathLengths[2]{0.F, 0.F};
  float derivatives[4][6]{};
  daughters[firstDocaTrack].GetDStoParticle(
      daughters[secondDocaTrack], pathLengths, derivatives);
  KFParticle firstAtDoca = daughters[firstDocaTrack];
  KFParticle secondAtDoca = daughters[secondDocaTrack];
  firstAtDoca.TransportToDS(pathLengths[0], derivatives[0]);
  secondAtDoca.TransportToDS(pathLengths[1], derivatives[3]);

  const Acts::Vector3 seed(
      0.5 * (firstAtDoca.GetX() + secondAtDoca.GetX()) * Acts::UnitConstants::cm,
      0.5 * (firstAtDoca.GetY() + secondAtDoca.GetY()) * Acts::UnitConstants::cm,
      0.5 * (firstAtDoca.GetZ() + secondAtDoca.GetZ()) * Acts::UnitConstants::cm);
  ActsPropagator propagator(m_actsGeometry);
  propagator.verbosity(Verbosity());
  const auto targetSurface = propagator.makeVertexSurface(seed);
  ActsTransformations transformations;

  propagated.reserve(tracks.size());
  for (size_t index = 0; index < tracks.size(); ++index)
  {
    auto inputParameters = propagator.makeTrackParams(tracks[index], m_svtxVertexMap);
    if (!inputParameters.ok())
    {
      propagated.clear();
      return false;
    }
    auto propagationResult =
        propagator.propagateTrack(inputParameters.value(), targetSurface);
    if (!propagationResult.ok())
    {
      propagated.clear();
      return false;
    }

    const auto &parameters = propagationResult.value().second;
    const auto position =
        parameters.position(m_actsGeometry->geometry().getGeoContext());
    const auto momentum = parameters.momentum();
    const auto globalCovariance =
        transformations.rotateActsCovToSvtxTrack(parameters);
    float state[6]{
        static_cast<float>(position.x() / Acts::UnitConstants::cm),
        static_cast<float>(position.y() / Acts::UnitConstants::cm),
        static_cast<float>(position.z() / Acts::UnitConstants::cm),
        static_cast<float>(momentum.x()),
        static_cast<float>(momentum.y()),
        static_cast<float>(momentum.z())};
    float covariance[21]{};
    unsigned int covarianceIndex = 0;
    for (unsigned int row = 0; row < 6; ++row)
    {
      for (unsigned int column = 0; column <= row; ++column)
      {
        covariance[covarianceIndex] = globalCovariance(row, column);
        ++covarianceIndex;
      }
    }

    KFParticle daughter;
    daughter.Create(state, covariance, tracks[index]->get_charge(), -1);
    daughter.NDF() = tracks[index]->get_ndf();
    daughter.Chi2() = tracks[index]->get_chisq();
    daughter.SetId(tracks[index]->get_id());
    propagated.push_back(daughter);
  }
  return true;
}

PHG4Particle *DFinderSphenix::matchTruthDecay(
    const DFinderDecayNode &node,
    const bool chargeConjugated,
    const std::vector<PHG4Particle *> &truthByProng,
    PHG4TruthInfoContainer *truthInfo) const
{
  const int expectedPdgId =
      chargeConjugated ? node.conjugatePdgId : node.pdgId;
  if (node.isFinalState())
  {
    if (node.prongIndex < 0 ||
        static_cast<size_t>(node.prongIndex) >= truthByProng.size())
    {
      return nullptr;
    }
    PHG4Particle *truth = truthByProng[node.prongIndex];
    return truth && truth->get_pid() == expectedPdgId ? truth : nullptr;
  }

  int truthParentId = 0;
  for (size_t index = 0; index < node.daughters.size(); ++index)
  {
    PHG4Particle *truthDaughter =
        matchTruthDecay(
            node.daughters[index],
            chargeConjugated,
            truthByProng,
            truthInfo);
    if (!truthDaughter)
    {
      return nullptr;
    }
    if (index == 0)
    {
      truthParentId = truthDaughter->get_parent_id();
    }
    else if (truthDaughter->get_parent_id() != truthParentId)
    {
      return nullptr;
    }
  }

  PHG4Particle *truthMother = truthInfo->GetParticle(truthParentId);
  return truthMother && truthMother->get_pid() == expectedPdgId
             ? truthMother
             : nullptr;
}
