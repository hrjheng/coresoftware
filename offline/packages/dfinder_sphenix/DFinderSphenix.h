#ifndef DFINDERSPHENIX_DFINDERSPHENIX_H
#define DFINDERSPHENIX_DFINDERSPHENIX_H

#include "DFinderCombinatorics.h"
#include "DFinderConfig.h"
#include "DFinderFitter.h"
#include "DFinderNtuple.h"

#include <kfparticle_sphenix/KFParticle_Tools.h>

#include <fun4all/SubsysReco.h>

#include <string>
#include <utility>
#include <vector>

class ActsGeometry;
class PHCompositeNode;
class PHG4Particle;
class PHG4TruthInfoContainer;
class SvtxTrack;
class SvtxVertexMap;

class DFinderSphenix : public SubsysReco, protected KFParticle_Tools
{
 public:
  explicit DFinderSphenix(const std::string &name = "DFinderSphenix",
                          const DFinderConfig &config = DFinderConfig{});
  ~DFinderSphenix() override = default;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void setDecayDescriptor(const std::string &value) { m_config.decayDescriptor = value; }
  void setTrackMapNodeName(const std::string &value) { m_config.trackMapNodeName = value; }
  void setVertexMapNodeName(const std::string &value) { m_config.vertexMapNodeName = value; }
  void useMbdVertex(bool value = true) { m_config.useMbdVertex = value; }
  void useGlobalVertex(bool value = true) { m_config.useGlobalVertex = value; }
  void setTrackMinPt(float value) { m_config.trackMinPt = value; }
  void setTrackMaxEta(float value) { m_config.trackMaxEta = value; }
  void setTrackMaxChi2ndof(float value) { m_config.trackMaxChi2ndof = value; }
  void setMinMvtxStates(int value) { m_config.minMvtxStates = value; }
  void setMinInttStates(int value) { m_config.minInttStates = value; }
  void setMinTpcStates(int value) { m_config.minTpcStates = value; }
  void requireCrossingZero(bool value = true) { m_config.requireCrossingZero = value; }
  void setMassWindow(float low, float high)
  {
    m_config.massWindowLo = low;
    m_config.massWindowHi = high;
  }
  void setMotherMinPt(float value) { m_config.motherMinPt = value; }
  void setMotherMaxRapidity(float value) { m_config.motherMaxRapidity = value; }
  void setMaxPairDoca(float value) { m_config.maxPairDoca = value; }
  void setMaxPairDocaXY(float value) { m_config.maxPairDocaXY = value; }
  void setMaxCandidatesPerEvent(unsigned int value) { m_config.maxCandidatesPerEvent = value; }
  void constrainIntermediateMasses(bool value = true) { m_config.constrainIntermediateMasses = value; }
  void setIntermediateMassWindows(const std::vector<std::pair<float, float>> &value)
  {
    m_config.intermediateMassWindows = value;
  }
  void setOutputFileName(const std::string &value) { m_config.outputFileName = value; }
  void doPropagateToSV(bool value = true) { m_config.doPropagateToSV = value; }
  void setMagneticFieldFile(const std::string &value) { m_config.magneticFieldFile = value; }

 private:
  struct TrackRecord
  {
    SvtxTrack *track{nullptr};
    KFParticle particle;
    int nMvtx{0};
    int nIntt{0};
    int nTpc{0};
    float dEdx{-1.F};
    PHG4Particle *truth{nullptr};
  };

  bool propagateToSecondaryVertex(const std::vector<KFParticle> &daughters,
                                  const std::vector<SvtxTrack *> &tracks,
                                  size_t firstDocaTrack,
                                  size_t secondDocaTrack,
                                  std::vector<KFParticle> &propagated) const;
  PHG4Particle *matchTruthDecay(const DFinderDecayNode &node,
                               bool chargeConjugated,
                               const std::vector<PHG4Particle *> &truthByProng,
                               PHG4TruthInfoContainer *truthInfo) const;

  DFinderConfig m_config;
  DFinderDecay m_decay;
  std::vector<DFinderHypothesis> m_pidPermutations;
  DFinderCombinatorics m_combinatorics;
  DFinderFitter m_fitter;
  DFinderNtuple m_ntuple;
  DFinderCutFlow m_cutFlow;
  ActsGeometry *m_actsGeometry{nullptr};
  SvtxVertexMap *m_svtxVertexMap{nullptr};
};

#endif
