#ifndef DFINDERSPHENIX_DFINDERFITTER_H
#define DFINDERSPHENIX_DFINDERFITTER_H

#include "DFinderCombinatorics.h"

#include <KFParticle.h>

#include <utility>
#include <vector>

#ifndef HomogeneousField
#error "dfinder_sphenix and KFParticle must be compiled with HomogeneousField"
#endif

struct DFinderFitResult
{
  KFParticle mother;
  std::vector<KFParticle> refittedDaughters;
  std::vector<KFParticle> intermediateParticles;
  std::vector<int> intermediateParentIndices;
  bool valid{false};
};

class DFinderFitter
{
 public:
  DFinderFitResult fit(
      const std::vector<KFParticle> &daughters,
      const DFinderDecay &decay,
      const DFinderHypothesis &hypothesis,
      bool constrainIntermediateMasses,
      const std::vector<std::pair<float, float>> &intermediateMassWindows) const;

 private:
  bool fitNode(
      const DFinderDecayNode &node,
      const std::vector<KFParticle> &leafParticles,
      bool chargeConjugated,
      bool constrainIntermediateMasses,
      const std::vector<std::pair<float, float>> &intermediateMassWindows,
      std::vector<KFParticle> &refittedByProng,
      std::vector<KFParticle> &intermediateParticles,
      std::vector<int> &intermediateParentIndices,
      KFParticle &fittedParticle) const;
};

#endif
