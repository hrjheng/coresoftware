#include "DFinderFitter.h"
#include "DFinderConfig.h"

#include <cmath>
#include <vector>

DFinderFitResult DFinderFitter::fit(
    const std::vector<KFParticle> &daughters,
    const DFinderDecay &decay,
    const DFinderHypothesis &hypothesis,
    const bool constrainIntermediateMasses,
    const std::vector<std::pair<float, float>> &intermediateMassWindows) const
{
  DFinderFitResult result;
  if (daughters.size() < DFinderConfig::minimumDaughters ||
      daughters.size() != hypothesis.prongs.size() ||
      daughters.size() != decay.prongs.size() ||
      (!intermediateMassWindows.empty() &&
       intermediateMassWindows.size() != decay.nIntermediateStates))
  {
    return result;
  }

  std::vector<KFParticle> leafParticles(daughters.size());
  std::vector<bool> assignedLeaves(daughters.size(), false);
  for (std::size_t trackIndex = 0; trackIndex < daughters.size(); ++trackIndex)
  {
    const DFinderProng &prong = hypothesis.prongs[trackIndex];
    if (prong.prongIndex < 0 ||
        static_cast<std::size_t>(prong.prongIndex) >= leafParticles.size() ||
        assignedLeaves[prong.prongIndex])
    {
      return result;
    }

    KFParticle inputDaughter = daughters[trackIndex];
    KFParticle daughter;
    daughter.Create(inputDaughter.Parameters(),
                    inputDaughter.CovarianceMatrix(),
                    static_cast<int>(inputDaughter.GetQ()),
                    prong.mass);
    daughter.NDF() = inputDaughter.GetNDF();
    daughter.Chi2() = inputDaughter.GetChi2();
    daughter.SetId(inputDaughter.Id());
    daughter.SetPDG(prong.pdgId);
    leafParticles[prong.prongIndex] = daughter;
    assignedLeaves[prong.prongIndex] = true;
  }

  std::vector<KFParticle> refittedByProng(daughters.size());
  result.intermediateParticles.resize(decay.nIntermediateStates);
  result.intermediateParentIndices.resize(decay.nIntermediateStates, -1);
  if (!fitNode(decay.mother,
               leafParticles,
               hypothesis.chargeConjugated,
               constrainIntermediateMasses,
               intermediateMassWindows,
               refittedByProng,
               result.intermediateParticles,
               result.intermediateParentIndices,
               result.mother))
  {
    return DFinderFitResult{};
  }

  result.refittedDaughters.resize(daughters.size());
  for (std::size_t trackIndex = 0; trackIndex < daughters.size(); ++trackIndex)
  {
    result.refittedDaughters[trackIndex] =
        refittedByProng[hypothesis.prongs[trackIndex].prongIndex];
  }
  result.valid = true;
  return result;
}

bool DFinderFitter::fitNode(
    const DFinderDecayNode &node,
    const std::vector<KFParticle> &leafParticles,
    const bool chargeConjugated,
    const bool constrainIntermediateMasses,
    const std::vector<std::pair<float, float>> &intermediateMassWindows,
    std::vector<KFParticle> &refittedByProng,
    std::vector<KFParticle> &intermediateParticles,
    std::vector<int> &intermediateParentIndices,
    KFParticle &fittedParticle) const
{
  if (node.isFinalState())
  {
    if (node.prongIndex < 0 ||
        static_cast<std::size_t>(node.prongIndex) >= leafParticles.size())
    {
      return false;
    }
    fittedParticle = leafParticles[node.prongIndex];
    return true;
  }

  std::vector<KFParticle> fittedDaughters(node.daughters.size());
  for (std::size_t index = 0; index < node.daughters.size(); ++index)
  {
    if (!fitNode(node.daughters[index],
                 leafParticles,
                 chargeConjugated,
                 constrainIntermediateMasses,
                 intermediateMassWindows,
                 refittedByProng,
                 intermediateParticles,
                 intermediateParentIndices,
                 fittedDaughters[index]))
    {
      return false;
    }
  }

  constexpr int kfConstructMethod = 2;
  KFParticle composite;
  composite.SetConstructMethod(kfConstructMethod);
  for (std::size_t index = 0; index < fittedDaughters.size(); ++index)
  {
    KFParticle daughterForFit = fittedDaughters[index];
    if (constrainIntermediateMasses && !node.daughters[index].isFinalState())
    {
      const DFinderDecayNode &daughterNode = node.daughters[index];
      KFParticle constrainedDaughter;
      constrainedDaughter.Create(daughterForFit.Parameters(),
                                 daughterForFit.CovarianceMatrix(),
                                 static_cast<int>(daughterForFit.GetQ()),
                                 daughterNode.mass);
      constrainedDaughter.NDF() = daughterForFit.GetNDF();
      constrainedDaughter.Chi2() = daughterForFit.GetChi2();
      constrainedDaughter.SetId(daughterForFit.Id());
      constrainedDaughter.SetPDG(
          chargeConjugated ? daughterNode.conjugatePdgId : daughterNode.pdgId);
      daughterForFit = constrainedDaughter;
    }
    composite.AddDaughter(daughterForFit);
    composite.AddDaughterId(daughterForFit.Id());
  }
  composite.SetPDG(chargeConjugated ? node.conjugatePdgId : node.pdgId);

  float mass = 0.F;
  float massError = 0.F;
  composite.GetMass(mass, massError);
  if (composite.GetNDF() <= 0 ||
      !std::isfinite(composite.GetChi2()) ||
      !std::isfinite(mass) ||
      !std::isfinite(massError))
  {
    return false;
  }

  if (node.intermediateIndex >= 0 && !intermediateMassWindows.empty())
  {
    const std::pair<float, float> &massWindow =
        intermediateMassWindows[node.intermediateIndex];
    if (mass < massWindow.first || mass > massWindow.second)
    {
      return false;
    }
  }

  for (std::size_t index = 0; index < fittedDaughters.size(); ++index)
  {
    KFParticle daughter = fittedDaughters[index];
    daughter.SetProductionVertex(composite);
    const DFinderDecayNode &daughterNode = node.daughters[index];
    if (daughterNode.isFinalState())
    {
      refittedByProng[daughterNode.prongIndex] = daughter;
    }
    else
    {
      intermediateParticles[daughterNode.intermediateIndex] = daughter;
      intermediateParentIndices[daughterNode.intermediateIndex] =
          node.intermediateIndex;
    }
  }

  fittedParticle = composite;
  return true;
}
