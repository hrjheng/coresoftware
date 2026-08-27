#include "DFinderCombinatorics.h"
#include "DFinderConfig.h"

#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <set>
#include <sstream>
#include <utility>

bool DFinderCombinatorics::parseDecayDescriptor(
    const std::string &descriptor,
    DFinderDecay &decay,
    std::string &error) const
{
  decay = DFinderDecay{};
  error.clear();

  std::string compact;
  compact.reserve(descriptor.size());
  for (const char character : descriptor)
  {
    if (!std::isspace(static_cast<unsigned char>(character)))
    {
      compact.push_back(character);
    }
  }

  if (compact.empty())
  {
    error = "descriptor is empty";
    return false;
  }

  if (compact.front() == '[')
  {
    const std::size_t closingBracket = compact.rfind(']');
    if (closingBracket == std::string::npos)
    {
      error = "descriptor has an opening '[' without a closing ']'";
      return false;
    }

    std::string suffix = compact.substr(closingBracket + 1);
    std::transform(suffix.begin(), suffix.end(), suffix.begin(),
                   [](const unsigned char character)
                   { return static_cast<char>(std::toupper(character)); });
    if (!suffix.empty() && suffix != "CC")
    {
      error = "only the optional suffix 'cc' may follow a bracketed descriptor";
      return false;
    }
    decay.includeChargeConjugate = suffix == "CC";
    compact = compact.substr(1, closingBracket - 1);
  }

  std::size_t position = 0;
  if (!parseDecayNode(compact,
                      position,
                      true,
                      -1,
                      decay.mother,
                      decay.prongs,
                      decay.nIntermediateStates,
                      error))
  {
    return false;
  }
  if (position != compact.size())
  {
    error = "unexpected trailing text in decay descriptor";
    return false;
  }

  if (decay.prongs.size() < DFinderConfig::minimumDaughters ||
      decay.prongs.size() > DFinderConfig::maximumDaughters)
  {
    error = "descriptor must contain between two and five charged final-state prongs";
    return false;
  }
  if (decay.nIntermediateStates > DFinderConfig::maximumIntermediateStates)
  {
    error = "descriptor has too many intermediate states for the configured daughter limit";
    return false;
  }

  decay.motherName = decay.mother.name;
  decay.motherPdgId = decay.mother.pdgId;
  decay.motherMass = decay.mother.mass;
  return true;
}

bool DFinderCombinatorics::parseDecayNode(
    const std::string &descriptor,
    std::size_t &position,
    const bool isMother,
    const int parentIntermediateIndex,
    DFinderDecayNode &node,
    std::vector<DFinderProng> &prongs,
    std::size_t &nIntermediateStates,
    std::string &error) const
{
  const std::string decayArrow{"->"};
  const std::size_t arrowPosition = descriptor.find(decayArrow, position);
  if (arrowPosition == std::string::npos || arrowPosition == position)
  {
    error = "each mother or intermediate state must contain a particle name followed by '->'";
    return false;
  }

  node.name = descriptor.substr(position, arrowPosition - position);
  TParticlePDG *particle = TDatabasePDG::Instance()->GetParticle(node.name.c_str());
  if (!particle)
  {
    error = "unknown mother or intermediate particle '" + node.name + "'";
    return false;
  }
  node.pdgId = particle->PdgCode();
  node.conjugatePdgId = particle->AntiParticle()
                            ? particle->AntiParticle()->PdgCode()
                            : particle->PdgCode();
  node.mass = particle->Mass();
  constexpr double rootChargeUnitsPerElementaryCharge = 3.;
  node.charge = static_cast<int>(
      std::lround(particle->Charge() / rootChargeUnitsPerElementaryCharge));
  node.parentIntermediateIndex = parentIntermediateIndex;
  if (!isMother)
  {
    node.intermediateIndex = static_cast<int>(nIntermediateStates);
    ++nIntermediateStates;
  }

  position = arrowPosition + decayArrow.size();
  const std::set<std::string> signedRootNames{"e", "mu", "pi", "K"};
  while (position < descriptor.size())
  {
    if (descriptor[position] == '}')
    {
      if (isMother)
      {
        error = "unexpected closing '}' in mother decay";
        return false;
      }
      ++position;
      break;
    }

    if (descriptor[position] == '{')
    {
      ++position;
      DFinderDecayNode daughterNode;
      if (!parseDecayNode(descriptor,
                          position,
                          false,
                          node.intermediateIndex,
                          daughterNode,
                          prongs,
                          nIntermediateStates,
                          error))
      {
        return false;
      }
      node.daughters.push_back(std::move(daughterNode));
      continue;
    }

    const std::size_t chargePosition = descriptor.find('^', position);
    if (chargePosition == std::string::npos || chargePosition == position ||
        chargePosition + 1 >= descriptor.size())
    {
      error = "each final-state prong must have a ^ charge";
      return false;
    }
    const std::size_t nextBrace = descriptor.find_first_of("{}", position);
    if (nextBrace != std::string::npos && nextBrace < chargePosition)
    {
      error = "malformed intermediate-state braces";
      return false;
    }

    const std::string baseName = descriptor.substr(position, chargePosition - position);
    const char chargeCharacter = descriptor[chargePosition + 1];
    int charge = 0;
    if (chargeCharacter == '+')
    {
      charge = 1;
    }
    else if (chargeCharacter == '-')
    {
      charge = -1;
    }
    else
    {
      error = "only charged reconstructed prongs (^+ or ^-) are supported";
      return false;
    }

    std::string rootName = baseName;
    if (signedRootNames.find(baseName) != signedRootNames.end())
    {
      rootName.push_back(chargeCharacter);
    }

    TParticlePDG *daughterParticle = TDatabasePDG::Instance()->GetParticle(rootName.c_str());
    if (!daughterParticle)
    {
      error = "unknown daughter particle '" + rootName + "'";
      return false;
    }
    const int databaseCharge = static_cast<int>(
        std::lround(daughterParticle->Charge() / rootChargeUnitsPerElementaryCharge));
    if (databaseCharge != charge)
    {
      std::ostringstream message;
      message << "charge " << chargeCharacter
              << " is inconsistent with particle '" << rootName << "'";
      error = message.str();
      return false;
    }

    DFinderDecayNode daughterNode;
    daughterNode.name = rootName;
    daughterNode.pdgId = daughterParticle->PdgCode();
    daughterNode.conjugatePdgId = daughterParticle->AntiParticle()
                                      ? daughterParticle->AntiParticle()->PdgCode()
                                      : daughterParticle->PdgCode();
    daughterNode.mass = daughterParticle->Mass();
    daughterNode.charge = charge;
    daughterNode.prongIndex = static_cast<int>(prongs.size());
    daughterNode.parentIntermediateIndex = node.intermediateIndex;
    node.daughters.push_back(daughterNode);
    prongs.push_back({daughterNode.name,
                      daughterNode.pdgId,
                      daughterNode.conjugatePdgId,
                      daughterNode.mass,
                      daughterNode.charge,
                      daughterNode.prongIndex,
                      daughterNode.parentIntermediateIndex});
    position = chargePosition + 2;
  }

  if (!isMother && (position == 0 || descriptor[position - 1] != '}'))
  {
    error = "intermediate state has an opening '{' without a closing '}'";
    return false;
  }
  if (node.daughters.size() < DFinderConfig::minimumDaughters)
  {
    error = "each mother or intermediate state must have at least two daughters";
    return false;
  }

  int daughterCharge = 0;
  for (const DFinderDecayNode &daughter : node.daughters)
  {
    daughterCharge += daughter.charge;
  }
  if (daughterCharge != node.charge)
  {
    std::ostringstream message;
    message << "daughter charge " << daughterCharge
            << " is inconsistent with particle '" << node.name
            << "' charge " << node.charge;
    error = message.str();
    return false;
  }
  return true;
}

std::vector<DFinderHypothesis> DFinderCombinatorics::makePidPermutations(
    const DFinderDecay &decay) const
{
  std::vector<DFinderHypothesis> hypotheses;
  std::set<std::vector<std::pair<int, int>>> seenAssignments;

  std::vector<DFinderHypothesis> chargeModes{{decay.prongs, false}};
  if (decay.includeChargeConjugate)
  {
    std::vector<DFinderProng> conjugate = decay.prongs;
    for (DFinderProng &prong : conjugate)
    {
      std::swap(prong.pdgId, prong.conjugatePdgId);
      prong.charge *= -1;
      if (TParticlePDG *particle = TDatabasePDG::Instance()->GetParticle(prong.pdgId))
      {
        prong.name = particle->GetName();
      }
    }
    chargeModes.push_back({std::move(conjugate), true});
  }

  for (const DFinderHypothesis &chargeMode : chargeModes)
  {
    std::vector<std::size_t> order(chargeMode.prongs.size());
    for (std::size_t index = 0; index < order.size(); ++index)
    {
      order[index] = index;
    }

    do
    {
      DFinderHypothesis hypothesis;
      hypothesis.chargeConjugated = chargeMode.chargeConjugated;
      std::vector<std::pair<int, int>> assignment;
      hypothesis.prongs.reserve(order.size());
      assignment.reserve(order.size());
      for (const std::size_t index : order)
      {
        const DFinderProng &prong = chargeMode.prongs[index];
        hypothesis.prongs.push_back(prong);
        assignment.emplace_back(prong.pdgId, prong.parentIntermediateIndex);
      }

      if (seenAssignments.insert(assignment).second)
      {
        hypotheses.push_back(std::move(hypothesis));
      }
    } while (std::next_permutation(order.begin(), order.end()));
  }

  return hypotheses;
}
