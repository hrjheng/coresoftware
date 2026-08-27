#ifndef DFINDERSPHENIX_DFINDERCOMBINATORICS_H
#define DFINDERSPHENIX_DFINDERCOMBINATORICS_H

#include <cstddef>
#include <string>
#include <vector>

struct DFinderProng
{
  std::string name;
  int pdgId{0};
  int conjugatePdgId{0};
  float mass{0.F};
  int charge{0};
  int prongIndex{-1};
  int parentIntermediateIndex{-1};
};

struct DFinderDecayNode
{
  std::string name;
  int pdgId{0};
  int conjugatePdgId{0};
  float mass{0.F};
  int charge{0};
  int prongIndex{-1};
  int intermediateIndex{-1};
  int parentIntermediateIndex{-1};
  std::vector<DFinderDecayNode> daughters;

  bool isFinalState() const { return daughters.empty(); }
};

struct DFinderDecay
{
  std::string motherName;
  int motherPdgId{0};
  float motherMass{0.F};
  bool includeChargeConjugate{false};
  std::vector<DFinderProng> prongs;
  DFinderDecayNode mother;
  std::size_t nIntermediateStates{0};
};

struct DFinderHypothesis
{
  std::vector<DFinderProng> prongs;
  bool chargeConjugated{false};
};

class DFinderCombinatorics
{
 public:
  bool parseDecayDescriptor(const std::string &descriptor, DFinderDecay &decay, std::string &error) const;

  std::vector<DFinderHypothesis> makePidPermutations(const DFinderDecay &decay) const;

 private:
  bool parseDecayNode(const std::string &descriptor,
                      std::size_t &position,
                      bool isMother,
                      int parentIntermediateIndex,
                      DFinderDecayNode &node,
                      std::vector<DFinderProng> &prongs,
                      std::size_t &nIntermediateStates,
                      std::string &error) const;
};

#endif
