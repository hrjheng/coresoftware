#ifndef DFINDERSPHENIX_DFINDERCONFIG_H
#define DFINDERSPHENIX_DFINDERCONFIG_H

#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

struct DFinderConfig
{
  // Vertex association strategy.
  //   AllVertices     - legacy: one row per candidate per vertex in the frame
  //   MatchedCrossing - one row per candidate, against the vertex whose beam
  //                     crossing equals the daughters' crossing
  enum class VertexAssociation
  {
    AllVertices,
    MatchedCrossing
  };

  static constexpr std::size_t minimumDaughters{2};
  static constexpr std::size_t maximumDaughters{5};
  static constexpr std::size_t maximumIntermediateStates{maximumDaughters - 2};

  std::string decayDescriptor{"[D0 -> K^- pi^+]cc"};
  std::string trackMapNodeName{"SvtxTrackMap"};
  std::string vertexMapNodeName{"SvtxVertexMap"};
  bool useMbdVertex{false};
  bool useGlobalVertex{true};
  VertexAssociation vertexAssociation{VertexAssociation::MatchedCrossing};

  // Require all daughters to carry the same, defined crossing before pairing.
  bool requireSameDaughterCrossing{true};

  // Keep unmatched candidates with vertex quantities set to NaN.
  bool keepUnmatchedCandidates{false};

  float trackMinPt{0.2F};
  float trackMaxEta{1.1F};
  float trackMaxChi2ndof{100.F};
  int minMvtxStates{0};
  int minInttStates{0};
  int minTpcStates{0};
  bool requireCrossingZero{false};

  // Also write same-sign track pairs, for use as an offline combinatorial
  // background reference. Same-sign pairs are processed identically to
  // opposite-sign pairs and are distinguished offline via rftk*_charge.
  bool keepSameSign{false};

  // Offsets from the descriptor mother mass, in GeV.
  float massWindowLo{-0.2F};
  float massWindowHi{0.2F};
  float motherMinPt{0.F};
  float motherMaxRapidity{std::numeric_limits<float>::max()};
  float maxPairDoca{std::numeric_limits<float>::max()};
  float maxPairDocaXY{std::numeric_limits<float>::max()};
  bool constrainIntermediateMasses{false};
  // Absolute GeV windows in descriptor preorder, one per intermediate.
  // An empty vector applies no intermediate-mass selection.
  std::vector<std::pair<float, float>> intermediateMassWindows;

  unsigned int maxCandidatesPerEvent{100000};
  std::string outputFileName{"dfinder.root"};
  bool doPropagateToSV{false};
  std::string magneticFieldFile{"FIELDMAP_TRACKING"};
};

struct DFinderCutFlow
{
  unsigned long long events{0};
  unsigned long long tracksInput{0};
  unsigned long long tracksSelected{0};
  unsigned long long tracksCrossingUndefined{0};
  unsigned long long tracksCrossingZero{0};
  unsigned long long crossingGroups{0};
  unsigned long long combinations{0};
  unsigned long long osPairs{0};
  unsigned long long ssPairs{0};
  unsigned long long hypothesesFormed{0};
  unsigned long long massPass{0};
  unsigned long long kinematicPass{0};
  unsigned long long docaPass{0};
  unsigned long long propagationFailed{0};
  unsigned long long fitsAttempted{0};
  unsigned long long fitsValid{0};
  unsigned long long noMatchedVertex{0};
  unsigned long long ambiguousVertex{0};
  unsigned long long rowsWritten{0};
  unsigned long long truncatedEvents{0};
};

#endif
