#ifndef DFINDERSPHENIX_DFINDERNTUPLE_H
#define DFINDERSPHENIX_DFINDERNTUPLE_H

#include "DFinderConfig.h"

#include <array>
#include <memory>
#include <string>

class TFile;
class TTree;

struct DFinderDaughterRecord
{
  float pt{0.F};
  float eta{0.F};
  float phi{0.F};
  float mass{0.F};
  int index{-1};
  int parentIntermediateIndex{-1};
  float massHypothesis{0.F};
  int charge{0};
  float chi2ndof{0.F};
  int nMvtx{0};
  int nIntt{0};
  int nTpc{0};
  int crossing{0};
  float dEdx{0.F};
  float ip3d{0.F};
  float ip3dSig{0.F};
  float ipxy{0.F};
  float ipxySig{0.F};
};

struct DFinderIntermediateRecord
{
  int pdgId{0};
  int parentIntermediateIndex{-1};
  float mass{0.F};
  float massErr{0.F};
  float pt{0.F};
  float eta{0.F};
  float phi{0.F};
  float px{0.F};
  float py{0.F};
  float pz{0.F};
  float vtxX{0.F};
  float vtxY{0.F};
  float vtxZ{0.F};
  float vtxchi2{0.F};
  int vtxdof{0};
};

struct DFinderCandidateRecord
{
  int run{-1};
  int event{-1};
  float centrality{-99.F};
  int nTracks{0};
  int nPV{0};
  float PVx{0.F};
  float PVy{0.F};
  float PVz{0.F};
  float PVchi2{0.F};
  int PVndf{0};
  int PVid{-1};
  int pvCrossing{std::numeric_limits<short int>::max()};
  int daughterCrossing{std::numeric_limits<short int>::max()};
  int pvMatched{0};
  int nVertexSameCrossing{0};
  int pvMatchedBy{-1};
  int pvNtracks{0};
  bool truncated{false};

  float mass{0.F};
  float unfittedMass{0.F};
  float pt{0.F};
  float unfittedPt{0.F};
  float eta{0.F};
  float phi{0.F};
  float rapidity{0.F};
  float px{0.F};
  float py{0.F};
  float pz{0.F};
  float vtxX{0.F};
  float vtxY{0.F};
  float vtxZ{0.F};
  float vtxXErr{0.F};
  float vtxYErr{0.F};
  float vtxZErr{0.F};
  float vtxYXErr{0.F};
  float vtxZXErr{0.F};
  float vtxZYErr{0.F};
  float vtxchi2{0.F};
  int vtxdof{0};
  float svpvDistance{0.F};
  float svpvDisErr{0.F};
  float svpvDistance2D{0.F};
  float svpvDisErr2D{0.F};
  float alpha{0.F};
  float maxDoca{0.F};
  float maxDocaXY{0.F};

  float fdchi2{0.F};
  float dira{0.F};
  float diraXY{0.F};
  float decayLength{0.F};
  float decayLengthErr{0.F};
  float decayLengthXY{0.F};
  float decayLengthErrXY{0.F};
  float decayTime{0.F};
  float decayTimeErr{0.F};
  float pseudoProperDecayTime{0.F};
  float PVdca{0.F};
  float PVdcaStddev{0.F};
  float PVdcaXY{0.F};
  float PVdcaStddevXY{0.F};
  float massErr{0.F};
  float ellipsoidVolume{0.F};

  int nIntermediate{0};
  std::array<DFinderDaughterRecord, DFinderConfig::maximumDaughters> daughters;
  std::array<DFinderIntermediateRecord, DFinderConfig::maximumIntermediateStates> intermediates;

  bool isTruthMatched{false};
  int truthPdgId{0};
  float truthPt{0.F};
  float truthY{0.F};
  float truthDecayLength{0.F};
  bool isPrompt{false};
};

struct DFinderEventRecord
{
  int run{-1};
  int event{-1};
  int nTracks{0};
  int nSelectedTracks{0};
  int nPV{0};
  unsigned long long combinations{0};
  unsigned long long fitsAttempted{0};
  unsigned long long fitsValid{0};
  unsigned long long rowsWritten{0};
  bool truncated{false};
};

class DFinderNtuple
{
 public:
  bool open(const std::string &fileName, bool keepSameSign);
  void fillCandidate(const DFinderCandidateRecord &candidate);
  void fillEvent(const DFinderEventRecord &event);
  void close(const DFinderCutFlow &cutFlow);

 private:
  void initializeCandidateBranches();
  void initializeEventBranches();
  void initializeCutFlowBranches();

  std::unique_ptr<TFile> m_file;
  TTree *m_candidateTree{nullptr};
  TTree *m_eventTree{nullptr};
  TTree *m_cutFlowTree{nullptr};
  DFinderCandidateRecord m_candidate;
  DFinderEventRecord m_event;
  DFinderCutFlow m_cutFlow;
};

#endif
