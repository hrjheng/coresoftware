#ifndef MBD_MBDEVENT_H
#define MBD_MBDEVENT_H

#include "MbdDefs.h"
#include "MbdSig.h"

#include <TFile.h>
#include <TTree.h>

#ifndef ONLINE
#include <fun4all/Fun4AllBase.h>
#endif

#include <vector>
#include <limits>

class PHCompositeNode;
class Event;
class MbdPmtContainer;
class MbdOut;
class MbdCalib;
class MbdGeom;
class CDBUtils;
class TF1;
class TCanvas;
#ifndef ONLINE
class CaloPacketContainer;
class Gl1Packet;
class PHG4TruthInfoContainer;
class PHG4VtxPoint;
#endif

class MbdEvent
{
 public:
  MbdEvent(const int cal_pass = 0);
  virtual ~MbdEvent();

  int SetRawData(Event *event, MbdPmtContainer *bbcpmts);
#ifndef ONLINE
  int SetRawData(CaloPacketContainer *mbdraw, MbdPmtContainer *bbcpmts, Gl1Packet *gl1raw);
#endif
  void PostProcessChannels(MbdPmtContainer *bbcpmts);
  int Calculate(MbdPmtContainer *bbcpmts, MbdOut *bbcout, PHCompositeNode *topNode = nullptr);
  int InitRun();
  int End();
  void Clear();

  void SetSim(const int s) { _simflag = s; }

  float get_bbcz() { return m_bbcz; }
  float get_bbczerr() { return m_bbczerr; }
  float get_bbct0() { return m_bbct0; }
  float get_bbct0err() { return m_bbct0err; }

  int   get_bbcn(const int iarm) { return m_bbcn[iarm]; }
  float get_bbcq(const int iarm) { return m_bbcq[iarm]; }
  float get_bbct(const int iarm) { return m_bbct[iarm]; }
  float get_bbcte(const int iarm) { return m_bbcte[iarm]; }

  int   get_pmtq(const int ipmt) { return m_pmtq[ipmt]; }
  float get_pmttt(const int ipmt) { return m_pmttt[ipmt]; }
  float get_pmttq(const int ipmt) { return m_pmttq[ipmt]; }

  int  get_EventNumber(void) const { return m_evt; }
  void set_EventNumber(int ievt) { m_evt = ievt; }

  void set_debug(const int d) { _debug = d; }

  MbdSig *GetSig(const int ipmt) { return &_mbdsig[ipmt]; }

  MbdCalib *GetCalib() { return _mbdcal; }

  int FillSampMaxCalib();

  int  calib_is_done() { return _calib_done; }

  int ProcessRawPackets(MbdPmtContainer *bbcpmts);

  int  Verbosity() { return _verbose; }
  void Verbosity(const int v) { _verbose = v; }

 private:
  static const int NCHPERPKT = 128;

  MbdGeom *_mbdgeom{nullptr};
  MbdCalib *_mbdcal{nullptr};

  int Read_Charge_Calib(const std::string &gainfname);
  int Read_TQ_T0_Offsets(const std::string &t0cal_fname);
  int Read_TQ_CLK_Offsets(const std::string &t0cal_fname);
  int Read_TT_CLK_Offsets(const std::string &t0cal_fname);
  //int DoQuickClockOffsetCalib();

  bool isbadtch(const int ipmtch);

  // Debugging variables
  int _debug{0};
#ifndef ONLINE
  PHG4TruthInfoContainer* _truth_container {nullptr};
  PHG4VtxPoint* _vtxp {nullptr};
  PHG4VtxPoint* GetPrimaryVtx(PHCompositeNode *topNode);
#endif
  int epmt[2]{-1, -1};  // pmt of earliest time
  // int lpmt[2] {-1,-1};        // pmt of latest time
  double tepmt[2]{1e9, 1e9};    // earliest time
  double tlpmt[2]{-1e9, -1e9};  // latest time
  void ReadSyncFile(const char *fname = "SYNC_INTTMBD.root");

  float gaincorr[MbdDefs::MBD_N_PMT]{};       // gain corrections
  float tq_t0_offsets[MbdDefs::MBD_N_PMT]{};  // t0 offsets in charge channels
  float tq_clk_offsets[MbdDefs::MBD_N_PMT]{};
  float tt_clk_offsets[MbdDefs::MBD_N_PMT]{};

  // float bz_offset{0.};

  int _verbose{0};
  int _runnum{0};
  int _simflag{0};
  int _nsamples{31};
  int _calib_done{0}; 
  unsigned int _no_sampmax{0};      //! sampmax calib doesn't exist
  int _is_online{0};                //! for OnlMon

  // alignment data
  Int_t   m_evt{0};
  Short_t m_clk{0};
  Short_t m_femclk{0};
  UInt_t  m_xmitclocks[2]{};     // [ipkt]
  UInt_t  m_femclocks[2][2]{};   // [ipkt][iadc]

  // raw data
  Float_t m_adc[MbdDefs::MBD_N_FEECH][MbdDefs::MAX_SAMPLES]{};   // raw waveform, adc values
  Float_t m_samp[MbdDefs::MBD_N_FEECH][MbdDefs::MAX_SAMPLES]{};  // raw waveform, sample values
  Float_t m_ampl[MbdDefs::MBD_N_FEECH]{};                        // raw amplitude

  std::vector<MbdSig> _mbdsig;

  Float_t m_pmtq[MbdDefs::MBD_N_PMT]{};   // npe in each arm
  Float_t m_pmttt[MbdDefs::MBD_N_PMT]{};  // time in each arm
  Float_t m_pmttq[MbdDefs::MBD_N_PMT]{};  // time in each arm

  int do_templatefit{1};

  // output data
  Short_t m_bbcn[2]{};                                            // num hits for each arm (north and south)
  Float_t m_bbcq[2]{};                                            // total charge (currently npe) in each arm
  Float_t m_bbct[2]{};                                            // time in arm
  Float_t m_bbcte[2]{};                                           // earliest hit time in arm
  Float_t m_bbctl[2]{};                                           // latest hit time in arm
  Float_t m_bbcz{std::numeric_limits<Float_t>::quiet_NaN()};      // z-vertex
  Float_t m_bbczerr{std::numeric_limits<Float_t>::quiet_NaN()};   // z-vertex error
  Float_t m_bbct0{std::numeric_limits<Float_t>::quiet_NaN()};     // start time
  Float_t m_bbct0err{std::numeric_limits<Float_t>::quiet_NaN()};  // start time error
  Float_t _tres{std::numeric_limits<Float_t>::quiet_NaN()};       // time resolution of one channel

  TH1 *hevt_bbct[2]{};  // time in each bbc, per event
  TF1 *gausfit[2]{nullptr, nullptr};

  float TRIG_SAMP[16]{};  // [board]

  // Calibration Data
  int _calpass{0};
  TString _caldir;
  //std::string _caldir;

  // sampmax and other calib stuff
  int CalcSampMaxCalib();
  std::unique_ptr<TFile> _calpass1_tfile{nullptr};
  std::unique_ptr<TFile> _calpass2_tfile{nullptr};
  TH1 *h_smax[256]{};     // [feech], max sample in event
  TH2 *h2_smax[2]{};      // [0 == time ch, 1 == chg ch], max sample in evt vs ch
  TH2 *h2_wave[2]{};      // [0 == time ch, 1 == chg ch], all samples in evt vs ch
  TH2 *h2_trange_raw{};   // raw tdc at maxsamp vs ch
  TH2 *h2_trange{};       // subtracted tdc at maxsamp vs ch
  //TH1 *h_trange[2]{};     // subtracted tdc at maxsamp, [S/N]

  // pedestals (hists are in MbdSig)
  int CalcPedCalib();

  //
  void ClusterEarliest(std::vector<float> &times, double& mean, double& rms, double& rmin, double& rmax);
 
  // debug stuff
  TCanvas *ac{nullptr};  // for plots used during debugging
  void PlotDebug();
  std::unique_ptr<TFile> _synctfile{nullptr};
  TTree *_syncttree{nullptr};
  Double_t _refz{ std::numeric_limits<double>::quiet_NaN() };
  std::vector<Int_t> bbevt;
  std::vector<UShort_t> bbclk;
  std::vector<Float_t> mybbz;
  std::vector<Long64_t> bco;
  std::vector<Double_t> intz;
  std::vector<Double_t> bbz;
};

#endif /* MBD_MBDEVENT_H */
