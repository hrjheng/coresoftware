#ifndef KSHORTRECONSTRUCTION_H
#define KSHORTRECONSTRUCTION_H

#include <fun4all/SubsysReco.h>

#include <trackbase/ActsTrackingGeometry.h>

#include <Eigen/Dense>

class TFile;
class TH1;
class TH2F;
class TH3D;
class TNtuple;

class ActsGeometry;
class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;

// Struct for Christof's track-pt correction histograms
struct corrMap
{
    TH3D *accpos3D = nullptr;
    TH3D *accneg3D = nullptr;
    TH2F *ptcorr2D = nullptr;
    TH2F *ringcorr2D = nullptr;
    TH2F *etacorr1D = nullptr;
    int verbosity{0};

    void setVerbosity(int verb) { verbosity = verb; }
    void loadCorrMaps(const std::string &file_accposneg3D, const std::string &file_ptcorr2D, const std::string &file_ringcorr2D, const std::string &file_etacorr1D);
    std::tuple<float, float, float> getCorrections(const float &pt, const float &eta, const float &phi, const int &charge) const;
};

class KshortReconstruction : public SubsysReco
{
  public:
    KshortReconstruction(const std::string &name = "KshortReconstruction");
    virtual ~KshortReconstruction() = default;

    int InitRun(PHCompositeNode *topNode) override;
    int process_event(PHCompositeNode *topNode) override;
    int End(PHCompositeNode * /*topNode*/) override;

    void setPtCut(double ptcut) { invariant_pt_cut = ptcut; }
    void setTrackPtCut(double ptcut) { track_pt_cut = ptcut; }
    void setTrackQualityCut(double cut) { _qual_cut = cut; }
    void setPairDCACut(double cut) { pair_dca_cut = cut; }
    void setTrackDCACut(double cut) { track_dca_cut = cut; }
    void setRequireMVTX(bool set) { _require_mvtx = set; }
    void setDecayMass(float decayMassSet) { decaymass = decayMassSet; } //(muons decaymass = 0.1057) (pions = 0.13957) (electron = 0.000511)
    void set_output_file(const std::string &outputfile) { filepath = outputfile; }
    void save_tracks(bool save = true) { m_save_tracks = save; }

    // set charge sign selection ()
    void setChargeSignSelection(std::string sign = "opposite")
    {
        if (sign == "opposite")
        {
            pair_charge_sign = -1;
        }
        else if (sign == "same")
        {
            pair_charge_sign = 1;
        }
        else
        {
            pair_charge_sign = 0;
        }
    }

    // Christof's track-pt correction
    void setCorrectionFiles(const std::string &accposneg3D = "/phenix/u/bogui/data/PtCorr/acc_qinv_hist.root",
                            const std::string &ptcorr2D = "/phenix/u/bogui/data/PtCorr/scale_corr_fullchi2_accqinv.root", //
                            const std::string &ringcorr2D = "/phenix/u/bogui/data/PtCorr/ring_corr.root",                 //
                            const std::string &etacorr1D = "/phenix/u/bogui/data/PtCorr/eta_corr.root"                    //
    )
    {
        m_file_accposneg3D = accposneg3D;
        m_file_ptcorr2D = ptcorr2D;
        m_file_ringcorr2D = ringcorr2D;
        m_file_etacorr1D = etacorr1D;
    }

  private:
    void fillNtp(SvtxTrack *track1, SvtxTrack *track2, Acts::Vector3 dcavals1, Acts::Vector3 dcavals2, Acts::Vector3 pca_rel1, Acts::Vector3 pca_rel2, double pair_dca, double invariantMass, double invariantPt, float invariantPhi, float rapidity, float pseudorapidity, Eigen::Vector3d projected_pos1, Eigen::Vector3d projected_pos2, Eigen::Vector3d projected_mom1, Eigen::Vector3d projected_mom2,
                 Acts::Vector3 pca_rel1_proj, Acts::Vector3 pca_rel2_proj, double pair_dca_proj, unsigned int track1_silicon_cluster_size, unsigned int track2_silicon_cluster_size, unsigned int track1_mvtx_cluster_size, unsigned int track1_mvtx_state_size, unsigned int track1_intt_cluster_size, unsigned int track1_intt_state_size, unsigned int track2_mvtx_cluster_size,
                 unsigned int track2_mvtx_state_size, unsigned int track2_intt_cluster_size, unsigned int track2_intt_state_size, int runNumber, int eventNumber);

    void fillHistogram(Eigen::Vector3d mom1, Eigen::Vector3d mom2, TH1 *massreco, double &invariantMass, double &invariantPt, float &invariantPhi, float &rapidity, float &pseudorapidity);

    // void findPcaTwoTracks(SvtxTrack *track1, SvtxTrack *track2, Acts::Vector3& pca1, Acts::Vector3& pca2, double& dca);
    void findPcaTwoTracks(const Acts::Vector3 &pos1, const Acts::Vector3 &pos2, Acts::Vector3 mom1, Acts::Vector3 mom2, Acts::Vector3 &pca1, Acts::Vector3 &pca2, double &dca) const;

    int getNodes(PHCompositeNode *topNode);

    Acts::Vector3 calculateDca(SvtxTrack *track, const Acts::Vector3 &momentum, Acts::Vector3 position);

    bool projectTrackToCylinder(SvtxTrack *track, double Radius, Eigen::Vector3d &pos, Eigen::Vector3d &mom);
    bool projectTrackToPoint(SvtxTrack *track, Eigen::Vector3d PCA, Eigen::Vector3d &pos, Eigen::Vector3d &mom);

    Acts::Vector3 getVertex(SvtxTrack *track);
    static std::vector<unsigned int> getTrackStates(SvtxTrack *track);

    TNtuple *ntp_reco_info{nullptr};
    ActsGeometry *_tGeometry{nullptr};
    SvtxTrackMap *m_svtxTrackMap{nullptr};
    SvtxVertexMap *m_vertexMap{nullptr};

    std::string filepath{""};
    float decaymass{0.13957}; // pion decay mass
    bool _require_mvtx{true};
    double _qual_cut{1000.0};
    double pair_dca_cut{0.05}; // kshort relative cut 500 microns
    double track_dca_cut{0.01};
    double invariant_pt_cut{0.1};
    double track_pt_cut{0.2};
    int pair_charge_sign{-1}; // 0: same- + opposite-sign, 1: same-sign, -1: opposite-sign
    TFile *fout{nullptr};
    TH1 *recomass{nullptr};

    bool m_save_tracks{false};
    SvtxTrackMap *m_output_trackMap{nullptr};
    std::string m_output_trackMap_node_name{"KshortReconstruction_SvtxTrackMap"};

    // Christof's track-pt correction files, default empty
    std::string m_file_accposneg3D = "";
    std::string m_file_ptcorr2D = "";
    std::string m_file_ringcorr2D = "";
    std::string m_file_etacorr1D = "";

    corrMap m_corrMap;
};

#endif // KSHORTRECONSTRUCTION_H
