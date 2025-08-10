#ifndef ICyDetG4ANALYZER_HXX 
#define ICyDetG4ANALYZER_HXX

#include <vector>

#include "IHandle.hxx"

class TTree;

namespace COMET{
    class IVInputFile;
    class ICOMETEvent;
    class IHitSelection;
}

class ICyDetG4Analyzer{
    public:
        ICyDetG4Analyzer():h_mc_tree(0)
                         ,eventId(-1)
                         ,nHitsFT(0)
                         ,nHitsAT(0)
                         ,maxLayerId(0)
                         ,nTurns(0)
                         ,CLEntrance(0)
                         ,CLExit(0)

                         ,trig(0)
                         ,trig_type(0)
                         ,trig_time(0)
                         ,trig_segId(0)

                         ,first_t(0)
                         ,first_px(0)
                         ,first_py(0)
                         ,first_pz(0)
                         ,first_x(0)
                         ,first_y(0)
                         ,first_z(0)
                         ,first_theta(0)
                         ,ini_px(0)
                         ,ini_py(0)
                         ,ini_pz(0)
                         ,ini_x(0)
                         ,ini_y(0)
                         ,ini_z(0)
                         ,ini_t(0)

                         ,cth_nHits(0)
                         ,cth_segId(0)
                         ,cth_hodId(0)
                         ,cth_cryType(0)
                         ,cth_pdgId(0)
                         ,cth_trackId(0)
                         ,cth_stepL(0)
                         ,cth_edep(0)
                         ,cth_beta(0)
                         ,cth_trig(0)
                         ,cth_t(0)
                         ,cth_px(0)
                         ,cth_py(0)
                         ,cth_pz(0)
                         ,cth_x(0)
                         ,cth_y(0)
                         ,cth_z(0)

                         ,cthmc_nHits(0)
                         ,cthmc_nCon(0)
                         ,cthmc_segId(0)
                         ,cthmc_hodId(0)
                         ,cthmc_cryType(0)
                         ,cthmc_pdgId(0)
                         ,cthmc_trackId(0)
                         ,cthmc_chargeP(0)
                         ,cthmc_chargeS(0)
                         ,cthmc_stepL(0)
                         ,cthmc_edep(0)
                         ,cthmc_beta(0)
                         ,cthmc_trig(0)
                         ,cthmc_tT(0)
                         ,cthmc_tP(0)
                         ,cthmc_px(0)
                         ,cthmc_py(0)
                         ,cthmc_pz(0)
                         ,cthmc_x(0)
                         ,cthmc_y(0)
                         ,cthmc_z(0)

                         ,cthTrig_nHits(0)
                         ,cthTrig_segId(0)
                         ,cthTrig_hodId(0)
                         ,cthTrig_cryType(0)
                         ,cthTrig_t(0)
                         ,cthTrig_ADC(0)

                         ,cdc_nHits(0)
                         ,cdc_cellId(0)
                         ,cdc_layerId(0)
                         ,cdc_pdgId(0)
                         ,cdc_trackId(0)
                         ,cdc_hittype(0)
                         ,cdc_edep(0)
                         ,cdc_trig(0)
                         ,cdc_length(0)
                         ,cdc_driftD(0)
                         ,cdc_pulseStartT(0)
                         ,cdc_pulseStopT(0)
                         ,cdc_IO(0)
                         ,cdc_nPairs(0)
                         ,cdc_DOCA(0)
                         ,cdc_t(0)
                         ,cdc_first_t(0)
                         ,cdc_last_t(0)
                         ,cdc_px(0)
                         ,cdc_py(0)
                         ,cdc_pz(0)
                         ,cdc_first_px(0)
                         ,cdc_first_py(0)
                         ,cdc_first_pz(0)
                         ,cdc_last_px(0)
                         ,cdc_last_py(0)
                         ,cdc_last_pz(0)
                         ,cdc_x(0)
                         ,cdc_y(0)
                         ,cdc_z(0)
                         ,cdc_first_x(0)
                         ,cdc_first_y(0)
                         ,cdc_first_z(0)
                         ,cdc_last_x(0)
                         ,cdc_last_y(0)
                         ,cdc_last_z(0)
                         ,cdc_wx(0)
                         ,cdc_wy(0)
                         ,cdc_wz(0)

                         ,fDisableMCPrimary(false)
                         ,fDisableMCTrigger(false)
                         ,fDisableMCCTHHits(false)
                         ,fDisableMCCDCHits(false)

                         ,fProtonMass(0)
                         ,fNeutronMass(0)
                         {}
        virtual ~ICyDetG4Analyzer()
        {}

        int Process(COMET::ICOMETEvent& ievent);
        bool SetOptions(std::map<std::string,std::string> options);
        int Init();
        int BeginFile(COMET::IVInputFile *const);
        int Finalize();

    private:
        void PrepareOutputs();
        int InitTree();
        void Read_MC_info(COMET::ICOMETEvent& event);
        void Read_MC_mom(COMET::ICOMETEvent& event);
        void Read_MC_trigger(COMET::ICOMETEvent& event);
        void Read_MC_cthHits(COMET::ICOMETEvent& event);
        void Read_MC_cdcHits(COMET::ICOMETEvent& event);
        void Analysis();
        void Delete();
        //// return particle mass in MeV
        double GetMass(const int& pdgCode,double& charge);
        double GetMass(const std::string& partName);

        TTree* h_mc_tree;

        int  eventId;

        bool fDisableMCPrimary;
        bool fDisableMCTrigger;
        bool fDisableMCCTHHits;
        bool fDisableMCCDCHits;

        // mass of p/n
        double fProtonMass;
        double fNeutronMass;

        // about output values
        int                   nHitsFT; // number signal hits of the first turn in CDC
        int                   nHitsAT; // number signal hits in CDC
        int                   maxLayerId; // maximum layerId hit by signal track in CDC
        int                   nTurns; // number of turns
        int                   CLEntrance; // maximum continuous number of layers hit by signal track upon entering CDC
        int                   CLExit; // maximum continuous number of layers hit by signal track upon leaving CDC

        int                   trig; // whether CTH is triggered
        int                   trig_type; // CTH hit pattern
        double                trig_time; // hit time in CTH by signal track [ns]
        int                   trig_segId; // volume Id in CTH hit by signal track

        double                first_t; // tof of the signal track upon entering CDC [ns]
        double                first_px; // momentum of the signal track upon entering CDC [GeV/c]
        double                first_py; // momentum of the signal track upon entering CDC [GeV/c]
        double                first_pz; // momentum of the signal track upon entering CDC [GeV/c]
        double                first_x; // position of the signal track upon entering CDC [cm]
        double                first_y; // position of the signal track upon entering CDC [cm]
        double                first_z; // position of the signal track upon entering CDC [cm]
        double                first_theta; // angle of the track with z axis upon entering CDC [rad]
        double                ini_px; // momentum of the signal track upon generation [GeV/c]
        double                ini_py; // momentum of the signal track upon generation [GeV/c]
        double                ini_pz; // momentum of the signal track upon generation [GeV/c]
        double                ini_x; // position of the signal track upon generation [cm]
        double                ini_y; // position of the signal track upon generation [cm]
        double                ini_z; // position of the signal track upon generation [cm]
        int                   ini_t; // time of the signal track upon generation [ns];

        int                   cth_nHits; // number of hits in CTH
        std::vector<int>     *cth_segId; // segment Id [0-N]: 0 from phi=0
        std::vector<int>     *cth_hodId; // hodoscope Id [0-N]: 0 upstream; 1 downstream
        std::vector<int>     *cth_cryType; // 0: scintillator; 1: cherenkov counter; 2: light guide
        std::vector<int>     *cth_pdgId; // PDG encoding
        std::vector<int>     *cth_trackId; // 1: signal; 
        std::vector<double>  *cth_stepL; // step length
        std::vector<double>  *cth_edep; // energy deposit
        std::vector<double>  *cth_beta; // beta of the track
        std::vector<int>     *cth_trig; // triggerd TDC
        std::vector<double>  *cth_t; // time when the track hit this cell [ns]
        std::vector<double>  *cth_px; // momentum of the track in this cell [MeV/c]
        std::vector<double>  *cth_py; // momentum of the track in this cell [MeV/c]
        std::vector<double>  *cth_pz; // momentum of the track in this cell [MeV/c]
        std::vector<double>  *cth_x; // position of the track point in this cell closest to this sense wire [cm]
        std::vector<double>  *cth_y; // position of the track point in this cell closest to this sense wire [cm]
        std::vector<double>  *cth_z; // position of the track point in this cell closest to this sense wire [cm]

        int                   cthmc_nHits; // number of hits in CTH
        std::vector<int>     *cthmc_nCon; // number of contributors
        std::vector<int>     *cthmc_segId; // segment Id [0-N]: 0 from phi=0
        std::vector<int>     *cthmc_hodId; // hodoscope Id [0-N]: 0 upstream; 1 downstream
        std::vector<int>     *cthmc_cryType; // 0: scintillator; 1: cherenkov counter; 2: light guide
        std::vector<int>     *cthmc_pdgId; // PDG encoding
        std::vector<int>     *cthmc_trackId; // 1: signal; 
        std::vector<double>  *cthmc_chargeP; // charge at peak [eplus]
        std::vector<double>  *cthmc_chargeS; // charge sum [eplus]
        std::vector<double>  *cthmc_stepL; // step length
        std::vector<double>  *cthmc_edep; // energy deposit
        std::vector<double>  *cthmc_beta; // beta of the track
        std::vector<int>     *cthmc_trig; // triggerd TDC
        std::vector<double>  *cthmc_tT; // time at the tick over threshold  [ns]
        std::vector<double>  *cthmc_tP; // time at the peak [ns]
        std::vector<double>  *cthmc_px; // momentum of the track in this cell [MeV/c]
        std::vector<double>  *cthmc_py; // momentum of the track in this cell [MeV/c]
        std::vector<double>  *cthmc_pz; // momentum of the track in this cell [MeV/c]
        std::vector<double>  *cthmc_x; // position of the track point in this cell closest to this sense wire [cm]
        std::vector<double>  *cthmc_y; // position of the track point in this cell closest to this sense wire [cm]
        std::vector<double>  *cthmc_z; // position of the track point in this cell closest to this sense wire [cm]

        int                   cthTrig_nHits;
        std::vector<int>     *cthTrig_segId;
        std::vector<int>     *cthTrig_hodId;
        std::vector<int>     *cthTrig_cryType;
        std::vector<int>     *cthTrig_t;
        std::vector<int>     *cthTrig_ADC;

        int                   cdc_nHits; // number of all CDC hits
        std::vector<int>     *cdc_cellId; // cell Id [0-N]
        std::vector<int>     *cdc_layerId; // layer Id [1-18] (0 and 19 are guard layers)
        std::vector<int>     *cdc_pdgId; // PDG encoding
        std::vector<int>     *cdc_trackId; // 1: signal; 
        std::vector<int>     *cdc_hittype; // 1: cell; 2: hit field wire; 3: hit sense wire
        std::vector<double>  *cdc_edep; // energy deposit [GeV]
        std::vector<int>     *cdc_trig; // triggerd TDC
        std::vector<double>  *cdc_length; // step length in this cell [cm]
        std::vector<double>  *cdc_driftD; // drift distance in this cell[cm]
        std::vector<double>  *cdc_pulseStartT; // time when the waveform starts [ns]
        std::vector<double>  *cdc_pulseStopT; // time when the waveform stops [ns]
        std::vector<double>  *cdc_IO; // going in or out?
        std::vector<int>     *cdc_nPairs; // number of electron-ion pairs generated during this cell
        std::vector<double>  *cdc_DOCA; // drift distance in this cell[cm]
        std::vector<double>  *cdc_t; // time when the track hit this cell [ns]
        std::vector<double>  *cdc_first_t; // time when the track entered this cell [ns]
        std::vector<double>  *cdc_last_t; // time when the track lasted this cell [ns]
        std::vector<double>  *cdc_px; // momentum of the track in this cell [MeV/c]
        std::vector<double>  *cdc_py; // momentum of the track in this cell [MeV/c]
        std::vector<double>  *cdc_pz; // momentum of the track in this cell [MeV/c]
        std::vector<double>  *cdc_first_px; // momentum of the track when entering this cell [MeV/c]
        std::vector<double>  *cdc_first_py; // momentum of the track when entering this cell [MeV/c]
        std::vector<double>  *cdc_first_pz; // momentum of the track when entering this cell [MeV/c]
        std::vector<double>  *cdc_last_px; // momentum of the track when lasting this cell [MeV/c]
        std::vector<double>  *cdc_last_py; // momentum of the track when lasting this cell [MeV/c]
        std::vector<double>  *cdc_last_pz; // momentum of the track when lasting this cell [MeV/c]
        std::vector<double>  *cdc_x; // position of the track point in this cell closest to this sense wire [cm]
        std::vector<double>  *cdc_y; // position of the track point in this cell closest to this sense wire [cm]
        std::vector<double>  *cdc_z; // position of the track point in this cell closest to this sense wire [cm]
        std::vector<double>  *cdc_first_x; // position of the track point in this cell as the track enters [cm]
        std::vector<double>  *cdc_first_y; // position of the track point in this cell as the track enters [cm]
        std::vector<double>  *cdc_first_z; // position of the track point in this cell as the track enters [cm]
        std::vector<double>  *cdc_last_x; // position of the track point in this cell as the track lasts [cm]
        std::vector<double>  *cdc_last_y; // position of the track point in this cell as the track lasts [cm]
        std::vector<double>  *cdc_last_z; // position of the track point in this cell as the track lasts [cm]
        std::vector<double>  *cdc_wx; // position of the point on this sense wire closest to the track in this cell [cm]
        std::vector<double>  *cdc_wy; // position of the point on this sense wire closest to the track in this cell [cm]
        std::vector<double>  *cdc_wz; // position of the point on this sense wire closest to the track in this cell [cm]
};
#endif
