#include <map>
#include <math.h>
#include <TTree.h>
#include <TRandom.h>

#include <TGeoNode.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include "IG4PrimaryVertex.hxx"
#include "ICyDetG4Analyzer.hxx"
#include "IVInputFile.hxx"
#include "ISIMG4Header.hxx"
#include "IMCHit.hxx"
#include "IHandle.hxx"
#include "ICDCGeom.hxx"
#include "ICOMETEvent.hxx"
#include "IGeomInfo.hxx"
#include "IGeomTools.hxx"
#include "IOARuntimeParameters.hxx"
#include "IDigitContainer.hxx"
#include "ICTHmcDigit.hxx"
#include "ICTHChannelId.hxx"
#include "IGeomIdManager.hxx"
#include "IGeometryDatabase.hxx"
#include "ICTHGeomId.hxx"
#include "IG4HitSegment.hxx"
#include "IG4HitGas.hxx"
#include "ITrigger.hxx"

int ICyDetG4Analyzer::Process(COMET::ICOMETEvent& event){
    PrepareOutputs();
    Read_MC_info(event);
    h_mc_tree->Fill();
    return 1;
}

/// Read information from Monte Carlo simulation
void ICyDetG4Analyzer::Read_MC_info(COMET::ICOMETEvent& event){
    eventId = event.GetEventId();
    if (!fDisableMCPrimary) Read_MC_mom(event);
    if (!fDisableMCTrigger) Read_MC_trigger(event);
    if (!fDisableMCCTHHits) Read_MC_cthHits(event);
    if (!fDisableMCCDCHits) Read_MC_cdcHits(event);
    return;
}

/// Read Monte Carlo CTH hits
void ICyDetG4Analyzer::Read_MC_cthHits(COMET::ICOMETEvent& event){
    COMETInfo("Started reading MC CTH hits");
    COMET::IHandle<COMET::IG4HitContainer> g4Hits = event.Get<COMET::IG4HitContainer>("truth/g4Hits/CTH");
    if (!g4Hits){
        COMETInfo("Cannot find g4Hit for CTH");
        return;
    }
    COMETInfo("    Got "<<g4Hits->size()<<" hits");
    COMET::IGeometryId pre_geomId;
    int pre_trackId = -1;
    for (COMET::IG4HitContainer::const_iterator hitItr = g4Hits->begin(); hitItr != g4Hits->end(); ++hitItr) {
        COMETInfo("   Starting");
        COMET::IG4HitSegment * g4Seg = dynamic_cast<COMET::IG4HitSegment*>(*hitItr);
        if (!g4Seg){
            COMETError("Cannot convert g4Hit to g4Seg");
            return;
        }

        int pdgEncoding = g4Seg->GetPDGEncoding();
        int trackId     = g4Seg->GetPrimaryId(); // primary Id should be same as trakc Id now..
        double stepL    = g4Seg->GetTrackLength();
        double edep     = g4Seg->GetEnergyDeposit();
        double tStart   = g4Seg->GetStartT();
        double xStart   = g4Seg->GetStartX();
        double yStart   = g4Seg->GetStartY();
        double zStart   = g4Seg->GetStartZ();
        double xStop    = g4Seg->GetStopX();
        double yStop    = g4Seg->GetStopY();
        double zStop    = g4Seg->GetStopZ();
        double px       = g4Seg->GetMomentumX();
        double py       = g4Seg->GetMomentumY();
        double pz       = g4Seg->GetMomentumZ();
        double pa       = g4Seg->GetMomentumMag();
        double charge   = 0;
        double mass     = this->GetMass(pdgEncoding,charge);
        double beta     = pa/sqrt(mass*mass+pa*pa);
        //COMETInfo("Particle info: "<<trackId<<", "<<pdgEncoding<< ", "<<mass<<", "<<charge<<", "<<TMath::Sqrt(px*px+py*py+pz*pz)<<", "<<beta);
        if (charge==0) beta = 0.0;

        TVector3 posStartLocal = TVector3(xStart,yStart,zStart);
        TVector3 momLocal      = TVector3(px,py,pz);

        COMET::IGeometryId geomId;
        if (!COMET::IOADatabase::Get().GeomId().GetGeometryId((xStart+xStop)/2,(yStart+yStop)/2,(zStart+zStop)/2, geomId)){
          COMETInfo("    cannot find geomId? "<<(xStart+xStop)/2<<" "<<(yStart+yStop)/2<<" "<<(zStart+zStop)/2);
          TGeoManager* geometry=COMET::IOADatabase::Get().GetGeometry();
          COMETInfo("    Got geometry @ "<<(void*)geometry);
          TString name = "NONE";
          TGeoNode * node = COMET::GeomTools::GetNode(geometry,TVector3((xStart+xStop)/2,(yStart+yStop)/2,(zStart+zStop)/2));
          if (node) name = node->GetName();
          COMETInfo("Cannot find geomId at "<<(xStart+xStop)/2<<" "<<(yStart+yStop)/2<<" "<<(zStart+zStop)/2<<" in node "<<name);
          continue;
        }
        int module = COMET::GeomId::CTH::GetModule(geomId);        /// upstream/downstream
        int scint  = COMET::GeomId::CTH::GetScintillator(geomId);  /// sinti layer or not
        int light  = COMET::GeomId::CTH::GetLightGuide(geomId);    /// light guide or not
        int number = COMET::GeomId::CTH::GetSegmentId(geomId);     /// id in phi segment
        if (pre_trackId!=trackId||pre_geomId!=geomId){ // create new hits for different geomId or different track Id
            COMETInfo("    New channel! "<<module<<", "<<scint<<", "<<number);
            cth_nHits++;
            cth_trig->push_back(1); // FIXME: for now take them as all triggered
            cth_stepL->push_back(stepL/unit::mm);
            cth_edep->push_back(edep/unit::MeV);
            cth_beta->push_back(beta);
            cth_segId->push_back(number);
            cth_hodId->push_back(module);
            cth_cryType->push_back(scint+light*10);
            cth_pdgId->push_back(pdgEncoding);
            cth_trackId->push_back(trackId);
            cth_t->push_back(tStart/unit::ns);
            cth_x->push_back(posStartLocal.X()/unit::mm);
            cth_y->push_back(posStartLocal.Y()/unit::mm);
            cth_z->push_back(posStartLocal.Z()/unit::mm);
            cth_px->push_back(momLocal.Px()/unit::MeV);
            cth_py->push_back(momLocal.Py()/unit::MeV);
            cth_pz->push_back(momLocal.Pz()/unit::MeV);
        }
        else{
            COMETInfo("    Same channel! "<<module<<", "<<scint<<", "<<number);
            cth_stepL->at(cth_nHits-1)+=stepL/unit::mm;
            cth_edep->at(cth_nHits-1)+=edep/unit::MeV;
        }
        pre_trackId = trackId;
        pre_geomId  = geomId;
        COMETInfo("   Finished");
    }

    COMET::IHandle<COMET::IHitSelection> hits = event.GetHitSelection("mcCTH");
    if(hits){
        COMETInfo("  We have mcCTH???");
        COMET::IHitSelection::const_iterator hitItr;
        // FIXME: now relying on the fact that hits from the same channel comes together
        int pre_number = -1;
        int pre_scint  = -1;
        int pre_module = -1;
        for (hitItr = hits->begin(); hitItr != hits->end(); ++hitItr){
            COMET::IHandle<COMET::IMCHit> hit = (*hitItr);
            const COMET::ICTHChannelId channelId = static_cast<const COMET::ICTHChannelId> (hit->GetChannelID());
            int module = channelId.GetModule();
            int scint  = channelId.GetScint();
            int light  = channelId.GetLightGuide();
            int number = channelId.GetSegmentId();
            double time = hit->GetTime();
            double charge = hit->GetCharge();
            COMET::IMCHit::ContributorContainer contributors = hit->GetContributors();
            int    trackId = 0;
            int    pdgEncoding = 0;
            double stepL = 0;
            double edep = 0;
            double tStart = 0;
            double xStart = 0;
            double yStart = 0;
            double zStart = 0;
            double xStop = 0;
            double yStop = 0;
            double zStop = 0;
            double px = 0;
            double py = 0;
            double pz = 0;
            double pa = 0;
            double pcharge = 0; /// charge of particle
            double mass = 0;
            double beta = 0;
            int    nContributors = contributors.size();
            for(int iCon = 0; iCon<nContributors; iCon++){
                COMET::IG4HitSegment * g4Seg = dynamic_cast<COMET::IG4HitSegment*> (contributors[iCon]);
                if (!g4Seg) continue;
                if (iCon==0){ 
                    trackId = g4Seg->GetPrimaryId();
                    pdgEncoding = g4Seg->GetPDGEncoding();
                    mass        = this->GetMass(pdgEncoding,pcharge);
                    tStart = g4Seg->GetStartT();
                    xStart = g4Seg->GetStartX();
                    yStart = g4Seg->GetStartY();
                    zStart = g4Seg->GetStartZ();
                    px = g4Seg->GetMomentumX();
                    py = g4Seg->GetMomentumY();
                    pz = g4Seg->GetMomentumZ();
                    pa = g4Seg->GetMomentumMag();
                    beta = pa/sqrt(mass*mass+pa*pa);
                    if (pcharge==0) beta = 0.0;
                }
                else if (iCon==nContributors-1){
                    xStop = g4Seg->GetStopX();
                    yStop = g4Seg->GetStopY();
                    zStop = g4Seg->GetStopZ();
                }
                stepL += g4Seg->GetTrackLength();
                edep  += g4Seg->GetEnergyDeposit();
            }

            TVector3 posStartLocal = TVector3(xStart,yStart,zStart);
            TVector3 momLocal = TVector3(px,py,pz);
            if (charge<1*unit::eplus){ // FIXME: taking 1 as the threshold at this moment
                pre_number = -1;
                pre_scint  = -1;
                pre_module = -1;
                continue;
            }
            else if (pre_module!=module||pre_scint!=scint||pre_number!=number){ // it's a new channel
                cthmc_nHits++;
                cthmc_nCon->push_back(nContributors);
                cthmc_trig->push_back(1); // FIXME: for now take them as all triggered
                cthmc_chargeP->push_back(charge/unit::eplus);
                cthmc_chargeS->push_back(charge/unit::eplus);
                cthmc_beta->push_back(beta);
                cthmc_edep->push_back(edep/unit::MeV);
                cthmc_stepL->push_back(stepL/unit::mm);
                cthmc_tT->push_back(time/unit::ns);
                cthmc_tP->push_back(time/unit::ns); // What is tP??
                cthmc_segId->push_back(number);
                cthmc_hodId->push_back(module);
                cthmc_cryType->push_back(scint+light*10);
                cthmc_pdgId->push_back(pdgEncoding);
                cthmc_trackId->push_back(trackId);
                cthmc_x->push_back(posStartLocal.X()/unit::mm);
                cthmc_y->push_back(posStartLocal.Y()/unit::mm);
                cthmc_z->push_back(posStartLocal.Z()/unit::mm);
                cthmc_px->push_back(momLocal.Px()/unit::MeV);
                cthmc_py->push_back(momLocal.Py()/unit::MeV);
                cthmc_pz->push_back(momLocal.Pz()/unit::MeV);
            }
            else{
                if (charge>cthmc_chargeP->at(cthmc_nHits-1)){
                    cthmc_chargeP->at(cthmc_nHits-1)=charge;
                    cthmc_tP->at(cthmc_nHits-1)=time/unit::nm; ///???
                }
                cthmc_chargeS->at(cthmc_nHits-1)+=charge;
            }
            pre_number = number;
            pre_scint  = scint;
            pre_module = module;
        }
    }
    else{
        COMETInfo("There is no mc cth hit in this event!");
    }

    COMET::IHandle<COMET::IDigitContainer> digits = event.GetDigits("CTHtrig");
    if (digits){
        COMET::IDigitContainer::const_iterator digItr;
        for (digItr = digits->begin(); digItr != digits->end(); ++digItr){
            COMET::ICTHmcDigit * digit = dynamic_cast<COMET::ICTHmcDigit*> (*digItr);
            const COMET::ICTHChannelId channelId = static_cast<const COMET::ICTHChannelId> (digit->GetChannelId());
            int module =channelId.GetModule();
            int scint  =channelId.GetScint();
            int light  =channelId.GetLightGuide();
            int number =channelId.GetSegmentId();
            const std::vector<UShort_t> & adcs = digit->GetADCs();
            COMET::ITrigger trigger = *digit->GetTrigger();
            if (trigger.GetTriggerStatus(channelId)!=COMET::ITrigger::kReadout) continue;
            int trigTick = trigger.GetTriggerTick(); ///
            if (trigTick<0 || trigTick>adcs.size())continue;
            cthTrig_segId->push_back(number);
            cthTrig_hodId->push_back(module);
            cthTrig_cryType->push_back(scint+light*10);
            cthTrig_t->push_back(trigTick);
            cthTrig_ADC->push_back((int)adcs[trigTick]);
            cthTrig_nHits++;
        }
    }
    else{
        COMETInfo("There is no mc cth digit in this event!");
    }
}

/// Read Monte Carlo CDC hits
void ICyDetG4Analyzer::Read_MC_cdcHits(COMET::ICOMETEvent& event){
    COMETInfo("Start reading MC CDC");
    COMET::IHandle<COMET::IG4HitContainer> hits = event.Get<COMET::IG4HitContainer>("truth/g4Hits/CDC");
    if(hits){
        const COMET::ICDCGeom& geoCDC = COMET::IGeomInfo::CDC();
        int nLayers = geoCDC.GetNumberOfLayers()+2;
        std::vector<bool> layerHitEnter(nLayers,false);
        std::vector<bool> layerHitExit(nLayers,false);
        first_t = 1e9;
        COMETInfo("    Got "<<hits->size()<<" hits");
        for (COMET::IG4HitContainer::const_iterator hitItr = hits->begin(); hitItr != hits->end(); ++hitItr) {
            COMET::IG4HitGas * this_mchit = dynamic_cast<COMET::IG4HitGas*>(*hitItr);
            if (!this_mchit){
                COMETError("Cannot convert hit from truth/g4Hits/CDC into a IG4HitGas");
                return;
            }
            int channelId = -1;
            COMET::IGeomInfo::CDC().GlobalPositionToChannel(
                    this_mchit->GetPosition(), channelId);
            int layerId = geoCDC.GetLayer(channelId);
            int cellId = geoCDC.GetCellId(channelId);
            int wireId = geoCDC.ChannelToWire(channelId);

            cdc_nHits++;

            TVector3 posLocal = COMET::IGeomInfo::Get().CDC().GlobalPosition_to_LocalPosition(this_mchit->GetPosition());
            TVector3 wirePosLocal = COMET::IGeomInfo::Get().CDC().GlobalPosition_to_LocalPosition(this_mchit->GetPOCAOnWire());
            TVector3 momLocal = COMET::IGeomInfo::Get().CDC().GlobalMomentum_to_LocalMomentum(this_mchit->GetMomentum());

            TVector3 firstPosLocal = COMET::IGeomInfo::Get().CDC().GlobalPosition_to_LocalPosition(this_mchit->GetFirstPosition());
            TVector3 lastPosLocal = COMET::IGeomInfo::Get().CDC().GlobalPosition_to_LocalPosition(this_mchit->GetLastPosition());
            TVector3 firstMomLocal = COMET::IGeomInfo::Get().CDC().GlobalMomentum_to_LocalMomentum(this_mchit->GetFirstMomentum());
            TVector3 lastMomLocal = COMET::IGeomInfo::Get().CDC().GlobalMomentum_to_LocalMomentum(this_mchit->GetLastMomentum());


            int PDGEncoding = this_mchit->GetPDGEncoding();
            int trackId = this_mchit->GetPrimaryId();
            double mchit_time = this_mchit->GetPosT();
            double mchit_first_time = this_mchit->GetFirstTime();
            double mchit_last_time = this_mchit->GetLastTime();
            double mchit_x = posLocal.X();
            double mchit_y = posLocal.Y();
            double mchit_z = posLocal.Z();
            double mchit_px = momLocal.X();
            double mchit_py = momLocal.Y();
            double mchit_pz = momLocal.Z();
            double mchit_DriftD = this_mchit->GetDOCA(); // fake driftD since we don't have that drift simulation in Geant
            double mchit_DOCA = this_mchit->GetDOCA();
            double mchit_edep = this_mchit->GetEnergyDeposit();
            double mchit_stepL = this_mchit->GetLength();
            double mchit_wx = wirePosLocal.X();
            double mchit_wy = wirePosLocal.Y();
            double mchit_wz = wirePosLocal.Z();

            if(mchit_time < first_t){
                first_t = mchit_time;
                first_x = posLocal.X();
                first_y = posLocal.Y();
                first_z = posLocal.Z();
                first_px = momLocal.Px();
                first_py = momLocal.Py();
                first_pz = momLocal.Pz();
                first_theta = atan2(sqrt(first_pz*first_pz+first_py*first_py),first_px);
            }

            cdc_layerId->push_back(layerId);
            cdc_cellId->push_back(cellId);
            cdc_pdgId->push_back(PDGEncoding);
            cdc_trackId->push_back(trackId);
            cdc_edep->push_back(mchit_edep/unit::MeV);
            cdc_trig->push_back(1); // Assuming everyone is triggered; Will be set in Read_DT_Hits()
            cdc_length->push_back(mchit_stepL/unit::mm);
            cdc_pulseStartT->push_back(1); // Default value is 1; Will be set in Read_DT_Hits()
            cdc_pulseStopT->push_back(1); // Default value is 1; Will be set in Read_DT_Hits()
            cdc_driftD->push_back(mchit_DriftD/unit::mm); // take DOCA as driftD for this moment; Will be updated by Read_DT_Hits()
            cdc_nPairs->push_back(1); // Default value is 1; FIXME: need to take a value from Ion hit
            cdc_DOCA->push_back(mchit_DOCA/unit::mm);
            cdc_t->push_back(mchit_time/unit::ns);
            cdc_first_t->push_back(mchit_first_time/unit::ns);
            cdc_last_t->push_back(mchit_last_time/unit::ns);
            cdc_x->push_back(posLocal.X()/unit::mm);
            cdc_y->push_back(posLocal.Y()/unit::mm);
            cdc_z->push_back(posLocal.Z()/unit::mm);
            cdc_first_x->push_back(firstPosLocal.X()/unit::mm);
            cdc_first_y->push_back(firstPosLocal.Y()/unit::mm);
            cdc_first_z->push_back(firstPosLocal.Z()/unit::mm);
            cdc_last_x->push_back(lastPosLocal.X()/unit::mm);
            cdc_last_y->push_back(lastPosLocal.Y()/unit::mm);
            cdc_last_z->push_back(lastPosLocal.Z()/unit::mm);

            cdc_px->push_back(momLocal.Px()/unit::MeV);
            cdc_py->push_back(momLocal.Py()/unit::MeV);
            cdc_pz->push_back(momLocal.Pz()/unit::MeV);
            cdc_first_px->push_back(firstMomLocal.Px()/unit::MeV);
            cdc_first_py->push_back(firstMomLocal.Py()/unit::MeV);
            cdc_first_pz->push_back(firstMomLocal.Pz()/unit::MeV);
            cdc_last_px->push_back(lastMomLocal.Px()/unit::MeV);
            cdc_last_py->push_back(lastMomLocal.Py()/unit::MeV);
            cdc_last_pz->push_back(lastMomLocal.Pz()/unit::MeV);
            cdc_wx->push_back(mchit_wx/unit::mm);
            cdc_wy->push_back(mchit_wy/unit::mm);
            cdc_wz->push_back(mchit_wz/unit::mm);
        }
        int currentLayersEnter = 0;
        for (unsigned int ilayer = 0; ilayer<layerHitEnter.size(); ilayer++){
            if (layerHitEnter[ilayer])
                currentLayersEnter++;
            else
                currentLayersEnter=0;
            if (CLEntrance<currentLayersEnter) CLEntrance = currentLayersEnter;
        }
        int currentLayersExit = 0;
        for (unsigned int ilayer = 0; ilayer<layerHitExit.size(); ilayer++){
            if (layerHitExit[ilayer])
                currentLayersExit++;
            else
                currentLayersExit=0;
            if (CLExit<currentLayersExit) CLExit = currentLayersExit;
        }
    }
    else{
        COMETInfo("There is no g4Hits/CDC in this event!");
        return;
    }
}

/// Read Monte Carlo vertices
void ICyDetG4Analyzer::Read_MC_mom(COMET::ICOMETEvent& event){
    COMET::IHandle<COMET::IG4PrimaryVertexContainer> fPrimary_vertex_container =
        event.Get<COMET::IG4PrimaryVertexContainer>("truth/G4PrimVertex00");
    if(fPrimary_vertex_container){
        for(COMET::IG4PrimaryVertexContainer::iterator i_vertex=fPrimary_vertex_container->begin();
                i_vertex!=fPrimary_vertex_container->end(); ++i_vertex)
        {   
            const TLorentzVector position = i_vertex->GetPosition();
            ini_x = position.X();
            ini_y = position.Y();
            ini_z = position.Z();
            ini_t = position.T();
            const COMET::TG4PrimaryParticleContainer& primary_particles=i_vertex->GetPrimaryParticles();
            for(COMET::TG4PrimaryParticleContainer::const_iterator i_particle=primary_particles.begin();
                    i_particle!=primary_particles.end(); ++i_particle) 
            {   
                // FIXME: for now only take the first one
//                if(fabs((*i_particle).GetPDGCode()) != 11) continue;
                const TLorentzVector momentum = i_particle->GetMomentum();
                ini_px = momentum.Px()/unit::MeV;
                ini_py = momentum.Py()/unit::MeV;
                ini_pz = momentum.Pz()/unit::MeV;
                break; // only the first one...
            }
            break; // only the first one...
        }
    } 
    return;
}

/// Read Monte Carlo trigger
void ICyDetG4Analyzer::Read_MC_trigger(COMET::ICOMETEvent& event){
    COMET::IHandle<COMET::ISIMG4Header> mcHeader = event.Get<COMET::ISIMG4Header>("truth/mcHeader");
    trig_type = mcHeader->GetTriggerInfo();
    if (trig_type>=14) trig = 1;
    else trig = 0;
    // FIXME: should get more trigger info in the future when the interface is ready
}

int ICyDetG4Analyzer::InitTree(){
    h_mc_tree = (TTree*) new TTree("mc","tree from mc truth");
    eventId = -1;
    Delete();

    h_mc_tree->Branch("eventId",&eventId,"eventId/I");

    if (!fDisableMCCDCHits) {
        h_mc_tree->Branch("nHitsFT",&nHitsFT);
        h_mc_tree->Branch("nHitsAT",&nHitsAT);
        h_mc_tree->Branch("maxLayerId",&maxLayerId);
        h_mc_tree->Branch("nTurns",&nTurns);
        h_mc_tree->Branch("CLEntrance",&CLEntrance);
        h_mc_tree->Branch("CLExit",&CLExit);
    }
    if (!fDisableMCTrigger) {
        h_mc_tree->Branch("trig",&trig);
        h_mc_tree->Branch("trig_type",&trig_type);
        h_mc_tree->Branch("trig_time",&trig_time);
        h_mc_tree->Branch("trig_segId",&trig_segId);
    }
    if (!fDisableMCCDCHits) {
        h_mc_tree->Branch("first_t",&first_t);
        h_mc_tree->Branch("first_px",&first_px);
        h_mc_tree->Branch("first_py",&first_py);
        h_mc_tree->Branch("first_pz",&first_pz);
        h_mc_tree->Branch("first_x",&first_x);
        h_mc_tree->Branch("first_y",&first_y);
        h_mc_tree->Branch("first_z",&first_z);
        h_mc_tree->Branch("first_theta",&first_theta);
    }
    if (!fDisableMCPrimary) {
        h_mc_tree->Branch("ini_px",&ini_px);
        h_mc_tree->Branch("ini_py",&ini_py);
        h_mc_tree->Branch("ini_pz",&ini_pz);
        h_mc_tree->Branch("ini_x",&ini_x);
        h_mc_tree->Branch("ini_y",&ini_y);
        h_mc_tree->Branch("ini_z",&ini_z);
        h_mc_tree->Branch("ini_t",&ini_t);
    }
    if (!fDisableMCCTHHits) {
        h_mc_tree->Branch("cth_nHits",&cth_nHits);
        h_mc_tree->Branch("cth_segId",&cth_segId);
        h_mc_tree->Branch("cth_hodId",&cth_hodId);
        h_mc_tree->Branch("cth_cryType",&cth_cryType);
        h_mc_tree->Branch("cth_pdgId",&cth_pdgId);
        h_mc_tree->Branch("cth_trackId",&cth_trackId);
        h_mc_tree->Branch("cth_stepL",&cth_stepL);
        h_mc_tree->Branch("cth_edep",&cth_edep);
        h_mc_tree->Branch("cth_beta",&cth_beta);
        h_mc_tree->Branch("cth_trig",&cth_trig);
        h_mc_tree->Branch("cth_t",&cth_t);
        h_mc_tree->Branch("cth_px",&cth_px);
        h_mc_tree->Branch("cth_py",&cth_py);
        h_mc_tree->Branch("cth_pz",&cth_pz);
        h_mc_tree->Branch("cth_x",&cth_x);
        h_mc_tree->Branch("cth_y",&cth_y);
        h_mc_tree->Branch("cth_z",&cth_z);
    
        h_mc_tree->Branch("cthmc_nHits",&cthmc_nHits);
        h_mc_tree->Branch("cthmc_nCon",&cthmc_nCon);
        h_mc_tree->Branch("cthmc_segId",&cthmc_segId);
        h_mc_tree->Branch("cthmc_hodId",&cthmc_hodId);
        h_mc_tree->Branch("cthmc_cryType",&cthmc_cryType);
        h_mc_tree->Branch("cthmc_pdgId",&cthmc_pdgId);
        h_mc_tree->Branch("cthmc_trackId",&cthmc_trackId);
        h_mc_tree->Branch("cthmc_chargeP",&cthmc_chargeP);
        h_mc_tree->Branch("cthmc_chargeS",&cthmc_chargeS);
        h_mc_tree->Branch("cthmc_stepL",&cthmc_stepL);
        h_mc_tree->Branch("cthmc_edep",&cthmc_edep);
        h_mc_tree->Branch("cthmc_beta",&cthmc_beta);
        h_mc_tree->Branch("cthmc_trig",&cthmc_trig);
        h_mc_tree->Branch("cthmc_tT",&cthmc_tT);
        h_mc_tree->Branch("cthmc_tP",&cthmc_tP);
        h_mc_tree->Branch("cthmc_px",&cthmc_px);
        h_mc_tree->Branch("cthmc_py",&cthmc_py);
        h_mc_tree->Branch("cthmc_pz",&cthmc_pz);
        h_mc_tree->Branch("cthmc_x",&cthmc_x);
        h_mc_tree->Branch("cthmc_y",&cthmc_y);
        h_mc_tree->Branch("cthmc_z",&cthmc_z);
    
        h_mc_tree->Branch("cthTrig_nHits",&cthTrig_nHits);
        h_mc_tree->Branch("cthTrig_segId",&cthTrig_segId);
        h_mc_tree->Branch("cthTrig_hodId",&cthTrig_hodId);
        h_mc_tree->Branch("cthTrig_cryType",&cthTrig_cryType);
        h_mc_tree->Branch("cthTrig_t",&cthTrig_t);
        h_mc_tree->Branch("cthTrig_ADC",&cthTrig_ADC);
    }
    if (!fDisableMCCDCHits) {
        h_mc_tree->Branch("cdc_nHits",&cdc_nHits);
        h_mc_tree->Branch("cdc_cellId",&cdc_cellId);
        h_mc_tree->Branch("cdc_layerId",&cdc_layerId);
        h_mc_tree->Branch("cdc_pdgId",&cdc_pdgId);
        h_mc_tree->Branch("cdc_trackId",&cdc_trackId);
        h_mc_tree->Branch("cdc_hittype",&cdc_hittype);
        h_mc_tree->Branch("cdc_edep",&cdc_edep);
        h_mc_tree->Branch("cdc_trig",&cdc_trig);
        h_mc_tree->Branch("cdc_length",&cdc_length);
        h_mc_tree->Branch("cdc_driftD",&cdc_driftD);
        h_mc_tree->Branch("cdc_pulseStartT",&cdc_pulseStartT);
        h_mc_tree->Branch("cdc_pulseStopT",&cdc_pulseStopT);
        h_mc_tree->Branch("cdc_IO",&cdc_IO);
        h_mc_tree->Branch("cdc_nPairs",&cdc_nPairs);
        h_mc_tree->Branch("cdc_DOCA",&cdc_DOCA);
        h_mc_tree->Branch("cdc_t",&cdc_t);
        h_mc_tree->Branch("cdc_px",&cdc_px);
        h_mc_tree->Branch("cdc_py",&cdc_py);
        h_mc_tree->Branch("cdc_pz",&cdc_pz);
        h_mc_tree->Branch("cdc_first_px",&cdc_first_px);
        h_mc_tree->Branch("cdc_first_py",&cdc_first_py);
        h_mc_tree->Branch("cdc_first_pz",&cdc_first_pz);
        h_mc_tree->Branch("cdc_first_x",&cdc_first_x);
        h_mc_tree->Branch("cdc_first_y",&cdc_first_y);
        h_mc_tree->Branch("cdc_first_z",&cdc_first_z);
        h_mc_tree->Branch("cdc_first_t",&cdc_first_t);
        h_mc_tree->Branch("cdc_last_px",&cdc_last_px);
        h_mc_tree->Branch("cdc_last_py",&cdc_last_py);
        h_mc_tree->Branch("cdc_last_pz",&cdc_last_pz);
        h_mc_tree->Branch("cdc_last_x",&cdc_last_x);
        h_mc_tree->Branch("cdc_last_y",&cdc_last_y);
        h_mc_tree->Branch("cdc_last_z",&cdc_last_z);
        h_mc_tree->Branch("cdc_last_t",&cdc_last_t);
        h_mc_tree->Branch("cdc_x",&cdc_x);
        h_mc_tree->Branch("cdc_y",&cdc_y);
        h_mc_tree->Branch("cdc_z",&cdc_z);
        h_mc_tree->Branch("cdc_wx",&cdc_wx);
        h_mc_tree->Branch("cdc_wy",&cdc_wy);
        h_mc_tree->Branch("cdc_wz",&cdc_wz);
    }

    cth_segId = new std::vector<int>;
    cth_hodId = new std::vector<int>;
    cth_cryType = new std::vector<int>;
    cth_pdgId = new std::vector<int>;
    cth_trackId = new std::vector<int>;
    cth_stepL = new std::vector<double>;
    cth_edep = new std::vector<double>;
    cth_beta = new std::vector<double>;
    cth_trig = new std::vector<int>;
    cth_t = new std::vector<double>;
    cth_px = new std::vector<double>;
    cth_py = new std::vector<double>;
    cth_pz = new std::vector<double>;
    cth_x = new std::vector<double>;
    cth_y = new std::vector<double>;
    cth_z = new std::vector<double>;

    cthmc_nCon = new std::vector<int>;
    cthmc_segId = new std::vector<int>;
    cthmc_hodId = new std::vector<int>;
    cthmc_cryType = new std::vector<int>;
    cthmc_pdgId = new std::vector<int>;
    cthmc_trackId = new std::vector<int>;
    cthmc_chargeP = new std::vector<double>;
    cthmc_chargeS = new std::vector<double>;
    cthmc_stepL = new std::vector<double>;
    cthmc_edep = new std::vector<double>;
    cthmc_beta = new std::vector<double>;
    cthmc_trig = new std::vector<int>;
    cthmc_tT = new std::vector<double>;
    cthmc_tP = new std::vector<double>;
    cthmc_px = new std::vector<double>;
    cthmc_py = new std::vector<double>;
    cthmc_pz = new std::vector<double>;
    cthmc_x = new std::vector<double>;
    cthmc_y = new std::vector<double>;
    cthmc_z = new std::vector<double>;

    cthTrig_segId = new std::vector<int>;
    cthTrig_hodId = new std::vector<int>;
    cthTrig_cryType = new std::vector<int>;
    cthTrig_t = new std::vector<int>;
    cthTrig_ADC = new std::vector<int>;

    cdc_cellId = new std::vector<int>;
    cdc_layerId = new std::vector<int>;
    cdc_pdgId = new std::vector<int>;
    cdc_trackId = new std::vector<int>;
    cdc_hittype = new std::vector<int>;
    cdc_edep = new std::vector<double>;
    cdc_trig = new std::vector<int>;
    cdc_length = new std::vector<double>;
    cdc_driftD = new std::vector<double>;
    cdc_pulseStartT = new std::vector<double>;
    cdc_pulseStopT = new std::vector<double>;
    cdc_IO = new std::vector<double>;
    cdc_nPairs = new std::vector<int>;
    cdc_DOCA = new std::vector<double>;
    cdc_t = new std::vector<double>;
    cdc_first_t = new std::vector<double>;
    cdc_last_t = new std::vector<double>;
    cdc_px = new std::vector<double>;
    cdc_py = new std::vector<double>;
    cdc_pz = new std::vector<double>;
    cdc_first_px = new std::vector<double>;
    cdc_first_py = new std::vector<double>;
    cdc_first_pz = new std::vector<double>;
    cdc_last_px = new std::vector<double>;
    cdc_last_py = new std::vector<double>;
    cdc_last_pz = new std::vector<double>;
    cdc_x = new std::vector<double>;
    cdc_y = new std::vector<double>;
    cdc_z = new std::vector<double>;
    cdc_first_x = new std::vector<double>;
    cdc_first_y = new std::vector<double>;
    cdc_first_z = new std::vector<double>;
    cdc_last_x = new std::vector<double>;
    cdc_last_y = new std::vector<double>;
    cdc_last_z = new std::vector<double>;
    cdc_wx = new std::vector<double>;
    cdc_wy = new std::vector<double>;
    cdc_wz = new std::vector<double>;
    return 1;
}

void ICyDetG4Analyzer::PrepareOutputs() {
    eventId = -1;

    nHitsFT = 0;
    nHitsAT = 0;
    maxLayerId = 0;
    nTurns = 0;
    CLEntrance = 0;
    CLExit = 0;

    trig = 0;
    trig_type = 0;
    trig_time = 0;
    trig_segId = 0;

    first_t = 0;
    first_px = 0;
    first_py = 0;
    first_pz = 0;
    first_x = 0;
    first_y = 0;
    first_z = 0;
    first_theta = 0;
    ini_px = 0;
    ini_py = 0;
    ini_pz = 0;
    ini_x = 0;
    ini_y = 0;
    ini_z = 0;
    ini_t = 0;

    cth_nHits = 0;
    cth_segId->clear();
    cth_hodId->clear();
    cth_cryType->clear();
    cth_pdgId->clear();
    cth_trackId->clear();
    cth_stepL->clear();
    cth_edep->clear();
    cth_beta->clear();
    cth_trig->clear();
    cth_t->clear();
    cth_px->clear();
    cth_py->clear();
    cth_pz->clear();
    cth_x->clear();
    cth_y->clear();
    cth_z->clear();

    cthmc_nHits = 0;
    cthmc_nCon->clear();
    cthmc_segId->clear();
    cthmc_hodId->clear();
    cthmc_cryType->clear();
    cthmc_pdgId->clear();
    cthmc_trackId->clear();
    cthmc_chargeP->clear();
    cthmc_chargeS->clear();
    cthmc_stepL->clear();
    cthmc_edep->clear();
    cthmc_beta->clear();
    cthmc_trig->clear();
    cthmc_tT->clear();
    cthmc_tP->clear();
    cthmc_px->clear();
    cthmc_py->clear();
    cthmc_pz->clear();
    cthmc_x->clear();
    cthmc_y->clear();
    cthmc_z->clear();

    cthTrig_nHits = 0;
    cthTrig_segId->clear();
    cthTrig_hodId->clear();
    cthTrig_cryType->clear();
    cthTrig_t->clear();
    cthTrig_ADC->clear();

    cdc_nHits = 0;
    cdc_cellId->clear();
    cdc_layerId->clear();
    cdc_pdgId->clear();
    cdc_trackId->clear();
    cdc_hittype->clear();
    cdc_edep->clear();
    cdc_trig->clear();
    cdc_length->clear();
    cdc_driftD->clear();
    cdc_pulseStartT->clear();
    cdc_pulseStopT->clear();
    cdc_IO->clear();
    cdc_nPairs->clear();
    cdc_DOCA->clear();
    cdc_t->clear();
    cdc_first_t->clear();
    cdc_last_t->clear();
    cdc_px->clear();
    cdc_py->clear();
    cdc_pz->clear();
    cdc_first_px->clear();
    cdc_first_py->clear();
    cdc_first_pz->clear();
    cdc_last_px->clear();
    cdc_last_py->clear();
    cdc_last_pz->clear();
    cdc_x->clear();
    cdc_y->clear();
    cdc_z->clear();
    cdc_first_x->clear();
    cdc_first_y->clear();
    cdc_first_z->clear();
    cdc_last_x->clear();
    cdc_last_y->clear();
    cdc_last_z->clear();
    cdc_wx->clear();
    cdc_wy->clear();
    cdc_wz->clear();
}

int ICyDetG4Analyzer::BeginFile(COMET::IVInputFile *const input) {
    TGeoManager::Import(input->GetFilename());
    return 1;
}

bool ICyDetG4Analyzer::SetOptions(std::map<std::string,std::string> options){
    std::map<std::string,std::string>::iterator it;
    for(it = options.begin(); it != options.end(); it++){
        std::string opt = (*it).first;
        std::string value = (*it).second;
        if ( opt == "Hello" ) {
            std::cout<<"Hello, user! Good luck!"<<std::endl;
        }
        else if(opt=="parameters") {
            COMET::IOARuntimeParameters::Get().ReadParamOverrideFile(value);
        }
        else if(opt=="disable-primary") {
            fDisableMCPrimary = true;
        }
        else if(opt=="disable-trigger") {
            fDisableMCTrigger = true;
        }
        else if(opt=="disable-cth") {
            fDisableMCCTHHits = true;
        }
        else if(opt=="disable-cdc") {
            fDisableMCCDCHits = true;
        }
        else {
            std::cout<<"Warning: This option is not defined here: "<<opt<<std::endl;
        }
    }
    return true;
}

int ICyDetG4Analyzer::Init() {
    fProtonMass  = this->GetMass("proton");
    COMETInfo( "     Proton  mass = "<<fProtonMass<<" MeV");
    fNeutronMass = this->GetMass("neutron");
    COMETInfo( "     Neutron mass = "<<fNeutronMass<<" MeV");
    int a = InitTree();  
    return a;
}

int ICyDetG4Analyzer::Finalize() {
    if(h_mc_tree!=NULL) h_mc_tree->Write();
    //Delete();/// FIXME: If some of vectors won't be saved in the output file,
             ///  they won't be automatically deleted by ROOT.
             ///  Should this be done after writing the ROOT file to avoid memory leak,
             ///  or can be done here?
    return 0;
}

void ICyDetG4Analyzer::Delete() {
    if(cth_segId) {delete cth_segId;} cth_segId = 0;
    if(cth_hodId) {delete cth_hodId;} cth_hodId = 0;
    if(cth_cryType) {delete cth_cryType;} cth_cryType = 0;
    if(cth_pdgId) {delete cth_pdgId;} cth_pdgId = 0;
    if(cth_trackId) {delete cth_trackId;} cth_trackId = 0;
    if(cth_stepL) {delete cth_stepL;} cth_stepL = 0;
    if(cth_edep) {delete cth_edep;} cth_edep = 0;
    if(cth_beta) {delete cth_beta;} cth_beta = 0;
    if(cth_trig) {delete cth_trig;} cth_trig = 0;
    if(cth_t) {delete cth_t;} cth_t = 0;
    if(cth_px) {delete cth_px;} cth_px = 0;
    if(cth_py) {delete cth_py;} cth_py = 0;
    if(cth_pz) {delete cth_pz;} cth_pz = 0;
    if(cth_x) {delete cth_x;} cth_x = 0;
    if(cth_y) {delete cth_y;} cth_y = 0;
    if(cth_z) {delete cth_z;} cth_z = 0;
    if(cthmc_nCon) {delete cthmc_nCon;} cthmc_nCon = 0;
    if(cthmc_segId) {delete cthmc_segId;} cthmc_segId = 0;
    if(cthmc_hodId) {delete cthmc_hodId;} cthmc_hodId = 0;
    if(cthmc_cryType) {delete cthmc_cryType;} cthmc_cryType = 0;
    if(cthmc_pdgId) {delete cthmc_pdgId;} cthmc_pdgId = 0;
    if(cthmc_trackId) {delete cthmc_trackId;} cthmc_trackId = 0;
    if(cthmc_stepL) {delete cthmc_stepL;} cthmc_stepL = 0;
    if(cthmc_chargeP) {delete cthmc_chargeP;} cthmc_chargeP = 0;
    if(cthmc_chargeS) {delete cthmc_chargeS;} cthmc_chargeS = 0;
    if(cthmc_edep) {delete cthmc_edep;} cthmc_edep = 0;
    if(cthmc_beta) {delete cthmc_beta;} cthmc_beta = 0;
    if(cthmc_trig) {delete cthmc_trig;} cthmc_trig = 0;
    if(cthmc_tT) {delete cthmc_tT;} cthmc_tT = 0;
    if(cthmc_tP) {delete cthmc_tP;} cthmc_tP = 0;
    if(cthmc_px) {delete cthmc_px;} cthmc_px = 0;
    if(cthmc_py) {delete cthmc_py;} cthmc_py = 0;
    if(cthmc_pz) {delete cthmc_pz;} cthmc_pz = 0;
    if(cthmc_x) {delete cthmc_x;} cthmc_x = 0;
    if(cthmc_y) {delete cthmc_y;} cthmc_y = 0;
    if(cthmc_z) {delete cthmc_z;} cthmc_z = 0;
    if(cthTrig_segId) {delete cthTrig_segId;} cthTrig_segId = 0;
    if(cthTrig_hodId) {delete cthTrig_hodId;} cthTrig_hodId = 0;
    if(cthTrig_cryType) {delete cthTrig_cryType;} cthTrig_cryType = 0;
    if(cthTrig_t) {delete cthTrig_t;} cthTrig_t = 0;
    if(cthTrig_ADC) {delete cthTrig_ADC;} cthTrig_ADC = 0;
    if(cdc_cellId) {delete cdc_cellId;} cdc_cellId = 0;
    if(cdc_layerId) {delete cdc_layerId;} cdc_layerId = 0;
    if(cdc_pdgId) {delete cdc_pdgId;} cdc_pdgId = 0;
    if(cdc_trackId) {delete cdc_trackId;} cdc_trackId = 0;
    if(cdc_hittype) {delete cdc_hittype;} cdc_hittype = 0;
    if(cdc_edep) {delete cdc_edep;} cdc_edep = 0;
    if(cdc_trig) {delete cdc_trig;} cdc_trig = 0;
    if(cdc_length) {delete cdc_length;} cdc_length = 0;
    if(cdc_driftD) {delete cdc_driftD;} cdc_driftD = 0;
    if(cdc_pulseStartT) {delete cdc_pulseStartT;} cdc_pulseStartT = 0;
    if(cdc_pulseStopT) {delete cdc_pulseStopT;} cdc_pulseStopT = 0;
    if(cdc_IO) {delete cdc_IO;} cdc_IO = 0;
    if(cdc_nPairs) {delete cdc_nPairs;} cdc_nPairs = 0;
    if(cdc_DOCA) {delete cdc_DOCA;} cdc_DOCA = 0;
    if(cdc_t) {delete cdc_t;} cdc_t = 0;
    if(cdc_px) {delete cdc_px;} cdc_px = 0;
    if(cdc_py) {delete cdc_py;} cdc_py = 0;
    if(cdc_pz) {delete cdc_pz;} cdc_pz = 0;
    if(cdc_first_px) {delete cdc_first_px;} cdc_first_px = 0;
    if(cdc_first_py) {delete cdc_first_py;} cdc_first_py = 0;
    if(cdc_first_pz) {delete cdc_first_pz;} cdc_first_pz = 0;
    if(cdc_first_x) {delete cdc_first_x;} cdc_first_x = 0;
    if(cdc_first_y) {delete cdc_first_y;} cdc_first_y = 0;
    if(cdc_first_z) {delete cdc_first_z;} cdc_first_z = 0;
    if(cdc_first_t) {delete cdc_first_t;} cdc_first_t = 0;
    if(cdc_last_px) {delete cdc_last_px;} cdc_last_px = 0;
    if(cdc_last_py) {delete cdc_last_py;} cdc_last_py = 0;
    if(cdc_last_pz) {delete cdc_last_pz;} cdc_last_pz = 0;
    if(cdc_last_x) {delete cdc_last_x;} cdc_last_x = 0;
    if(cdc_last_y) {delete cdc_last_y;} cdc_last_y = 0;
    if(cdc_last_z) {delete cdc_last_z;} cdc_last_z = 0;
    if(cdc_last_t) {delete cdc_last_t;} cdc_last_t = 0;
    if(cdc_x) {delete cdc_x;} cdc_x = 0;
    if(cdc_y) {delete cdc_y;} cdc_y = 0;
    if(cdc_z) {delete cdc_z;} cdc_z = 0;
    if(cdc_wx) {delete cdc_wx;} cdc_wx = 0;
    if(cdc_wy) {delete cdc_wy;} cdc_wy = 0;
    if(cdc_wz) {delete cdc_wz;} cdc_wz = 0;
}

double ICyDetG4Analyzer::GetMass(const int& pdgCode,double& charge) {
    //// PDGCode larger than 1000000000 means this is a nuclear
    if (pdgCode<1000000000) {
        TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(pdgCode);
        charge = part->Charge()/3.0;
        return part->Mass()*unit::GeV; /// Convert GeV to MeV
    } else {
        int AAA = (pdgCode% 10000) / 10;
        int ZZZ = (pdgCode/10000)% 1000;
        charge = 0.0;
        return (fProtonMass*(double)ZZZ+fNeutronMass*(double)(AAA-ZZZ));
    }
}

double ICyDetG4Analyzer::GetMass(const std::string& partName) {
    TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(partName.c_str());
    return part->Mass()*unit::GeV; /// Convert GeV to MeV
}
