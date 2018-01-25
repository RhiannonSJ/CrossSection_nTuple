////////////////////////////////////////////////////////////////////////
// Class:       AnalysisNTuple
// Plugin Type: analyzer (art v2_08_04)
// File:        AnalysisNTuple_module.cc
//
// Generated at Tue Jan 23 18:32:51 2018 by Rhiannon Jones using 
// cetskelgen from cetlib version v3_01_01.
//
////////////////////////////////////////////////////////////////////////
//
// This module constructs a root file which contains 
//  MCParticle, 
//  Reco Track, 
//  Reco Shower,
//  True/Reco event information 
//  for use in topology-based analyses
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include <sstream>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "TROOT.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"


namespace pndr {
  class AnalysisNTuple;
}


class pndr::AnalysisNTuple : public art::EDAnalyzer {
public:
  explicit AnalysisNTuple(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalysisNTuple(AnalysisNTuple const &) = delete;
  AnalysisNTuple(AnalysisNTuple &&) = delete;
  AnalysisNTuple & operator = (AnalysisNTuple const &) = delete;
  AnalysisNTuple & operator = (AnalysisNTuple &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  // Declare member data here.
  // Member data to check there is a numuCC interaction within the fiducial
  // volume
  std::map< std::vector< int >, int > m_selection;
  float m_detectorHalfLengthX;
  float m_detectorHalfLengthY;
  float m_detectorHalfLengthZ;
  float m_coordinateOffsetX;
  float m_coordinateOffsetY;
  float m_coordinateOffsetZ;
  float m_selectedBorderX;
  float m_selectedBorderY;
  float m_selectedBorderZ;

  // Counters
  int all_events, fiducial_contained_events;
  int all_tracks, fiducial_contained_tracks;

  // ROOT
  TTree *event_tree, *mcparticle_tree, *recotrack_tree, *recoshower_tree;

  // Variables associated with event_tree
  // Universal
  int event_id;

  // Truth
  int t_nu_pdgcode, t_interaction;
  int t_particles;
  int t_protons, t_muons, t_charged_pions, t_kaons, t_neutral_pions;
  double t_inv_mass, t_nu_lepton_angle, t_vertex_energy;
  double t_vertex[3];
  
  // Reco
  int r_particles, r_tracks, r_showers;
  double r_vertex[3];

  // Variables associated with mcparticle_tree
  // Variables associated with recotrack_tree
  // Variables associated with recoshower_tree
  
};


pndr::AnalysisNTuple::AnalysisNTuple(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{

  this->reconfigure(p);

}

void pndr::AnalysisNTuple::analyze(art::Event const & e)
{
  
  all_events++;

  //`typedef::std::map< std::vector< int >, int > topology_map;
  
  // Implementation of required member function here.
  bool contained = true;

  // Get the MCTruth information 
  art::Handle< std::vector< simb::MCTruth > > mct_handle;
  e.getByLabel("generator", mct_handle );
  int mct_size = mct_handle->size();
 
  // Get tracks and check they remain in the fiducial volume 
  art::Handle< std::vector< recob::Track > > trk_handle;
  e.getByLabel("pmalgtrackmaker", trk_handle );
  int trk_size = trk_handle->size();
 
  all_tracks += trk_size;

  /*
  art::Handle< std::vector< recob::Shower > > shw_handle;
  e.getByLabel("pandoraNu", shw_handle );
  int shw_size = shw_handle->size();
  */

  if(mct_handle.isValid() && mct_size) { 
  
    // Loop over the truth info
    for(auto const& mct : (*mct_handle)) {
 
      // Check the neutrino came from the beam
      if(mct.Origin() != simb::kBeamNeutrino) continue;
 
      // Check the neutrino interaction vertex is within the fiducial volume
      float nu_vertex_x = mct.GetNeutrino().Lepton().Vx();
      float nu_vertex_y = mct.GetNeutrino().Lepton().Vy();
      float nu_vertex_z = mct.GetNeutrino().Lepton().Vz();
 
      if (    (nu_vertex_x > (m_detectorHalfLengthX - m_coordinateOffsetX - m_selectedBorderX)) 
           || (nu_vertex_x < (-m_coordinateOffsetX + m_selectedBorderX)) 
           || (nu_vertex_y > (m_detectorHalfLengthY - m_coordinateOffsetY - m_selectedBorderY)) 
           || (nu_vertex_y < (-m_coordinateOffsetY + m_selectedBorderY)) 
           || (nu_vertex_z > (m_detectorHalfLengthZ - m_coordinateOffsetZ - m_selectedBorderZ)) 
           || (nu_vertex_z < (-m_coordinateOffsetZ + m_selectedBorderZ))) contained = false;
    
    } 
  }
  // Proceed with the nTuple-filling if we are in the fiducial volume 
  // and all tracks are contained
  if(contained){
  
    fiducial_contained_events++;
    
    // Get the MCTruth information 
    art::Handle< std::vector< simb::MCTruth > > mct_handle;
    e.getByLabel("generator", mct_handle );
    int mct_size = mct_handle->size();
  
    art::Handle< std::vector< recob::Track > > trk_handle;
    e.getByLabel("pmalgtrackmaker", trk_handle );
    int trk_size = trk_handle->size();
   
    art::Handle< std::vector< recob::Shower > > shw_handle;
    e.getByLabel("pandoraNu", shw_handle );
    int shw_size = shw_handle->size();

    art::Handle< std::vector< recob::PFParticle > > pfp_handle;
    e.getByLabel("pandoraNu", pfp_handle );
    int pfp_size = pfp_handle->size();
    
// ------------------------------------------------------------------------------
//                           EVENT-TREE INFORMATION
// ------------------------------------------------------------------------------
    if(mct_handle.isValid() && mct_size && 
       shw_handle.isValid() && 
       pfp_handle.isValid() && pfp_size && 
       trk_handle.isValid() ) {
  
// ------------------------------------------------------------------------------
//                     EVENT-TREE RECONSTRUCTED INFORMATION
// ------------------------------------------------------------------------------
      
      // Initialise counter for the number of reconstructed primary particles
      int n_primaries           = 0;
      int n_primary_tracks     = 0;
      int n_primary_showers    = 0;
      unsigned int neutrino_id = 0;
     
      bool neutrino_found = false;
      bool muon_neutrino = false;

      // Loop over PFParticles and find how many are primary and whether 
      // only 1 neutrino was found
      for(int i = 0; i < pfp_size; ++i) {
    
        art::Ptr< recob::PFParticle > pfp( pfp_handle, i );

        // Only look for a single muon neutrino event
        if(!neutrino_found && pfp->IsPrimary()){
          neutrino_found = true;
          if(pfp->PdgCode() == 14) {
            muon_neutrino = true;
            neutrino_id = pfp->Self();
          }
        }
        else if(pfp->IsPrimary()) return;
      
      }

      if(!muon_neutrino) return;

      // Get vertex association
      art::FindMany< recob::Vertex  > fvtx( pfp_handle, e, "pandoraNu" );
      std::vector<const recob::Vertex*> vtx_assn = fvtx.at(neutrino_id);

      if(vtx_assn.size()  > 1) return;
      if(vtx_assn.size() == 0) return;

      // Set array to be current vertex position
      vtx_assn[0]->XYZ(r_vertex);
      
      // Find the number of reconstructed primary final state particles 
      for(int i = 0; i < pfp_size; ++i) {
    
        art::Ptr< recob::PFParticle > pfp( pfp_handle, i );

        if(neutrino_id == pfp->Parent()) n_primaries++;
      
      }

      if(trk_handle.isValid() && trk_size) {
      
        // Loop over PMA tracks and find any within 2 cm of the primary vertex
        // count them
        for(int i = 0; i < trk_size; ++i) {
        
          art::Ptr< recob::Track > trk( trk_handle, i );
         
          float track_vtx_x = trk->Vertex()[0];
          float track_vtx_y = trk->Vertex()[1];
          float track_vtx_z = trk->Vertex()[2];
          float track_end_x = trk->End()[0];
          float track_end_y = trk->End()[1];
          float track_end_z = trk->End()[2];
          
          // Find the distance of the current track's vertex and end point from 
          // the primary vertex location
          // If the end point is closer, the track was reconstructed the wrong way around
          double s_vtx = sqrt(pow(track_vtx_x - r_vertex[0],2) + pow(track_vtx_y - r_vertex[1],2) + pow(track_vtx_z - r_vertex[2],2));
          double s_end = sqrt(pow(track_end_x - r_vertex[0],2) + pow(track_end_y - r_vertex[1],2) + pow(track_end_z - r_vertex[2],2));

          // If the end is closer than the vertex, flip the track
          bool should_flip = s_end < s_vtx;

          // If the track is the right way around and the vertex is within 2cm of 
          // the primary vertex location, add a counter to the number of primary tracks
          // Or if the track is the wrong way around and the end is within 10cm of the
          // primary vertex location, add a counter to the number of primary tracks
          if((!should_flip && s_vtx > 10) || (should_flip  && s_end > 10)) continue;
            
          // Check that the primary track's start and end position is within the detector volume
          if (    (track_vtx_x > (m_detectorHalfLengthX - m_coordinateOffsetX - m_selectedBorderX)) 
               || (track_vtx_x < (-m_coordinateOffsetX + m_selectedBorderX)) 
               || (track_vtx_y > (m_detectorHalfLengthY - m_coordinateOffsetY - m_selectedBorderY)) 
               || (track_vtx_y < (-m_coordinateOffsetY + m_selectedBorderY)) 
               || (track_vtx_z > (m_detectorHalfLengthZ - m_coordinateOffsetZ - m_selectedBorderZ)) 
               || (track_vtx_z < (-m_coordinateOffsetZ + m_selectedBorderZ))
               || (track_end_x > (m_detectorHalfLengthX - m_coordinateOffsetX - m_selectedBorderX)) 
               || (track_end_x < (-m_coordinateOffsetX + m_selectedBorderX)) 
               || (track_end_y > (m_detectorHalfLengthY - m_coordinateOffsetY - m_selectedBorderY)) 
               || (track_end_y < (-m_coordinateOffsetY + m_selectedBorderY)) 
               || (track_end_z > (m_detectorHalfLengthZ - m_coordinateOffsetZ - m_selectedBorderZ)) 
               || (track_end_z < (-m_coordinateOffsetZ + m_selectedBorderZ))) continue;
            
          // Try and eliminate shower fragments
          if(trk->Length() > 0.5) n_primary_tracks++;
        } 
      }
      if(shw_handle.isValid() && shw_size) {
      
        // Loop over PMA tracks and find any within 2 cm of the primary vertex
        // count them
        for(int i = 0; i < shw_size; ++i) {
        
          art::Ptr< recob::Shower > shw( shw_handle, i );
         
          // Distance s of the shower start from the reconstructed 
          // neutrino vertex
          double s_vtx = sqrt(pow(shw->ShowerStart()[0] - r_vertex[0],2) + pow(shw->ShowerStart()[1] - r_vertex[1],2) + pow(shw->ShowerStart()[2] - r_vertex[2],2));
         
          // Cut of 40 cm for showers to accommodate photon conversion after
          // neutral pion decay
          if(s_vtx < 40) n_primary_showers++;
        }
      }

      // Number of reconstructed primary particles, primary tracks and 
      // all reconstructed showers
      // Location of the primary reconstructed vertex
      r_tracks    = n_primary_tracks;
      r_showers   = n_primary_showers; 
      r_particles = n_primaries;
      
// ------------------------------------------------------------------------------
//                     EVENT-TREE TRUTH INFORMATION
// ------------------------------------------------------------------------------

      art::FindMany< simb::MCParticle  > fmcp( mct_handle, e, "largeant" );
      
      // Loop over the truth handle and get everything we can
      for(int i = 0; i < mct_size; ++i) {
    
        std::cout << "----------------------------------------------------------------" << std::endl;
        art::Ptr< simb::MCTruth > mct( mct_handle, i );

        std::cout << " Interaction type : " << mct->GetNeutrino().InteractionType() << std::endl;

        // Number of true final state particles
        t_particles = mct->NParticles() - 1;

        for( int j = 0; j < mct->NParticles(); ++j ) {
        
          simb::MCParticle part = mct->GetParticle(j);
          int motherPdg(mct->GetParticle(part.Mother()).PdgCode());
          
          std::cout << "----------------------------------------------------------------" << std::endl;
          std::cout << " ID : " << part.TrackId() << ", PdgCode : " << part.PdgCode() << ", Mother PdgCode : " << motherPdg << ", Mother ID : " << part.Mother();
          std::cout << ", Process : " << part.Process() << std::endl;
         
          /*
          if(part.Mother() != -1){
            int grandmotherPdg(mct->GetParticle(mct->GetParticle(part.Mother()).Mother()).PdgCode());

            if((motherPdg == 14 || grandmotherPdg == 1000180400) && part.PdgCode() != 1000180390){
            
              std::cout << "----------------------------------------------------------------" << std::endl;
              std::cout << " ID : " << part.TrackId() << ", PdgCode : " << part.PdgCode() << ", Mother PdgCode : " << motherPdg << ", Mother ID : " << part.Mother();
              std::cout << ", Process : " << part.Process() << std::endl;
            
            }
          }*/
        }
        
        std::cout << "====================================================================" << std::endl;

        std::vector<const simb::MCParticle*> mcp_assn = fmcp.at(i);
        
        for(const simb::MCParticle* part : mcp_assn) {

          int motherPdg(mct->GetParticle(part->Mother()).PdgCode());
          
          std::cout << "----------------------------------------------------------------" << std::endl;
          std::cout << " ID : " << part->TrackId() << ", PdgCode : " << part->PdgCode() << ", Mother PdgCode : " << motherPdg << ", Mother ID : " << part->Mother();
          std::cout << ", Process : " << part->Process() << std::endl;
         
          /*
          if(part.Mother() != -1){
            int grandmotherPdg(mct->GetParticle(mct->GetParticle(part.Mother()).Mother()).PdgCode());

            if((motherPdg == 14 || grandmotherPdg == 1000180400) && part.PdgCode() != 1000180390){
            
              std::cout << "----------------------------------------------------------------" << std::endl;
              std::cout << " ID : " << part.TrackId() << ", PdgCode : " << part.PdgCode() << ", Mother PdgCode : " << motherPdg << ", Mother ID : " << part.Mother();
              std::cout << ", Process : " << part.Process() << std::endl;
            
            }
          }*/
        
        }
        /*// Loop over particles starting from the second since the first is 
        // the neutrino
        for(unsigned int j = 1; < mct->NParticles(); ++j){
      
          simb::MCParticle part = mct->GetParticle(j);

        
        }*/
      }
      event_id += 1;
  
      // Fill the event tree once everything has been set
      event_tree->Fill();

    }
  }
}

void pndr::AnalysisNTuple::beginJob()
{
  // Implementation of optional member function here.
  // Initialise the counters
  all_events                = 0;
  all_tracks                = 0;
  fiducial_contained_events = 0;

  event_id = 0;

  // Initiate tree
  event_tree      = new TTree("event_tree",      "Event tree: True and reconstructed SBND event information");
  mcparticle_tree = new TTree("particle_tree",   "MCParticle tree: True SBND initial and final state topology information");
  recotrack_tree  = new TTree("recotrack_tree",  "Track tree: Reconstructed final state track information");
  recoshower_tree = new TTree("recoshower_tree", "Shower tree: Reconstructed final state shower information");

  // Event tree branches
  event_tree->Branch("event_id",          &event_id,          "event_id/I");
  event_tree->Branch("t_nu_pdgcode",      &t_nu_pdgcode,      "t_nu_pdgcode/I");
  event_tree->Branch("t_interaction",     &t_interaction,     "t_interaction/I");
  event_tree->Branch("t_particles",       &t_particles,       "t_particles/I");
  event_tree->Branch("t_protons",         &t_protons,         "t_protons/I");
  event_tree->Branch("t_muons",           &t_muons,           "t_muons/I");
  event_tree->Branch("t_charged_pions",   &t_charged_pions,   "t_charged_pions/I");
  event_tree->Branch("t_neutral_pions",   &t_neutral_pions ,  "t_neutral_pions/I");
  event_tree->Branch("t_kaons",           &t_kaons,           "t_kaons/I");
  event_tree->Branch("t_inv_mass",        &t_inv_mass,        "t_inv_mass/D");
  event_tree->Branch("t_nu_lepton_angle", &t_nu_lepton_angle, "t_nu_lepton_angle/D");
  event_tree->Branch("t_vertex_energy",   &t_vertex_energy,   "t_vertex_energy/D");
  event_tree->Branch("t_vertex",          &t_vertex,          "t_vertex[3]/D");
  event_tree->Branch("r_particles",       &r_particles,       "r_particles/I");
  event_tree->Branch("r_tracks",          &r_tracks,          "r_tracks/I");
  event_tree->Branch("r_showers",         &r_showers,         "r_showers/I");
  event_tree->Branch("r_vertex",          &r_vertex,          "r_vertex[3]/D");

  // MCParticle tree branches

  // Reco Track tree branches
  
  // Reco Shower tree branches

  
  event_tree->SetDirectory(0);
  mcparticle_tree->SetDirectory(0);
  recotrack_tree->SetDirectory(0);
  recoshower_tree->SetDirectory(0);

}

void pndr::AnalysisNTuple::endJob()
{
  // Implementation of optional member function here.

  std::cout << "=================================================================================" << std::endl;
  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  std::cout << std::endl;

  std::cout << " Total number of events                            : ";
  std::cout << all_events;
  std::cout << std::endl;

  std::cout << " Fraction of fiducial events with contained tracks : ";
  std::cout << fiducial_contained_events / double(all_events);
  std::cout << std::endl;

  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  std::cout << "=================================================================================" << std::endl;

  // Print the tree, write the file, close
  TFile file("/sbnd/app/users/rsjones/LArSoft_v06_63_00/LArSoft-v06_63_00/srcs/recoparameters/recoparameters/output_files/tree_test.root", "RECREATE");
  event_tree->Write();
  mcparticle_tree->Write();
  file.Write();
  file.Close();

  delete event_tree;
  delete mcparticle_tree;
  delete recotrack_tree;
  delete recoshower_tree;

}

void pndr::AnalysisNTuple::reconfigure(fhicl::ParameterSet const & p)
{
   std::vector< int > blankVect;
   std::vector< std::vector< int > > input;
 
   std::vector< int > selection1 = p.get< std::vector< int > >("Selection1",        blankVect);
   if ( selection1.size() != 0 ) input.push_back(selection1);
 
   std::vector< int > selection2 = p.get< std::vector< int > >("Selection2",        blankVect);
   if ( selection2.size() != 0 ) input.push_back(selection2);
 
   for ( auto & inputVect : input ) {
     if ( inputVect.size() < 2 ) {
       std::cerr << " Error: Selection vector must have at least 2 elements " <<    std::endl;
       std::cerr << "        First element:     Number of particles of PDG code(s)  specified " << std::endl;
       std::cerr << "        Remaining element: PDG codes to filter on " << std::   endl;
       exit(1);
     }
 
     int count = inputVect[0];
     inputVect.erase( inputVect.begin() );
 
     m_selection.insert( std::make_pair( inputVect, count ) );
   }
 
  // Implementation of optional member function here. 
   m_detectorHalfLengthX = p.get<float>("DetectorHalfLengthX");
   m_detectorHalfLengthY = p.get<float>("DetectorHalfLengthY");
   m_detectorHalfLengthZ = p.get<float>("DetectorHalfLengthZ");
   m_coordinateOffsetX   = p.get<float>("CoordinateOffsetX");
   m_coordinateOffsetY   = p.get<float>("CoordinateOffsetY");
   m_coordinateOffsetZ   = p.get<float>("CoordinateOffsetZ");
   m_selectedBorderX     = p.get<float>("SelectedBorderX");
   m_selectedBorderY     = p.get<float>("SelectedBorderY");
   m_selectedBorderZ     = p.get<float>("SelectedBorderZ");
}

DEFINE_ART_MODULE(pndr::AnalysisNTuple)
