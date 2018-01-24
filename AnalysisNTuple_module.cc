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
  // Member data to check the neutrino interaction occurs within the fiducial
  // volume
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
  TTree *event_tree, *particle_tree;
  float test_e, test_p;
  
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

  if(mct_handle.isValid() && mct_size && 
     trk_handle.isValid() && trk_size ) { 
  
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
   
      else{
        // Check that any reconstructed tracks lie within the fiducial volume
        for(auto const& trk : (*trk_handle)) {

          float track_vtx_x = trk.Vertex()[0];
          float track_vtx_y = trk.Vertex()[1];
          float track_vtx_z = trk.Vertex()[2];
          float track_end_x = trk.End()[0];
          float track_end_y = trk.End()[1];
          float track_end_z = trk.End()[2];
          
          // Check that the track start and end position is within the detector volume
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
               || (track_end_z < (-m_coordinateOffsetZ + m_selectedBorderZ))) contained = false;
          else{
            fiducial_contained_tracks++;
          }
        }
      }
    } 
  }
  // Proceed with the nTuple-filling if we are in the fiducial volume 
  // and all tracks are contained
  if(contained){
  
    fiducial_contained_events++;
    
    //static EVENT event;
    //static PARTICLE particle;

    // TTree example
    test_e = 1.5;

    event_tree->Fill();

    for(int i = 0; i < 100; ++i){
    
      test_p = 3.4*i;
      particle_tree->Fill();

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
  fiducial_contained_events = 0;

  // Initiate tree
  event_tree    = new TTree("event_tree",    "Event tree holding information about SBND final state topologies");
  particle_tree = new TTree("particle_tree", "Particle tree holding information about SBND final state topologies");

  event_tree->Branch(   "test_e", &test_e, "test_e/F");
  particle_tree->Branch("test_p", &test_p, "test_p/F");
  
  event_tree->SetDirectory(0);
  particle_tree->SetDirectory(0);

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
  std::cout << std::endl;
  
  std::cout << " Total number of tracks                            : ";
  std::cout << all_tracks;
  std::cout << std::endl;

  std::cout << " Fraction of fiducial, contained tracks            : ";
  std::cout << fiducial_contained_tracks / double(all_tracks);
  std::cout << std::endl;

  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  std::cout << "=================================================================================" << std::endl;

  // Print the tree, write the file, close
  TFile file("/sbnd/app/users/rsjones/LArSoft_v06_63_00/LArSoft-v06_63_00/srcs/recoparameters/recoparameters/output_files/tree_test.root", "RECREATE");
  event_tree->Write();
  particle_tree->Write();
  file.Write();
  file.Close();

}

void pndr::AnalysisNTuple::reconfigure(fhicl::ParameterSet const & p)
{
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
