////////////////////////////////////////////////////////////////////////
// Class:       GunParticlePID
// Plugin Type: analyzer (art v2_08_04)
// File:        GunParticlePID_module.cc
//
// Generated at Mon Jan 22 06:06:15 2018 by Rhiannon Jones using cetskelgen
// from cetlib version v3_01_01.
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

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace pndr {
  class GunParticlePID;
}


class pndr::GunParticlePID : public art::EDAnalyzer {
public:
  explicit GunParticlePID(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GunParticlePID(GunParticlePID const &) = delete;
  GunParticlePID(GunParticlePID &&) = delete;
  GunParticlePID & operator = (GunParticlePID const &) = delete;
  GunParticlePID & operator = (GunParticlePID &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  // Declare member data here.
  // Member data for the type of gun particle, for the file name
  std::string m_particleType;
  
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
  int contained_track;

  // NTuple
  TNtuple *fNt_pid;

};


pndr::GunParticlePID::GunParticlePID(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{

  this->reconfigure(p);

}

void pndr::GunParticlePID::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  // Get the MCTruth information 
  art::Handle< std::vector< simb::MCTruth > > mct_handle;
  e.getByLabel("generator", mct_handle );
  int mct_size = mct_handle->size();

  // Get tracks and check they remain in the fiducial volume 
  art::Handle< std::vector< recob::Track > > trk_handle;
  e.getByLabel("pmalgtrackmaker", trk_handle );
  int trk_size = trk_handle->size();

  bool contained = true;

  if(mct_handle.isValid() && mct_size && trk_handle.isValid() && trk_size) {

    // Loop over the truth info
    for(auto const& mct : (*mct_handle)) {

      // Check the particle vertex is within the fiducial volume
      for(int i = 0; i < mct.NParticles(); ++i) {
    
        float truth_vtx_x = mct.GetParticle(i).Vx();
        float truth_vtx_y = mct.GetParticle(i).Vy();
        float truth_vtx_z = mct.GetParticle(i).Vz();
        
        // Check that the true vertex is within the fiducial volume
        if (    (truth_vtx_x > (m_detectorHalfLengthX - m_coordinateOffsetX - m_selectedBorderX)) 
             || (truth_vtx_x < (-m_coordinateOffsetX + m_selectedBorderX)) 
             || (truth_vtx_y > (m_detectorHalfLengthY - m_coordinateOffsetY - m_selectedBorderY)) 
             || (truth_vtx_y < (-m_coordinateOffsetY + m_selectedBorderY)) 
             || (truth_vtx_z > (m_detectorHalfLengthZ - m_coordinateOffsetZ - m_selectedBorderZ)) 
             || (truth_vtx_z < (-m_coordinateOffsetZ + m_selectedBorderZ))) contained = false;
    
      }
    }
    
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
  
    }
  }
  
  if( contained ){
  
    contained_track++;

    // Access PFParticles and their associated tracks. 
    // See if the associated tracks have associated calorimetry
    
    // Start to look at PID
    /*
    art::Handle< std::vector< anab::ParticleID > > pid_handle;
    e.getByLabel("chi2pid", pid_handle );
    int pid_size = pid_handle->size();
    */
    
    // Tracks
    art::Handle< std::vector< recob::Track > > trk_handle;
    e.getByLabel("pmalgtrackmaker", trk_handle );
    int trk_size = trk_handle->size();

    // Initialise PID handle associated with the tracks
    art::FindMany< anab::ParticleID > fmpid( trk_handle, e, "chi2pid" );
        
    // Check that we can access the PID information
    if( trk_size && trk_handle.isValid() ){
  
      for( int i = 0; i < trk_size; ++i ){

        //art::Ptr< recob::Track > trk( trk_handle, i );

        // Define PID handle
        std::vector<const anab::ParticleID*> pid_assn = fmpid.at(i);
 
        // Loop over PID association
        for ( size_t j = 0; j < pid_assn.size(); ++j ){

          if (!pid_assn[j]) continue;
          if (!pid_assn[j]->PlaneID().isValid) continue;
            
          // Get the plane number
          int planenum = pid_assn[j]->PlaneID().Plane;

          // Should only be 1,2 or 3
          if (planenum<0||planenum>2) continue;

          // Only look at the collection plane, since this is where the dEdx
          // is acquired and we need this for the PIDA values
          if( planenum == 2 ){
            
            int ndf       = pid_assn[j]->Ndf();
            double chi2Pr = pid_assn[j]->Chi2Proton();
            double chi2Mu = pid_assn[j]->Chi2Muon();
            double chi2Pi = pid_assn[j]->Chi2Pion();
            double chi2Ka = pid_assn[j]->Chi2Kaon();
            double pida   = pid_assn[j]->PIDA();

            // Fill the nTuple
            fNt_pid->Fill(ndf, chi2Pr, chi2Mu, chi2Pi, chi2Ka, pida);

          }
        }
      }
    }
  }
}

void pndr::GunParticlePID::beginJob()
{
  // Implementation of optional member function here.
  contained_track = 0;
  
  fNt_pid = new TNtuple( "fNt_pid",  "Particle gun PID information", "Ndf:Chi2Proton:Chi2Muon:Chi2Pion:Chi2Kaon:PIDA" );
  fNt_pid->SetDirectory(0);

}

void pndr::GunParticlePID::endJob()
{
  // Implementation of optional member function here.
  std::cout << " ================================================================================= " << std::endl;
  std::cout << " Number of contained tracks : " << contained_track << std::endl;
  std::cout << " Particle type              : " << m_particleType  << std::endl;

  std::stringstream name;
  name.clear();

  char file_path[1024];
  std::string temp;

  name << "/sbnd/app/users/rsjones/LArSoft_v06_63_00/LArSoft-v06_63_00/srcs/recoparameters/recoparameters/output_files/" << m_particleType << "_pid_information.root";
  temp = name.str();

  strcpy( file_path, temp.c_str() );  
  
  TFile *f = new TFile( file_path, "RECREATE" );

  fNt_pid->Write();

  f->Close();

  delete f;
  delete fNt_pid;

}

void pndr::GunParticlePID::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  // Get the particle type
  m_particleType        = p.get<std::string>("ParticleType");

  // Get the fiducial volume 
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

DEFINE_ART_MODULE(pndr::GunParticlePID)
