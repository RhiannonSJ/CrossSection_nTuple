////////////////////////////////////////////////////////////////////////
// Class:       ParameterFinding
// Plugin Type: analyzer (art v2_08_04)
// File:        ParameterFinding_module.cc
//
// Generated at Tue Jan  9 13:43:33 2018 by Rhiannon Jones using cetskelgen
// from cetlib version v3_01_01.
//
// To start with, find out what parameters it is possible to obtain from 
// Pandora, PMA and any associations
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

namespace xsec {
  class ParameterFinding;
}


class xsec::ParameterFinding : public art::EDAnalyzer {
public:
  explicit ParameterFinding(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ParameterFinding(ParameterFinding const &) = delete;
  ParameterFinding(ParameterFinding &&) = delete;
  ParameterFinding & operator = (ParameterFinding const &) = delete;
  ParameterFinding & operator = (ParameterFinding &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  // Declare member data here.
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
  int primary;

};


xsec::ParameterFinding::ParameterFinding(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  
  this->reconfigure(p);

}

void xsec::ParameterFinding::analyze(art::Event const & e)
{
  
  // Implementation of required member function here.
  // Get the MCTruth information 
  art::Handle< std::vector< simb::MCTruth > > mct_handle;
  e.getByLabel("generator", mct_handle );
  int mct_size = mct_handle->size();

  if(mct_handle.isValid() && mct_size) {
  
    // Loop over the truth info
    for(auto const& mct : (*mct_handle)) {
 
      // Check the neutrino came from the beam
      if(mct.Origin() != simb::kBeamNeutrino) continue;

      // True neutrino vertex ( primary vertex position)
      double nu_x, nu_y, nu_z;

      nu_x = mct.GetNeutrino().Lepton().Vx();
      nu_y = mct.GetNeutrino().Lepton().Vy();
      nu_z = mct.GetNeutrino().Lepton().Vz();
        
      // Access PFParticles and their associated tracks. 
      // See if the associated tracks have associated calorimetry
      art::Handle< std::vector< recob::PFParticle > > pfp_handle;
      e.getByLabel("pandoraNu", pfp_handle );
      int pfp_size = pfp_handle->size();

      art::Handle< std::vector< recob::Vertex > > vtx_handle;
      e.getByLabel("pmalgtrackmaker", vtx_handle );
      int vtx_size = vtx_handle->size();

      // x,y,z of the primary reconstructed vertex
      double x = 0.;
      double y = 0.;
      double z = 0.;

      // Only loop over the pma vertices if the primary does exist
      bool primary_exists = false;

      // Get the location of the primary vertex from pandora 
      // Then find that vertex using PMA and get track and calorimetry associations
      if( pfp_handle.isValid() && pfp_size ){

        bool primary_found = false;

        // Get vertex associations
        art::FindMany< recob::Vertex  > fvtx( pfp_handle, e, "pandoraNu" );

        // Loop over PFParticles
        for( int i = 0; i < pfp_size; ++i ){ 

          art::Ptr< recob::PFParticle > pfp( pfp_handle, i );
        
          if( !primary_found ){
          
            // Find the primary and get the associated vertex
            if( pfp->IsPrimary() ){

              // Define vertex handle
              std::vector<const recob::Vertex*> vtx_assn = fvtx.at(i);
            
              for( unsigned int j = 0; j < vtx_assn.size(); ++j ){
              
                double xyz[3];

                // Set array to be current vertex position
                vtx_assn[j]->XYZ(xyz);

                // Set x, y, z of the primary
                x = xyz[0];
                y = xyz[1];
                z = xyz[2];

                primary_exists = true;
              
              }
            }
          }
        }
      }
      // Starting from the vertices given by PMA, find the primary and get 
      // track/calorimetry associations
      if( primary_exists && vtx_handle.isValid() && vtx_size ){
  
        std::vector< double > vertex_distances, true_distances;
        vertex_distances.clear();

        // Loop over vertices and find the primary
        for( int i = 0; i < vtx_size; ++ i ){
        
          art::Ptr< recob::Vertex > vtx( vtx_handle, i );
        
          double curr_x, curr_y, curr_z;
          double xyz[3];

          // Set array to be current vertex position
          vtx->XYZ(xyz);

          // Set x, y, z of the primary
          curr_x = xyz[0];
          curr_y = xyz[1];
          curr_z = xyz[2];
          vertex_distances.push_back( sqrt( pow( curr_x - x, 2 )    + pow( curr_y - y, 2 )    + pow( curr_z - z, 2 ) ) );
          true_distances.push_back(   sqrt( pow( curr_x - nu_x, 2 ) + pow( curr_y - nu_y, 2 ) + pow( curr_z - nu_z, 2 ) ) );

        }

        std::vector< double >::iterator min_r = std::min_element( std::begin( vertex_distances ), std::end( vertex_distances ) );
        double pma_pndr_distance               = vertex_distances[ std::distance( std::begin( vertex_distances ),  min_r ) ];
     
        // If the recostruction packages are reasonably consistent
        if( pma_pndr_distance < 5 && primary_exists && vtx_handle.isValid() && vtx_size ){
        
          // Loop over the vertices and get the track-calo associations
          for( int i = 0; i < vtx_size; ++ i ){
          
            art::Ptr< recob::Vertex > vtx( vtx_handle, i );
          
            double curr_x, curr_y, curr_z, curr_r;
            double xyz[3];

            // Set array to be current vertex position
            vtx->XYZ(xyz);

            // Set x, y, z of the primary
            curr_x = xyz[0];
            curr_y = xyz[1];
            curr_z = xyz[2];
            curr_r = sqrt( pow( curr_x - x, 2 )    + pow( curr_y - y, 2 )    + pow( curr_z - z, 2 ) );

            if( curr_r == pma_pndr_distance ){
            
              std::cout << " PMA - Pandora distance : " << curr_r << std::endl;
           
              primary++;

            }
          }
        }
      }
    }
  }
}

void xsec::ParameterFinding::beginJob()
{
  // Implementation of optional member function here.
  primary = 0;

}

void xsec::ParameterFinding::endJob()
{
  // Implementation of optional member function here.

  std::cout << " Number of times the primary is found nicely : " << primary;
  std::cout << std::endl;
}

void xsec::ParameterFinding::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
   std::vector< int > blankVect;
   std::vector< std::vector< int > > input;
 
   std::vector< int > selection1 = p.get< std::vector< int > >("Selection1",        blankVect);
   if ( selection1.size() != 0 ) input.push_back(selection1);
 
   std::vector< int > selection2 = p.get< std::vector< int > >("Selection2",        blankVect);
   if ( selection2.size() != 0 ) input.push_back(selection2);
 
   std::vector< int > selection3 = p.get< std::vector< int > >("Selection3",        blankVect);
   if ( selection3.size() != 0 ) input.push_back(selection3);
 
 
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

DEFINE_ART_MODULE(xsec::ParameterFinding)
