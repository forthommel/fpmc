#include "Fpmc/HepMCWrapper.h"
#include "Fpmc/FpmcParameters.h"

#include "HepMC/GenEvent.h"
#define HEPMC_HEPEVT_NMXHEP 4000
#include "HepMC/HEPEVT_Wrapper.h"

#include <string>
#include <sstream>
#include <stdexcept>
#include <iostream>

#ifdef HEPMC_VERSION2
# include "HepMC/PdfInfo.h"
#else // HepMC v>=3
# include "HepMC/WriterAscii.h"
# include "HepMC/GenPdfInfo.h"
# include "HepMC/GenParticle.h"
# include "HepMC/GenVertex.h"
# include <set>
# include <map>
extern "C" struct HEPEVT hepevt_;
#endif

namespace fpmc
{
  HepMCWrapper::HepMCWrapper( const char* card ) :
    Fpmc( card ), hepMCEvt_( new HepMC::GenEvent ), hepMCVerbosity_( true )
  {
#ifndef HEPMC_VERSION2
    HepMC::HEPEVT_Wrapper::set_hepevt_address( (char*)&hepevt_ );
#endif
    //params_.dump();
  }

  HepMCWrapper::~HepMCWrapper()
  {}

  const HepMC::GenEvent&
  HepMCWrapper::event()
  {
    //----- start by generating the next event with FPMC

    if ( !Fpmc::next() || hwevnt_.IERROR )
      throw std::runtime_error( "Failed to generate the next event!" );

#ifdef HEPMC_VERSION2
    hepMCEvt_.reset( *conv_.read_next_event() );
#else
    //HepMC::HEPEVT_Wrapper::zero_everything();
    HepMC::HEPEVT_Wrapper::print_hepevt();
    for ( unsigned short i = 8; i <= HepMC::HEPEVT_Wrapper::number_entries(); ++i ) {
      HepMC::HEPEVT_Wrapper::set_position( i-7, HepMC::HEPEVT_Wrapper::x( i ), HepMC::HEPEVT_Wrapper::y( i ), HepMC::HEPEVT_Wrapper::z( i ), HepMC::HEPEVT_Wrapper::t( i ) );
      HepMC::HEPEVT_Wrapper::set_momentum( i-7, HepMC::HEPEVT_Wrapper::px( i ), HepMC::HEPEVT_Wrapper::py( i ), HepMC::HEPEVT_Wrapper::pz( i ), HepMC::HEPEVT_Wrapper::e( i ) );
      HepMC::HEPEVT_Wrapper::set_mass( i-7, HepMC::HEPEVT_Wrapper::m( i ) );
      HepMC::HEPEVT_Wrapper::set_id( i-7, HepMC::HEPEVT_Wrapper::id( i ) );
      HepMC::HEPEVT_Wrapper::set_status( i-7, HepMC::HEPEVT_Wrapper::status( i ) );
      HepMC::HEPEVT_Wrapper::set_parents( i-7, std::max( HepMC::HEPEVT_Wrapper::first_parent( i )-7, 0 ), std::max( HepMC::HEPEVT_Wrapper::last_parent( i )-7, 0 ) );
    }
    HepMC::HEPEVT_Wrapper::set_number_entries( HepMC::HEPEVT_Wrapper::number_entries()-7 );
    //--- fix parentage for incoming 'partons'
    int part_status = 101;
    for ( unsigned short i = 1; i <= 2; ++i ) {
      //HepMC::HEPEVT_Wrapper::set_status( i, part_status++ );
      HepMC::HEPEVT_Wrapper::set_status( i, 3 );
      HepMC::HEPEVT_Wrapper::set_parents( i, 0, 0 );
      HepMC::HEPEVT_Wrapper::set_children( i, 3, 0 );
    }
    //--- fix parentage for two-'parton' system
    //HepMC::HEPEVT_Wrapper::set_status( 3, part_status++ );
    HepMC::HEPEVT_Wrapper::set_children( 3, 4, HepMC::HEPEVT_Wrapper::number_entries() );
    //--- central event
    for ( unsigned short i = 4; i <= HepMC::HEPEVT_Wrapper::number_entries(); ++i ) {
      HepMC::HEPEVT_Wrapper::set_parents( i, 3, 0 );
      //HepMC::HEPEVT_Wrapper::set_status( i, 1 );
    }
    HepMC::HEPEVT_Wrapper::print_hepevt();
    //if ( !HepMC::HEPEVT_Wrapper::fix_daughters() )
      //throw std::runtime_error( "Failed to fix particles parentage!" );
    //  std::cout << "prout" << std::endl;
    for ( unsigned short i = 1; i <= HepMC::HEPEVT_Wrapper::number_entries(); ++i )
      std::cout << HepMC::HEPEVT_Wrapper::event_number() << "|" << i << ": " << HepMC::HEPEVT_Wrapper::status( i ) << "|" << HepMC::HEPEVT_Wrapper::number_parents( i ) << "|" << HepMC::HEPEVT_Wrapper::number_children( i ) << "|" << HepMC::HEPEVT_Wrapper::number_children_exact( i ) << std::endl;

    HepMC::GenEvent evt;
//    if ( !HepMC::HEPEVT_Wrapper::HEPEVT_to_GenEvent( hepMCEvt_.get() ) )
    if ( !HepMC::HEPEVT_Wrapper::HEPEVT_to_GenEvent( &evt ) )
      throw std::runtime_error( "Failed to fetch the HEPEVT block!" );
    //std::cout << HepMC::HEPEVT_Wrapper::number_entries() << " >> " << hepMCEvt_->vertices().size() << std::endl;

/*std::map<HepMC::GenParticlePtr,int > hepevt_particles;
std::map<int,HepMC::GenParticlePtr > particles_index;
std::map<HepMC::GenVertexPtr,std::pair<std::set<int>,std::set<int> > > hepevt_vertices;
std::map<int,HepMC::GenVertexPtr > vertex_index;
for ( int i = 1; i <= HepMC::HEPEVT_Wrapper::number_entries(); i++ ) {
  HepMC::GenParticlePtr p=std::make_shared<HepMC::GenParticle>();
  p->set_momentum(HepMC::FourVector( HepMC::HEPEVT_Wrapper::px(i), HepMC::HEPEVT_Wrapper::py(i), HepMC::HEPEVT_Wrapper::pz(i), HepMC::HEPEVT_Wrapper::e(i) ));
  p->set_status(HepMC::HEPEVT_Wrapper::status(i));
  p->set_pid(HepMC::HEPEVT_Wrapper::id(i)); //Confusing!
  p->set_generated_mass( HepMC::HEPEVT_Wrapper::m(i));
  hepevt_particles[p]=i;
  particles_index[i]=p;
  HepMC::GenVertexPtr v=std::make_shared<HepMC::GenVertex>();
  v->set_position(HepMC::FourVector( HepMC::HEPEVT_Wrapper::x(i), HepMC::HEPEVT_Wrapper::y(i), HepMC::HEPEVT_Wrapper::z(i), HepMC::HEPEVT_Wrapper::t(i)));
  v->add_particle_out(p);
  std::set<int> in, out;
  out.insert(i);
  vertex_index[i]=v;
  hepevt_vertices[v]= { in, out };
}
for (auto& it1 : hepevt_particles)
  for (auto& it2 : hepevt_particles)
    if (HepMC::HEPEVT_Wrapper::first_parent(it2.second)<=it1.second&&it1.second<=HepMC::HEPEVT_Wrapper::last_parent(it2.second)) //I'm you father, Luck!
      hepevt_vertices[it2.first->production_vertex()].first.insert(it1.second);
for ( int i = 1; i <= HepMC::HEPEVT_Wrapper::number_entries(); i++ ) vertex_index[i]->remove_particle_out(particles_index[i]);
std::map<HepMC::GenVertexPtr,std::pair<std::set<int>,std::set<int> > > final_vertices_map;
for (auto& vs : hepevt_vertices) {
  if ((final_vertices_map.size()==0)||(vs.second.first.size()==0&&vs.second.second.size()!=0)) { final_vertices_map.insert(vs);  continue; } //Always insert particles out of nowhere
  std::map<HepMC::GenVertexPtr,std::pair<std::set<int>,std::set<int> > >::iterator  v2;
  for (v2=final_vertices_map.begin(); v2!=final_vertices_map.end(); v2++) if (vs.second.first==v2->second.first) {v2->second.second.insert(vs.second.second.begin(),vs.second.second.end()); break;}
  if (v2==final_vertices_map.end()) final_vertices_map.insert(vs);
}
std::vector<HepMC::GenParticlePtr> final_particles;
std::set<int> used;
for (auto& it : final_vertices_map) {
  HepMC::GenVertexPtr v=it.first;
  std::set<int> in=it.second.first, out=it.second.second;
  used.insert(in.begin(),in.end());
  used.insert(out.begin(),out.end());
  for (const auto& el : in) v->add_particle_in(particles_index[el]);
  if (in.size()!=0) for (const auto& el : out) v->add_particle_out(particles_index[el]);
}
for (auto& el : used) final_particles.push_back(particles_index[el]);
evt.add_tree( final_particles );*/

    std::cout << HepMC::HEPEVT_Wrapper::number_entries() << " >> " << hepMCEvt_->particles().size() << "|" << hepMCEvt_->vertices().size() << std::endl;
    std::cout << HepMC::HEPEVT_Wrapper::number_entries() << " >> " << evt.particles().size() << "|" << evt.vertices().size() << std::endl;

#endif

#ifdef HEPMC_VERSION2
    hepMCEvt_->set_event_number( event_-1 );
    hepMCEvt_->set_signal_process_id( params_.processId() );
    hepMCEvt_->set_event_scale( -1. );
#endif
    hepMCEvt_->weights().push_back( hwevnt_.EVWGT );

#ifdef HEPMC_VERSION2
    HepMC::PdfInfo pdfInfo;
    pdfInfo.set_x1( hwhard_.XX[0] );
    pdfInfo.set_x2( hwhard_.XX[1] );
    pdfInfo.set_scalePDF( hwhard_.EMSCA ); //FIXME
    hepMCEvt_->set_pdf_info( pdfInfo );
#else
    auto pdf_info = std::make_shared<HepMC::GenPdfInfo>();
    pdf_info->set( 0, 0, hwhard_.XX[0], hwhard_.XX[1], hwhard_.EMSCA, 0, 0 ); //FIXME
    hepMCEvt_->set_pdf_info( pdf_info );
#endif

    /*HepMC::GenParticle* incomingParton = NULL;
    HepMC::GenParticle* targetParton = NULL;
    // find incoming parton (first entry with IST=121)
    for(HepMC::GenEvent::particle_const_iterator it = event()->particles_begin(); (it != event()->particles_end() && incomingParton==NULL); it++)
      if((*it)->status()==121) incomingParton = (*it);

    // find target parton (first entry with IST=122)
    for(HepMC::GenEvent::particle_const_iterator it = event()->particles_begin(); (it != event()->particles_end() && targetParton==NULL); it++)
      if((*it)->status()==122) targetParton = (*it);*/

#ifdef HEPMC_VERSION2
    //******** Verbosity ********
    if ( event_ <= maxEventsToPrint_ && hepMCVerbosity_ ) {
      // Prints HepMC event
      dbg_ << "\n----------------------" << std::endl
           << "Event process id = " << hepMCEvt_->signal_process_id() << std::endl;
      hepMCEvt_->print();
    }
#endif

    return *hepMCEvt_;
  }

  void
  HepMCWrapper::write( const char* out )
  {
    if ( !hepMCEvt_.get() ) return;

#ifdef HEPMC_VERSION2
    std::ofstream output( out );
    hepMCEvt_->write( output );
    output.close();
#else // HepMC v>=3
    HepMC::WriterAscii output( out );
    output.write_event( *hepMCEvt_ );
    output.close();
#endif
  }
}
