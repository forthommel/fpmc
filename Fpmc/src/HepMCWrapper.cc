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
    if ( !HepMC::HEPEVT_Wrapper::HEPEVT_to_GenEvent( hepMCEvt_.get() ) )
      throw std::runtime_error( "Failed to fetch the HEPEVT block!" );
    std::cout << HepMC::HEPEVT_Wrapper::number_entries() << " >> " << hepMCEvt_->vertices().size() << std::endl;
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
