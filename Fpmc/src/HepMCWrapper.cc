#include "Fpmc/HepMCWrapper.h"
#include "Fpmc/FpmcParameters.h"

#include "HepMC/GenEvent.h"

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
# include "HepMC/HEPEVT_Wrapper.h"
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

using namespace HepMC;
  const HepMC::GenEvent&
  HepMCWrapper::event()
  {
    //----- start by generating the next event with FPMC

    if ( !Fpmc::next() || hwevnt_.IERROR )
      throw std::runtime_error( "Failed to generate the next event!" );

#ifdef HEPMC_VERSION2
    hepMCEvt_.reset( *conv_.read_next_event() );
#else
    for ( unsigned short i = 8; i <= HepMC::HEPEVT_Wrapper::number_entries(); ++i ) {
      HepMC::HEPEVT_Wrapper::set_position( i-7,
        HepMC::HEPEVT_Wrapper::x( i ), HepMC::HEPEVT_Wrapper::y( i ), HepMC::HEPEVT_Wrapper::z( i ),
        HepMC::HEPEVT_Wrapper::t( i ) );
      HepMC::HEPEVT_Wrapper::set_momentum( i-7,
        HepMC::HEPEVT_Wrapper::px( i ), HepMC::HEPEVT_Wrapper::py( i ), HepMC::HEPEVT_Wrapper::pz( i ),
        HepMC::HEPEVT_Wrapper::e( i ) );
      HepMC::HEPEVT_Wrapper::set_mass( i-7, HepMC::HEPEVT_Wrapper::m( i ) );
      HepMC::HEPEVT_Wrapper::set_id( i-7, HepMC::HEPEVT_Wrapper::id( i ) );
      HepMC::HEPEVT_Wrapper::set_status( i-7, HepMC::HEPEVT_Wrapper::status( i ) );
      HepMC::HEPEVT_Wrapper::set_parents( i-7,
        std::max( HepMC::HEPEVT_Wrapper::first_parent( i )-7, 0 ),
        std::max( HepMC::HEPEVT_Wrapper:: last_parent( i )-7, 0 ) );
    }
    HepMC::HEPEVT_Wrapper::set_number_entries( HepMC::HEPEVT_Wrapper::number_entries()-7 );
    //--- fix parentage for incoming 'partons'
    for ( unsigned short i = 1; i <= 2; ++i ) {
      HepMC::HEPEVT_Wrapper::set_status( i, 3 );
      HepMC::HEPEVT_Wrapper::set_parents( i, 0, 0 );
      HepMC::HEPEVT_Wrapper::set_children( i, 3, 0 );
    }
    //--- fix parentage for two-'parton' system
    HepMC::HEPEVT_Wrapper::set_children( 3, 4, HepMC::HEPEVT_Wrapper::number_entries() );
    //--- central event
    for ( unsigned short i = 4; i <= HepMC::HEPEVT_Wrapper::number_entries(); ++i )
      HepMC::HEPEVT_Wrapper::set_parents( i, 3, 0 );

    if ( !HepMC::HEPEVT_Wrapper::HEPEVT_to_GenEvent( hepMCEvt_.get() ) )
      throw std::runtime_error( "Failed to fetch the HEPEVT block!" );

    HepMC::HEPEVT_Wrapper::zero_everything();
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
}

