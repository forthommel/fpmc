#ifndef Fpmc_ArgsParser_h
#define Fpmc_ArgsParser_h

#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <algorithm>
#include <stdexcept>

namespace fpmc
{
  class ArgsParser
  {
    public:
      struct Parameter {
        Parameter() {}
        Parameter( std::string key, std::string descr ) : key( key ), description( descr ) {}
        const bool operator=( const Parameter& rhs ) const { return ( key == rhs.key ); }
        std::string key, description;
      };

    public:
      ArgsParser( int argc, char* argv[], const std::vector<Parameter>& required_parameters, const std::vector<Parameter>& optional_parameters ) :
        help_str_( { "-h", "--help" } ),
        required_params_( required_parameters ), optional_params_( optional_parameters ) {
        command_name_ = argv[0];
        if ( argc >1  ) {
          args_.resize( argc-1 );
          std::copy( argv+1, argv+argc, args_.begin() );
        }
        for ( const auto& str : help_str_ ) {
          if ( find( args_.begin(), args_.end(), str ) != args_.end() ) {
            print_help();
            exit( 0 );
          }
        }
      }
      /// Read required parameters
      std::map<std::string,std::string> required_parameters() const {
        std::map<std::string,std::string> out;
        std::ostringstream oss; oss << "The following parameter(s) was/were not set:\n";
        bool valid = true;
        for ( const auto& par : required_params_ ) {
          std::ostringstream par_ss; par_ss << "--" << par.key;
          const auto key = find( args_.begin(), args_.end(), par_ss.str() );
          if ( key == args_.end() ) {
            valid = false;
            oss << "(*) " << par.key << ": " << par.description << "\n";
            continue;
          }

          const auto value = key + 1;
          if ( value == args_.end() ) {
            std::ostringstream oss; oss << "Invalid value for parameter: " << par.key << std::endl;
            throw Exception( oss.str() );
          }
          out[par.key] = *value;
        }
        if ( !valid ) {
          print_help();
          throw Exception( oss.str() );
        }
        return out;
      }
      /// Read optional parameters
      std::map<std::string,std::string> optional_parameters() const {
        std::map<std::string,std::string> out;
        for ( const auto& par : optional_params_ ) {
          std::ostringstream par_ss; par_ss << "--" << par.key;
          const auto key = find( args_.begin(), args_.end(), par_ss.str() );
          if ( key == args_.end() ) continue; // Parameter not set

          const auto value = key + 1;
          if ( value == args_.end() ) {
            std::ostringstream oss; oss << "Invalid value for parameter: " << par.key << std::endl;
            throw Exception( oss.str() );
          }
          out[par.key] = *value;
        }
        return out;
      }
      /// Show usage
      void print_help() const {
        std::ostringstream oss;
        oss << "Usage: " << command_name_ << " ";
        for ( const auto& par : required_params_ ) {
          oss << "--" <<  par.key << " <" << par.key << "> ";
        }
        for ( const auto& par : optional_params_ ) {
          oss << "--" <<  par.key << " [" << par.key << "] ";
        }
        oss << std::endl;
        std::cout << oss.str(); 
      }
    private:
      struct Exception : public std::runtime_error
      {
        Exception( const char* text ) : std::runtime_error( text ) {}
        Exception( const std::string& text ) : std::runtime_error( text ) {}
      };
      std::string command_name_;
      const std::vector<std::string> help_str_;
      std::vector<Parameter> required_params_, optional_params_;
      std::vector<std::string> args_;
  };
}

#endif

