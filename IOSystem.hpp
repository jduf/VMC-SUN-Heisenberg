#ifndef DEF_IOSYSTEM
#define DEF_IOSYSTEM

#include "RSTFile.hpp"
#include "Vector.hpp"
#include "List.hpp"

class IOSystem{
	public:
		IOSystem(std::string const& filename, std::string const& sim, std::string const& info, std::string const& analyse, std::string const& path, std::string const& dir, RSTFile* const rst_file);
		virtual ~IOSystem() = default;
		/*{Forbidden*/
		IOSystem(IOSystem const&) = delete;
		IOSystem(IOSystem&&) = delete;
		IOSystem& operator=(IOSystem) = delete;
		/*}*/

		/*!Set this to the value contained in t*/
		void set_IOSystem(IOSystem const* const ios);

		/*!Returns the filename (only useful for mc via CreateSystem)*/
		std::string const& get_filename() const { return filename_; }
		/*!Returns the path (only useful for mc via CreateSystem)*/
		std::string const& get_path() const { return path_; }
		/*Returns the total path from info*/
		std::string get_info_path() const { return info_+path_+dir_; }
		/*Returns information about the System stored in system_info_*/
		RST const& get_system_info() const { return system_info_; }

		/*!Call the corresponding virtual std::string extract_...()*/
		std::string analyse(unsigned int const& level);

	protected:
		IOSystem(std::string const& filename, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc, Vector<unsigned int> const& ref);
		IOSystem() = default;

		std::string sim_	  = "sim/";		//!< sim directory name
		std::string info_	  = "info/";	//!< directory name where the .rst and .html files will be saved
		std::string analyse_  = "analyse/";	//!< directory name where the .gp dat .dat files will be saved
		std::string path_	  = "";			//!< temporary variable of the path leading to dir_
		std::string dir_	  = "";			//!< temporary variable containing filename_
		std::string filename_ = "";			//!< whole filename of the .jdbin file

		IOFiles* read_		  = NULL;//!< binary file from which data are extracted
		IOFiles* jd_write_	  = NULL;//!< binary file in which data are written in sim_/path_
		IOFiles* data_write_  = NULL;//!< textfile in which data are written saved in analyse_/path_
		RSTFile* rst_file_	  = NULL;//!< rst file saved in info_/path_/dir_

		RST system_info_;//!< information about the wavefunction

		virtual std::string extract_level_1(){ return filename_; }
		virtual std::string extract_level_2(){ return filename_; }
		virtual std::string extract_level_3(){ return filename_; }
		virtual std::string extract_level_4(){ return filename_; }
		virtual std::string extract_level_5(){ return filename_; }
		virtual std::string extract_level_6(){ return filename_; }
		virtual std::string extract_level_7(){ return filename_; }
		virtual std::string extract_level_8(){ return filename_; }
		virtual std::string extract_level_9(){ return filename_; }
};
#endif
