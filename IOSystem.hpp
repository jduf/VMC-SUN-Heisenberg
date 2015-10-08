#ifndef DEF_IOSYSTEM
#define DEF_IOSYSTEM

#include "RSTFile.hpp"

class IOSystem{
	public:
		IOSystem(std::string const& filename, std::vector<std::string> names);
		IOSystem() = default;
		virtual ~IOSystem() = default;
		/*{Forbidden*/
		IOSystem(IOSystem const&) = delete;
		IOSystem(IOSystem&&) = delete;
		IOSystem& operator=(IOSystem ios) = delete;
		/*}*/

		/*!Set this to the value contained in t*/
		void set_IOSystem(IOSystem const* const t);

		/*!Returns the filename (only usefull for mc via CreateSystem)*/
		std::string const& get_filename() const { return filename_; };
		/*!Returns the path (only usefull for mc via CreateSystem)*/
		std::string const& get_path() const { return path_; }

		/*!Call the corresponding virtual std::string extract_...()*/
		std::string analyse(unsigned int const& level);

	protected:
		IOFiles* read_		  = NULL;//!< binary file from which data are extracted
		IOFiles* jd_write_	  = NULL;//!< binary file in which data are written in sim_/path_
		IOFiles* data_write_  = NULL;//!< textfile in which data are written saved in analyse_/path_
		RSTFile* rst_file_	  = NULL;//!< rst file saved in info_/path_/dir_

		std::string sim_	  = "sim/";		//!< sim directory name
		std::string info_	  = "info/";	//!< directory name where the .rst and .html files will be saved
		std::string analyse_  = "analyse/";	//!< directory name where the .gp dat .dat files will be saved 
		std::string path_	  = "";			//!< temporary variable of the path leading to dir_
		std::string dir_	  = "";			//!< temporary variable containing filename_
		std::string filename_ = "";			//!< whole filename of the .jdbin file

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
