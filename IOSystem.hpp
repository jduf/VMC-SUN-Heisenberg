#ifndef DEF_IOSYSTEM
#define DEF_IOSYSTEM

#include "RSTFile.hpp"

class IOSystem{
	public:
		IOSystem(std::string const& filename = "");
		IOSystem() = default;
		virtual ~IOSystem(){};
		/*!Forbids constructors*/
		IOSystem(IOSystem const&) = delete;
		IOSystem(IOSystem&&) = delete;

		/*!Set this to the value contained in t*/
		void set_IOSystem(IOSystem* t);
		/*!Set .jdbin output file and initialize its header*/
		void init_output_file(IOFiles& output);

		/*!Returns the filename (only usefull for mc via CreateSystem)*/
		std::string const& get_filename() const { return filename_; };
		/*!Returns the path (only usefull for mc via CreateSystem)*/
		std::string const& get_path() const { return path_; }

		/*!Call the corresponding virtual std::string extract_...()*/
		std::string analyse(unsigned int const& level);

	protected:
		std::string sim_;		//!< sim directory name
		std::string info_;		//!< directory name where the .rst and .html files will be saved
		std::string analyse_;	//!< directory name where the .gp dat .dat files will be saved 
		std::string path_;		//!< temporary variable of the path leading to dir_
		std::string dir_;		//!< temporary variable containing filename_
		std::string filename_;	//!< whole filename of the .jdbin file

		IOFiles* read_;			//!< binary file from which data are extracted
		IOFiles* jd_write_;		//!< binary file in which data are written in sim_/path_
		IOFiles* data_write_;	//!< textfile in which data are written saved in analyse_/path_
		RSTFile* rst_file_;		//!< rst file saved in info_/path_/dir_

		virtual std::string extract_level_1(){return filename_;}
		virtual std::string extract_level_2(){return filename_;}
		virtual std::string extract_level_3(){return filename_;}
		virtual std::string extract_level_4(){return filename_;}
		virtual std::string extract_level_5(){return filename_;}
		virtual std::string extract_level_6(){return filename_;}
		virtual std::string extract_level_7(){return filename_;}

	private:
		/*!Forbids assigment*/
		IOSystem& operator=(IOSystem const& ios);
};
#endif
