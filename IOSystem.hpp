#ifndef DEF_IOSYSTEM
#define DEF_IOSYSTEM

#include "RSTFile.hpp"

class IOSystem{
	public:
		IOSystem(std::string const& filename);
		virtual ~IOSystem(){}

		void set_IOSystem(IOSystem* t);
		std::string analyse(unsigned int const& level);

		/*!Returns the filename with the path and create the directory tree*/
		std::string get_filename() const;
		void init_output_file(IOFiles& output);

	protected:
		std::string sim_;
		std::string info_;
		std::string analyse_;
		std::string path_;
		std::string dir_;
		std::string filename_;

		IOFiles* read_;
		IOFiles* jd_write_;
		IOFiles* data_write_;
		RSTFile* rst_file_;

		virtual std::string extract_level_1(){return filename_;}
		virtual std::string extract_level_2(){return filename_;}
		virtual std::string extract_level_3(){return filename_;}
		virtual std::string extract_level_4(){return filename_;}
		virtual std::string extract_level_5(){return filename_;}
		virtual std::string extract_level_6(){return filename_;}
		virtual std::string extract_level_7(){return filename_;}

	private:
		/*!Forbids copy*/
		IOSystem(IOSystem const& ios);
		/*!Forbids assigment*/
		IOSystem& operator=(IOSystem a);
};
#endif
