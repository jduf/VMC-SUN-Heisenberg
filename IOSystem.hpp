#ifndef DEF_IOSYSTEM
#define DEF_IOSYSTEM

#include "RSTFile.hpp"

class IOSystem{
	public:
		IOSystem(std::string const& filename);
		IOSystem(IOSystem const& a);
		virtual ~IOSystem(){}

		std::string analyse(unsigned int const& level, IOSystem* t);
		IOFiles* open_and_get_jd_write();
		void close_jd_write(){delete jd_write_;} 

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

		virtual std::string extract_level_1();
		virtual std::string extract_level_2();
		virtual std::string extract_level_3();
		virtual std::string extract_level_4();
		virtual std::string extract_level_5();
		virtual std::string extract_level_6();

	private:
		IOSystem& operator=(IOSystem a);
		/*!Copy-And-Swap Idiom*/
		void swap_to_assign(IOSystem& t1, IOSystem& t2);
};
#endif
