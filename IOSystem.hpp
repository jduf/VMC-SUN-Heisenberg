#ifndef DEF_IOSYSTEM
#define DEF_IOSYSTEM

#include "RSTFile.hpp"

class IOSystem{
	public:
		IOSystem(std::string const& filename);
		IOSystem(IOSystem const& a);
		virtual ~IOSystem(){};

		virtual std::string analyse(IOSystem const& t);

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

	private:
		IOSystem& operator=(IOSystem a);
		/*!Copy-And-Swap Idiom*/
		void swap_to_assign(IOSystem& t1, IOSystem& t2);
};
#endif
