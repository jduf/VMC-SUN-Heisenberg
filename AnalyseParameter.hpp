#ifndef DEF_ANALYSEPARAMETER
#define DEF_ANALYSEPARAMETER

#include "Analyse.hpp"

class AnalyseParameter : public Analyse{
	public:
		AnalyseParameter(){}
		~AnalyseParameter(){}

	protected:
		void open_files(std::string const& jdfile, std::string const& datafile, Directory const& d);
		void close_files();

		void extract_level_5();
		void extract_level_4();
		void extract_level_3();
		void extract_level_2();
};
#endif
