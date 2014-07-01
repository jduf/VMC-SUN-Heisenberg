#ifndef DEF_ANALYSEMAGNETIZATION
#define DEF_ANALYSEMAGNETIZATION

#include "Analyse.hpp"

class AnalyseMagnetization : public Analyse{
	public:
		AnalyseMagnetization(){}
		~AnalyseMagnetization(){}

	protected:
		void close_files();
		void open_files(std::string const& jdfile, std::string const& datafile, Directory const& d);
		void extract_level_5();
		void extract_level_4();
		void extract_level_3();
};
#endif
