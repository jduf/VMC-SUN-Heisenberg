#ifndef DEF_ANALYSEMAGNETIZATION
#define DEF_ANALYSEMAGNETIZATION

#include "Analyse.hpp"

class AnalyseMagnetization : public Analyse{
	public:
		AnalyseMagnetization(std::string const& sim);
		~AnalyseMagnetization(){}

	protected:
		void close_files();
		void open_files();

		std::string extract_level_6();
		std::string extract_level_5();
		std::string extract_level_4();
		std::string extract_level_3();
		std::string extract_level_2();
};
#endif
