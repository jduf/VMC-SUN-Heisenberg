#ifndef DEF_ANALYSEPARAMETER
#define DEF_ANALYSEPARAMETER

#include "Analyse.hpp"

class AnalyseParameter : public Analyse{
	public:
		AnalyseParameter(std::string const& sim);
		~AnalyseParameter(){}

	protected:
		void open_files();
		void close_files();

		std::string extract_level_7();
		std::string extract_level_6();
		std::string extract_level_5();
		std::string extract_level_4();
		std::string extract_level_3();
		std::string extract_level_2();
};
#endif
