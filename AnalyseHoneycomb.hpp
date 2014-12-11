#ifndef DEF_ANALYSEHONEYCOMB 
#define DEF_ANALYSEHONEYCOMB 

#include "Analyse.hpp"

class AnalyseHoneycomb : public Analyse{
	public:
		AnalyseHoneycomb(std::string const& sim);
		~AnalyseHoneycomb(){}

	protected:
		void open_files();
		void close_files();

		std::string extract_level_6();
};
#endif
