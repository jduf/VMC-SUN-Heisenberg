#ifndef DEF_ANALYSELONGRANGECORRELATION 
#define DEF_ANALYSELONGRANGECORRELATION 

#include "Analyse.hpp"

class AnalyseLongRangeCorrelation : public Analyse{
	public:
		AnalyseLongRangeCorrelation(std::string const& sim);
		~AnalyseLongRangeCorrelation(){}

	protected:
		void open_files();
		void close_files();

		std::string extract_level_7();
};
#endif
