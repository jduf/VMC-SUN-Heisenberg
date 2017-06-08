#ifndef DEF_ANALYSEENERGY
#define DEF_ANALYSEENERGY

#include "Analyse.hpp"

class AnalyseEnergy : public Analyse{
	public:
		AnalyseEnergy(std::string const& sim, unsigned int const& max_level, unsigned int const& bash_file, unsigned int const& ref);
		/*!Default destructor*/
		~AnalyseEnergy();
		/*{Forbidden*/
		AnalyseEnergy() = delete;
		AnalyseEnergy(AnalyseEnergy const&) = delete;
		AnalyseEnergy(AnalyseEnergy&&) = delete;
		AnalyseEnergy& operator=(AnalyseEnergy) = delete;
		/*}*/

	protected:
		void open_files();
		void close_files();

		std::string extract_level_8();
		std::string extract_level_7();
		std::string extract_level_6();
		std::string extract_level_5();
		std::string extract_level_4();
		std::string extract_level_3();
		std::string extract_level_2();
};
#endif
