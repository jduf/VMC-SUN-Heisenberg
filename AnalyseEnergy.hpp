#ifndef DEF_ANALYSEENERGY
#define DEF_ANALYSEENERGY

#include "Analyse.hpp"

class AnalyseEnergy : public Analyse{
	public:
		AnalyseEnergy(std::string const& sim, std::string const& path, unsigned int const& max_level, bool const& run_cmd);
		/*!Default destructor*/
		~AnalyseEnergy() = default;
		/*{Forbidden*/
		AnalyseEnergy() = delete;
		AnalyseEnergy(AnalyseEnergy const&) = delete;
		AnalyseEnergy(AnalyseEnergy&&) = delete;
		AnalyseEnergy& operator=(AnalyseEnergy) = delete;
		/*}*/

	protected:
		void open_files();
		void close_files();

		std::string extract_level_6();
		std::string extract_level_5();
		std::string extract_level_4();
		std::string extract_level_3();
		std::string extract_level_2();

	private:
		std::string extract_default();
};
#endif
