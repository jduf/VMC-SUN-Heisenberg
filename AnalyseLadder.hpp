#ifndef DEF_ANALYSELADDER
#define DEF_ANALYSELADDER

#include "Analyse.hpp"

class AnalyseLadder : public Analyse{
	public:
		AnalyseLadder(std::string const& sim, unsigned int const& max_level, unsigned int const& bash_file);
		/*!Default destructor*/
		~AnalyseLadder();
		/*{Forbidden*/
		AnalyseLadder() = delete;
		AnalyseLadder(AnalyseLadder const&) = delete;
		AnalyseLadder(AnalyseLadder&&) = delete;
		AnalyseLadder& operator=(AnalyseLadder) = delete;
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
