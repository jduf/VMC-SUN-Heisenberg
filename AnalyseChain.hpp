#ifndef DEF_ANALYSECHAIN
#define DEF_ANALYSECHAIN

#include "Analyse.hpp"

class AnalyseChain : public Analyse{
	public:
		AnalyseChain(std::string const& sim, std::string const& path, unsigned int const& max_level, unsigned int const& run_cmd);
		~AnalyseChain();

	protected:
		IOFiles* outfile_;//!< file summarizing important datas

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
