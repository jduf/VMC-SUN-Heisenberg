#ifndef DEF_ANALYSELADDER
#define DEF_ANALYSELADDER

#include "Analyse.hpp"
#include "VMCMinimization.hpp"

class AnalyseLadder : public Analyse{
	public:
		AnalyseLadder(std::string const& path, unsigned int const& max_level);
		~AnalyseLadder() = default;

	protected:
		void open_files();
		void close_files();

		std::string extract_level_9();
		std::string extract_level_8();
		std::string extract_level_7();
		std::string extract_level_6();
		std::string extract_level_5();
};
#endif
