#ifndef DEF_ANALYSELADDERPSO
#define DEF_ANALYSELADDERPSO

#include "Analyse.hpp"
#include "VMCExtract.hpp"

class AnalyseLadderPSO : public Analyse{
	public:
		AnalyseLadderPSO(std::string const& sim, unsigned int const& max_level, bool const& run_cmd, unsigned int const& display_results);
		~AnalyseLadderPSO() = default;

	protected:
		unsigned int display_results_;

		void open_files();
		void close_files();

		std::string extract_level_9();
		std::string extract_level_8();
		std::string extract_level_7();
		std::string extract_level_6();
		std::string extract_level_5();
};
#endif
