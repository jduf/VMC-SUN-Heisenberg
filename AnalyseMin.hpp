#ifndef DEF_ANALYSEMIN
#define DEF_ANALYSEMIN

#include "Analyse.hpp"
#include "VMCMinimization.hpp"

class AnalyseMin : public Analyse{
	public:
		AnalyseMin(std::string const& sim, std::string const& path, unsigned int const& max_level, bool const& run_cmd);
		~AnalyseMin() = default;

	protected:
		IOFiles complete_jobs_;

		void open_files();
		void close_files();

		std::string extract_level_9();
		std::string extract_level_8();
		std::string extract_level_7();
		std::string extract_level_6();
		std::string extract_level_5();
};
#endif
