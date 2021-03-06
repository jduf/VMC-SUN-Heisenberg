#ifndef DEF_ANALYSEMIN
#define DEF_ANALYSEMIN

#include "Analyse.hpp"
#include "VMCAnalyse.hpp"

class AnalyseMin : public Analyse{
	public:
		AnalyseMin(std::string const& sim, unsigned int const& max_level, unsigned int const& bash_file, unsigned int const& ref);
		~AnalyseMin() = default;

	protected:
		void open_files();
		void close_files();

		std::string extract_level_9();
		std::string extract_level_8();
		std::string extract_level_7();
		std::string extract_level_6();
		std::string extract_level_5();
		std::string extract_level_4();
		std::string extract_level_3();
		std::string extract_level_2();
};
#endif
