#ifndef DEF_ANALYSEEXTRACT
#define DEF_ANALYSEEXTRACT

#include "Analyse.hpp"
#include "VMCExtract.hpp"

/*!This class must be used in a tree that contains different types of phase
 * space in level 9 (could contain only one type of phase space but should not
 * handle the same phase space at different stage of the optimisation because
 * it can't be used with consider_only_most_recent_jdbin_ = true;) */
class AnalyseExtract : public Analyse{
	public:
		AnalyseExtract(std::string const& sim, unsigned int const& max_level, unsigned int const& bash_file, unsigned int const& display_results, unsigned int const& ref);
		~AnalyseExtract();
		/*{Forbidden*/
		AnalyseExtract() = delete;
		AnalyseExtract(AnalyseExtract const&) = delete;
		AnalyseExtract(AnalyseExtract&&) = delete;
		AnalyseExtract& operator=(AnalyseExtract) = delete;
		/*}*/

	protected:
		unsigned int const display_results_;
		List<MCSim> kept_samples_;

		void open_files();
		void close_files();

		std::string extract_level_9();
		std::string extract_level_8();
		std::string extract_level_7();
		std::string extract_level_6();
		std::string extract_level_4();
		std::string extract_level_3();
		std::string extract_level_2();
};
#endif
