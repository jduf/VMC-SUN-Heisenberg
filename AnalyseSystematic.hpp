#ifndef DEF_ANALYSESYSTEMATIC
#define DEF_ANALYSESYSTEMATIC

#include "Analyse.hpp"
#include "VMCSystematic.hpp"

class AnalyseSystematic : public Analyse{
	public:
		AnalyseSystematic(std::string const& sim, unsigned int const& max_level, unsigned int const& bash_file);
		/*!Default destructor*/
		~AnalyseSystematic();
		/*{Forbidden*/
		AnalyseSystematic() = delete;
		AnalyseSystematic(AnalyseSystematic const&) = delete;
		AnalyseSystematic(AnalyseSystematic&&) = delete;
		AnalyseSystematic& operator=(AnalyseSystematic) = delete;
		/*}*/

	protected:
		List<MCSim> kept_samples_;

		void open_files();
		void close_files();

		std::string extract_level_9();
		std::string extract_level_2();
};
#endif
