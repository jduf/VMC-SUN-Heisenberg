#ifndef DEF_ANALYSESYSTEMATIC
#define DEF_ANALYSESYSTEMATIC

#include "Analyse.hpp"
#include "VMCSystematic.hpp"

class AnalyseSystematic : public Analyse{
	public:
		AnalyseSystematic(std::string const& sim, std::string const& path, unsigned int const& max_level, unsigned int const& run_cmd, bool const& display_results);
		/*!Default destructor*/
		~AnalyseSystematic();
		/*{Forbidden*/
		AnalyseSystematic() = delete;
		AnalyseSystematic(AnalyseSystematic const&) = delete;
		AnalyseSystematic(AnalyseSystematic&&) = delete;
		AnalyseSystematic& operator=(AnalyseSystematic) = delete;
		/*}*/

	protected:
		bool display_results_;
		List<MCSim> kept_samples_;

		void open_files();
		void close_files();

		std::string extract_level_9();
		std::string extract_level_2();
};
#endif
