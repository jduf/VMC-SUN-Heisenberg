#ifndef DEF_ANALYSESYSTEMATICLADDER
#define DEF_ANALYSESYSTEMATICLADDER

#include "Analyse.hpp"
#include "VMCSystematic.hpp"

class AnalyseSystematicLadder : public Analyse{
	public:
		AnalyseSystematicLadder(std::string const& sim, std::string const& path, unsigned int const& max_level, unsigned int const& run_cmd, bool const& display_results);
		/*!Default destructor*/
		~AnalyseSystematicLadder() = default;
		/*{Forbidden*/
		AnalyseSystematicLadder() = delete;
		AnalyseSystematicLadder(AnalyseSystematicLadder const&) = delete;
		AnalyseSystematicLadder(AnalyseSystematicLadder&&) = delete;
		AnalyseSystematicLadder& operator=(AnalyseSystematicLadder) = delete;
		/*}*/

	protected:
		bool display_results_;
		List<MCSim> kept_samples_all_;
		List<MCSim> kept_samples_J_;
		List<MCSim> best_J_;
		List<MCSim> best_wf_;
		List<std::string> best_sim_;
		std::string date_;

		void open_files();
		void close_files();

		std::string extract_level_9();
};
#endif
