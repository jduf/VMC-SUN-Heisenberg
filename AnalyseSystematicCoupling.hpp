#ifndef DEF_ANALYSESYSTEMATICCOUPLING
#define DEF_ANALYSESYSTEMATICCOUPLING

#include "Analyse.hpp"
#include "VMCSystematic.hpp"

class AnalyseSystematicCoupling : public Analyse{
	public:
		AnalyseSystematicCoupling(std::string const& sim, unsigned int const& max_level, unsigned int const& run_cmd, bool const& display_results, unsigned int const& ref);
		/*!Default destructor*/
		~AnalyseSystematicCoupling() = default;
		/*{Forbidden*/
		AnalyseSystematicCoupling() = delete;
		AnalyseSystematicCoupling(AnalyseSystematicCoupling const&) = delete;
		AnalyseSystematicCoupling(AnalyseSystematicCoupling&&) = delete;
		AnalyseSystematicCoupling& operator=(AnalyseSystematicCoupling) = delete;
		/*}*/

	protected:
		bool display_results_;
		List<MCSim> kept_samples_all_;
		List<MCSim> kept_samples_J_;
		List<MCSim> best_J_;

		void open_files();
		void close_files();

		std::string extract_level_9();
};
#endif
