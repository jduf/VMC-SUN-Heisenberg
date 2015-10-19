#ifndef DEF_ANALYSEHONEYCOMB
#define DEF_ANALYSEHONEYCOMB

#include "Analyse.hpp"

class AnalyseHoneycomb : public Analyse{
	public:
		AnalyseHoneycomb(std::string const& sim, unsigned int const& max_level);
		/*!Default destructor*/
		~AnalyseHoneycomb() = default;
		/*{Forbidden*/
		AnalyseHoneycomb() = delete;
		AnalyseHoneycomb(AnalyseHoneycomb const&) = delete;
		AnalyseHoneycomb(AnalyseHoneycomb&&) = delete;
		AnalyseHoneycomb& operator=(AnalyseHoneycomb) = delete;
		/*}*/

	protected:
		void open_files();
		void close_files();

		std::string extract_level_6();
};
#endif
