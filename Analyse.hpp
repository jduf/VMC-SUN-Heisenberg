#ifndef DEF_ANALYSE
#define DEF_ANALYSE

#include "Directory.hpp"
#include "CreateSystem.hpp"

class Analyse: public IOSystem{
	public:
		Analyse(std::string const& sim, std::string const& path, unsigned int const& max_level, unsigned int const& run_cmd);
		/*Default destructor*/
		virtual ~Analyse();
		/*{Forbidden*/
		Analyse() = delete;
		Analyse(Analyse const&) = delete;
		Analyse(Analyse&&) = delete;
		Analyse& operator=(Analyse const&) = delete;
		/*}*/

	protected:
		List<std::string> all_link_names_;
		List<std::string> all_link_files_;
		List<RSTFile> list_rst_;

		unsigned int const max_level_;
		unsigned int nof_ 		  = 0;
		unsigned int level_ 	  = 0;
		std::string rel_level_ 	  = "";
		bool consider_only_most_recent_jdbin_ = false;

		void do_analyse();
		/*!At each level, this method is called to create usefull output files*/
		virtual void open_files()  = 0;
		/*!At each level, this method is called to close open files*/
		virtual void close_files() = 0;

		std::string extract_default();
		std::string extract_best_of_previous_level();
		std::string fit_therm_limit();

	private:
		unsigned int study_;
		unsigned int const run_cmd_;
		
		void recursive_search();
		void search_jdbin();
		void extract_jdbin();
};
#endif
