#ifndef DEF_ANALYSE
#define DEF_ANALYSE

#include "Directory.hpp"
#include "CreateSystem.hpp"
#include "List.hpp"

class Analyse: public IOSystem{
	public:
		Analyse(std::string const& path);
		/*Default destructor*/
		virtual ~Analyse() = default;
		/*{Forbidden*/
		Analyse() = delete;
		Analyse(Analyse const& a) = delete;
		Analyse(Analyse&& a) = delete;
		Analyse& operator=(Analyse const&) = delete;
		/*}*/

	protected:
		List<std::string> all_link_names_;
		List<std::string> all_link_files_;
		List<RSTFile> rst_file_;

		unsigned int level_;
		unsigned int nof_;

		void do_analyse();
		virtual void open_files()=0;
		virtual void close_files()=0;
		std::string extract_level_7();

	private:
		unsigned int study_;
		
		void recursive_search();
		void search_jdbin();
		void extract_jdbin();
};
#endif
