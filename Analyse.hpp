#ifndef DEF_ANALYSE
#define DEF_ANALYSE

#include "Directory.hpp"
#include "CreateSystem.hpp"
#include "List.hpp"

class Analyse: public IOSystem{
	public:
		Analyse(std::string const& path, unsigned int const& max_level);
		/*Default destructor*/
		virtual ~Analyse() = default;
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

		std::string rel_level_;
		unsigned int const max_level_;
		unsigned int level_;
		unsigned int nof_;

		void do_analyse();
		virtual void open_files()  = 0;
		virtual void close_files() = 0;

	private:
		unsigned int study_;
		
		void recursive_search();
		void search_jdbin();
		void extract_jdbin();
};
#endif
