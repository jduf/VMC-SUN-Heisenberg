#ifndef DEF_ANALYSE
#define DEF_ANALYSE

#include "Directory.hpp"
#include "CreateSystem.hpp"

class Analyse: public IOSystem{
	public:
		Analyse();
		virtual ~Analyse(){}
		void go(std::string search_in);

	protected:
		List<std::string> all_link_names_;
		List<std::string> all_link_files_;
		List<RSTFile> rst_file_;

		unsigned int level_;
		unsigned int nof_;

		void recursive_search();
		void search_jdbin();
		void extract_jdbin();

		virtual void open_files()=0;
		virtual void close_files()=0;

	private:
		/*Forbids copy*/
		Analyse(Analyse const& a);
		/*Forbids assignment*/
		Analyse& operator=(Analyse const& a);
};
#endif
