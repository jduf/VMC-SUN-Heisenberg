#ifndef DEF_ANALYSE
#define DEF_ANALYSE

#include "Directory.hpp"
#include "List.hpp"
#include "CreateSystem.hpp"
#include "Parseur.hpp"
#include "RSTFile.hpp"

class Analyse{
	public:
		Analyse(std::string search_in);

	private:
		//GenericSystem<double>* GSR_;
		//GenericSystem<std::complex<double> >* GSC_;

		void print();

		std::string sim_;
		std::string info_;
		std::string analysis_;
		std::string path_;
		std::string dir_;
		std::string filename_;
		List<std::string> all_link_names_;
		List<std::string> all_link_files_;
		List<RSTFile> rst_;

		IOFiles* write_;

		unsigned int level_;

		void recursive_search();
		void search_jdbin();
		void extract_jdbin();
		void extract_sim();
		void extract_level_4();
		void extract_level_3();
		void extract_level_2();
		void extract_level_1();
};
#endif
