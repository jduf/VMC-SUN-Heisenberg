#ifndef DEF_ANALYSE
#define DEF_ANALYSE

#include "Directory.hpp"
#include "CreateSystem.hpp"

class Analyse{
	public:
		Analyse();
		virtual ~Analyse(){};
		void go(std::string search_in);

	protected:
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
		IOFiles* data_write_;

		unsigned int level_;

		Vector<unsigned int> ref_;
		unsigned int N_;
		unsigned int m_;
		unsigned int n_;
		Vector<unsigned int> M_;
		int bc_;
		Data<double> E_;
		DataSet<double> corr_;
		DataSet<double> long_range_corr_;

		void recursive_search();
		void search_jdbin();
		void extract_jdbin();

		virtual void open_files(std::string const& jdfile, std::string const& datafile, Directory const& d)=0;
		virtual void close_files()=0;
		virtual void extract_level_5()=0;
		virtual void extract_level_4();
		virtual void extract_level_3();
		virtual void extract_level_2();
		virtual void extract_level_1();
};
#endif
