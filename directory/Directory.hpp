#ifndef DEF_DIRECTORY
#define DEF_DIRECTORY

#include<dirent.h>
#include<iostream>
#include<vector>

class Directory {
	public:
		Directory();

		inline std::string get_name(unsigned int i) const {return fname[i];};
		inline std::string get_path(unsigned int i) const {return path[i];};
		inline std::string get_ext(unsigned int i) const {return ext[i];};
		inline unsigned int size() const {return path.size();};

		void print();
		void search_ext(std::string extension, std::string curr_dir=".");

	private:
		std::string dir;
		std::vector<std::string> fname;
		std::vector<std::string> path;
		std::vector<std::string> ext;

		void split_fname(std::string f);

};

#endif
