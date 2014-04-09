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
		inline std::string operator[](unsigned int i) const {return path[i]+fname[i]+ext[i];};
		inline unsigned int size() const {return path.size();};

		void print();
		void list_dir(std::string curr_dir);
		void search_file(std::string const& keyword, std::string curr_dir, bool follow_link, bool recursive);
		void search_file_ext(std::string const& ext, std::string curr_dir, bool follow_link, bool recursive);
		void sort();

		void clear();

	private:
		std::string dir;
		std::vector<std::string> fname;
		std::vector<std::string> path;
		std::vector<std::string> ext;

		void split_ext(std::string f);
};
#endif
