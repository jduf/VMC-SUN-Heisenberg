#ifndef DEF_DIRECTORY
#define DEF_DIRECTORY

#include<dirent.h>
#include<iostream>
#include<vector>

class Directory {
	public:
		/*!Constructor*/
		Directory(){}
		/*!Destructor*/
		~Directory(){}
		/*!Removes all the path_+filename_+extention_ stored*/
		void set();

		/*!Returns the ith filename_*/
		inline std::string get_name(unsigned int i) const {return filename_[i];};
		/*!Returns the ith path_*/
		inline std::string get_path(unsigned int i) const {return path_[i];};
		/*!Returns the ith extension_*/
		inline std::string get_ext(unsigned int i) const {return ext_[i];};
		/*!Returns the ith path_+filename_+extention_*/
		inline std::string operator[](unsigned int i) const {return path_[i]+filename_[i]+ext_[i];};
		/*!Returns the number of files*/
		inline unsigned int size() const {return path_.size();};

		/*!Prints all the path_+filename_+extention_*/
		void print();
		/*!Lists all directories stored in curr_dir*/
		void list_dir(std::string curr_dir);
		/*!Lists all files matching part ok keyword*/
		void search_file(std::string const& keyword, std::string curr_dir, bool follow_link, bool recursive);
		/*!Lists all files whose enxtention is*/
		void search_file_ext(std::string const& ext, std::string curr_dir, bool follow_link, bool recursive);
		/*!Sorts by alphabetical order all the path_+filename_+extention_*/
		void sort();

	private:
		/*!Forbids copy*/
		Directory(Directory const& d);
		/*!Forbids assignment*/
		Directory& operator=(Directory d);

		std::vector<std::string> path_;		//!< stores the paths
		std::vector<std::string> filename_;	//!< stores the filenames
		std::vector<std::string> ext_;		//!< stores the extentions

		/*!Splits a filname into filename_+ext_*/
		void split_ext(std::string f);
};
#endif
