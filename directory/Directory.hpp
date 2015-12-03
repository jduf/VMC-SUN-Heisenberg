#ifndef DEF_DIRECTORY
#define DEF_DIRECTORY

#include <dirent.h>
#include <sys/stat.h>
#include <iostream>
#include <vector>

class Directory{
	public:
		/*!Default constructor*/
		Directory() = default;
		/*!Default destructor*/
		~Directory() = default;
		/*!Removes all the path_+filename_+extension_ stored*/
		void set();
		/*{Forbidden*/
		Directory(Directory const&) = delete;
		Directory(Directory&&) = delete;
		Directory& operator=(Directory) = delete;
		/*}*/

		/*!Returns the ith filename_*/
		inline std::string const& get_name(unsigned int i) const { return filename_[i]; }
		/*!Returns the ith const& path_*/
		inline std::string const& get_path(unsigned int i) const { return path_[i]; }
		/*!Returns the ith const& extension_*/
		inline std::string const& get_ext(unsigned int i) const { return ext_[i]; }
		/*!Returns the ith const& path_+filename_+extension_*/
		inline std::string operator[](unsigned int i) const { return path_[i]+filename_[i]+ext_[i]; }
		/*!Returns the number of files*/
		inline unsigned int size() const { return path_.size(); }

		/*!Prints all the path_+filename_+extension_*/
		void print() const;
		/*!Lists all directories stored in curr_dir*/
		void list_dir(std::string curr_dir);
		/*!Lists all files matching part ok keyword*/
		void search_file(std::string const& keyword, std::string curr_dir, bool follow_link, bool recursive);
		/*!Lists all files whose extension is*/
		void search_file_ext(std::string const& ext, std::string curr_dir, bool follow_link, bool recursive);
		/*!Sorts by alphabetical order all the path_+filename_+extension_*/
		void sort();

	private:
		std::vector<std::string> path_;		//!< stores the paths
		std::vector<std::string> filename_;	//!< stores the filenames
		std::vector<std::string> ext_;		//!< stores the extensions

		/*!Splits a filename into filename_+ext_*/
		void split_ext(std::string f);
};
#endif
