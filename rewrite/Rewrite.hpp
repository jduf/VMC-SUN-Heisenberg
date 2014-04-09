#ifndef DEF_REWRITE
#define DEF_REWRITE

#include "Header.hpp"

#include<iostream>
#include<fstream>
#include<string>
#include<stdio.h>

class Rewrite{
	public:
		/*!Opens a file named "filename", by default open a binary file*/
		Rewrite(std::string filename);
		/*!Closes the file*/
		~Rewrite();

		void rewrite_header(std::string s);

	private:
		/*!Forbids copy constructor*/
		Rewrite(Rewrite const& s);
		/*!Forbids assertion operator*/
		Rewrite& operator=(Rewrite const&);

		std::string filename; //!< name of the file to Rewrite in
		FILE *bfile; //!< pointer on the binery file to Rewrite in
		bool unlocked; //!< true if the file is ready to be Rewrite in
		size_t reading_point_; //!< last bit read
};
#endif
