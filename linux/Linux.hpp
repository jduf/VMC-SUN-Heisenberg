#ifndef DEF_LINUX
#define DEF_LINUX

/*{Description*/
/*!The following macro are essential. They need to be set to the correct
 * path/executable in order to run.*/
/*}*/
#define MY_RST2HTML_STYLESHEET "/home/jdufour/travail/cpp-dev/rst/css/best.css"
#define MY_BIN_RST2HTML "rst2html"
#define MY_BIN_RST2LATEX "rst2latex"
#define MY_BIN_PDFLATEX "pdflatex"
#define MY_BIN_GNUPLOT "gnuplot"
#define MY_BIN_LATEX "latex"
#define MY_BIN_DVIPDF "dvipdf"
#define MY_BIN_PDFCROP "pdfcrop"
#define MY_BIN_PDF2PNG "convert"
#define MY_BIN_HTMLBROWSER "/usr/bin/firefox"

#include <cstdlib>
#include <string>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include "sys/stat.h"
#include "string.h"

class Linux{
	public:
		/*!Default constructor*/
		Linux() = default;
		/*!Default destructor*/
		~Linux() = default;
		/*{Forbidden*/
		Linux(Linux const&) = delete;
		Linux(Linux&&) = delete;
		Linux& operator=(Linux) = delete;
		/*}*/

		/*!Execute a UNIX command and get its exit value*/
		void operator()(std::string cmd, bool silent);
		static void open(std::string const& filename);
		static void close(bool const& run_now);

		/*!Returns exit value of the last command*/
		int status(){ return ev_; }
		/*!Returns a string containing the current path*/
		std::string pwd(){ return std::string(get_current_dir_name()) + '/'; }

		/*!Creates a directory*/
		void mkdir(const char *directory, mode_t mode = 0755);
		/*!Creates a path*/
		void mkpath(const char *path, mode_t mode = 0755);

		static std::string latex(std::string const& path, std::string const& filename);
		static std::string pdflatex(std::string const& path, std::string const& filename);
		static std::string dvipdf(std::string const& path, std::string const& filename);
		static std::string pdfcrop(std::string const& path, std::string const& filename);
		static std::string pdf2png(std::string const& infile, std::string const& outfile);
		/*{*//*!Creates gnuplot plots
			   Using a simple gnuplot file (with extension .gp) creates and
			   .eps picture and .tex file which can be used to create .pdf
			   files via Linux::pdflatex
			   A default size size is set such that its ratio equals the
			   golden number and when reduced by 70%, fits perfectly in one
			   column in revtex-4.1 articles
			   To specify a personal ratio, the first line of the .gp file
			   must start with "#latex_size x,y" where x,y are dimension of
			   the desired picture *//*}*/
		static std::string gp2latex(std::string const& texfile, std::string const& path, std::string const& gpfile);
		static std::string rst2latex(std::string const& texfile, std::string const& path, std::string const& filename);
		static std::string rst2html(std::string const& path, std::string const& filename);
		static std::string html_browser(std::string const& html);

	private:
		class Bash{
			public:
				Bash() = default;
				~Bash();
				std::ofstream file_;
				std::string filename_;
				/*{Forbidden*/
				Bash(Bash const&) = delete;
				Bash(Bash&&) = delete;
				Bash& operator=(Bash) = delete;
				/*}*/
		};

		int ev_ = 0;//!< exit value of the last UNIX command
		static Bash bash_;
};
#endif
