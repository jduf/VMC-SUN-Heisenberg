#ifndef DEF_LINUX
#define DEF_LINUX

/*{Description*/
/*!The following macro a essential. They need to be set to the correct
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
#define MY_BIN_MKDIR "/bin/mkdir"

#include <cstdlib> 
#include <string>
#include <unistd.h>
#include <fstream>
#include <iostream>

class Linux {
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
		void operator()(std::string cmd, bool silent){ ev_=system((cmd+(silent?" > /dev/null 2> /dev/null":"")).c_str()); }
		/*!Returns exit value of the last command*/
		int status(){return ev_;}
		/*!Returns a string containing the current path*/
		std::string pwd(){ return std::string(get_current_dir_name()) + '/'; }

		/*!Creates a directory with -p option*/
		void mkdir(std::string const& directory){ 
			std::string cmd(MY_BIN_MKDIR);
			(*this)(cmd + " -p " + directory, true); 
		}

		static std::string latex(std::string const& path, std::string const& filename){
			std::string cmd(MY_BIN_LATEX);
			cmd+= " -output-directory " +path + " ";
			cmd+= path+filename + ".tex";
			return cmd;
		}

		static std::string pdflatex(std::string const& path, std::string const& filename){
			std::string cmd(MY_BIN_PDFLATEX);
			cmd+= " -shell-escape";
			cmd+= " -output-directory " +path + " ";
			cmd+= filename + ".tex";
			return cmd;
		}

		static std::string dvipdf(std::string const& path, std::string const& filename){
			std::string cmd(MY_BIN_DVIPDF);
			cmd+= " " + path  + filename + ".dvi ";
			cmd+= path + filename + ".pdf ";
			return cmd;
		}

		static std::string pdfcrop(std::string const& path, std::string const& filename){
			std::string cmd(MY_BIN_PDFCROP);
			cmd+= " "  + path + filename + ".pdf ";
			cmd+=       path + filename + ".pdf > /dev/null"; 
			return cmd;
		}

		static std::string pdf2png(std::string const& infile, std::string const& outfile){
			std::string cmd(MY_BIN_PDF2PNG);
			cmd+= " -density 500 -resize 20% " + infile + ".pdf " + outfile + ".png";
			return cmd;
		}

		/*{Description*/
		/*!Using a simple gnuplot file (with extension .gp) creates and .eps
		 * picture and .tex file which can be used to create .pdf files via
		 * Linux::pdflatex
		 * A default size size is set such that its ratio equals the golden
		 * number and when reduced by 70%, fits perfectly in one column in
		 * revtex-4.1 articles
		 * To specify a personal ratio, the first line of the .gp file must
		 * start with "#latex_size x,y" where x,y are dimension of the desired
		 * picture
		 * */
		/*}*/
		static std::string gp2latex(std::string const& texfile, std::string const& path, std::string const& gpfile){
			std::string cmd(MY_BIN_GNUPLOT);
			cmd = "cd " + path + "; " + cmd;
			std::ifstream file(path+gpfile+".gp",std::ifstream::in);
			std::string size;
			if(file.is_open() && std::getline(file,size) && size.find("#latex_size") != std::string::npos){ 
				size = size.substr(12); 
				std::cerr<<__PRETTY_FUNCTION__<<" : set size "<<size<<std::endl;
			} else { size = "12.15cm,7.54"; } 
			cmd+= " -e \"set terminal epslatex color size "+size+" standalone lw 2 header \'\\\\usepackage{amsmath,amssymb}\'; set output \'" + texfile + ".tex\'\" ";
			cmd+= path + gpfile + ".gp";
			return cmd;
		}

		static std::string rst2latex(std::string const& path, std::string const& filename){
			std::string cmd(MY_BIN_RST2LATEX);
			cmd+= " " + path + filename + ".rst ";
			cmd+=       path + filename + ".tex ";
			return cmd;
		}

		static std::string rst2html(std::string const& path, std::string const& filename){
			std::string cmd(MY_BIN_RST2HTML);
			cmd+= " --stylesheet=" + std::string(MY_RST2HTML_STYLESHEET);
			cmd+= " --field-name-limit=0 "; 
			cmd+= path + filename + ".rst ";
			cmd+= path + filename + ".html ";
			return cmd;
		}

	private:
		int ev_ = 0;//!< exit value of the last UNIX command
};
#endif
