#define MY_RST2HTML_STYLESHEET "/home/jdufour/travail/cpp-dev/rst/css/best.css "
#define MY_BIN_RST2HTML "rst2html "
#define MY_BIN_RST2LATEX "rst2latex "
#define MY_BIN_PDFLATEX "pdflatex "
#define MY_BIN_GNUPLOT "gnuplot "
#define MY_BIN_LATEX "latex "
#define MY_BIN_DVIPDF "dvipdf "
#define MY_BIN_PDFCROP "pdfcrop "
#define MY_BIN_PDF2PNG "convert "

#ifndef DEF_LINUX
#define DEF_LINUX

#include <cstdlib> 
#include <string>
#include <unistd.h>

class Linux {
	public:
		/*!Constructor*/
		Linux(){}
		/*!Destructor*/
		~Linux(){}

		/*!Execute a UNIX command and get its exit value*/
		void operator()(std::string cmd, bool silent=false){ ev_=system((cmd+(silent?"> /dev/null 2> /dev/null":"")).c_str()); }
		/*!Returns exit value of the last command*/
		int status(){return ev_;}
		/*!Returns a string containing the current path*/
		std::string pwd(){ return std::string(get_current_dir_name()) + '/'; }
		
		static std::string rst2html(std::string const& path, std::string const& filename){
			std::string cmd(MY_BIN_RST2HTML);
			cmd+= "--stylesheet=" + std::string(MY_RST2HTML_STYLESHEET);
			cmd+= "--field-name-limit=0 "; 
			cmd+= path+filename+".rst ";
			cmd+= path+filename+".html ";
			return cmd;
		}

		static std::string rst2latex(std::string const& path, std::string const& filename){
			std::string cmd(MY_BIN_RST2LATEX);
			cmd+= path+filename+".rst ";
			cmd+= path+filename+".tex ";
			return cmd;
		}

		static std::string pdflatex(std::string const& path, std::string const& filename){
			std::string cmd(MY_BIN_PDFLATEX);
			cmd+= "-shell-escape -output-directory " + path + " "  + filename + ".tex";
			return cmd;
		}

		static std::string latex(std::string const& filename){
			std::string cmd(MY_BIN_LATEX);
			cmd+= "-shell-escape ";
			cmd+= filename + ".tex ";
			return cmd;
		}

		static std::string dvipdf(std::string const& path, std::string const& filename){
			std::string cmd(MY_BIN_DVIPDF);
			cmd+= filename + ".dvi ";
			cmd+= path+filename + ".pdf ";
			return cmd;
		}

		static std::string pdfcrop(std::string const& path, std::string const& filename){
			std::string cmd(MY_BIN_PDFCROP);
			cmd+= path + filename + ".pdf ";
			cmd+= path + filename + ".pdf > /dev/null"; 
			return cmd;
		}

		static std::string gp2latex(std::string const& texfile, std::string const& path, std::string const& gpfile){
			std::string cmd(MY_BIN_GNUPLOT);
			cmd = "cd " + path + "; " + cmd;
			cmd+= "-e \"set terminal epslatex color size 12.5cm,7.73 standalone lw 2 header \'\\\\usepackage{amsmath,amssymb}\'; set output \'" + texfile + ".tex\'\" ";
			cmd+= path + gpfile + ".gp";
			return cmd;
		}

		static std::string pdf2png(std::string const& infile, std::string const& outfile){
			std::string cmd(MY_BIN_PDF2PNG);
			cmd+= "-density 500 -resize 20% " + infile + ".pdf " + outfile + ".png";
			return cmd;
		}

	private:
		/*!Forbids copy*/
		Linux(Linux const& l);
		/*!Forbids assignment*/
		Linux& operator=(Linux const& l);

		int ev_;//!< exit value of the last UNIX command
};
#endif
