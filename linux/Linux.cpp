#include "Linux.hpp"

void Linux::operator()(std::string cmd, bool silent){
	cmd += silent?" > /dev/null 2>&1":"";
	if(Linux::bash_.file_.is_open()){ Linux::bash_.file_<<cmd<<std::endl; }
	else { ev_=system(cmd.c_str()); }
}

Linux::Bash::~Bash(){
	if(file_.is_open()){ file_.close(); }
	std::ofstream file_;
}

Linux::Bash Linux::bash_;

void Linux::open(std::string const& filename){
	std::cerr<<__PRETTY_FUNCTION__<<" : open Bash file ("<<filename<<")"<<std::endl;
	if(bash_.file_.is_open()){
		std::cerr<<__PRETTY_FUNCTION__<<" : can't open more than one Bash file at the same time"<<std::endl;
	} else {
		bash_.file_.open(filename.c_str(), std::ios::out);
		if(bash_.file_.is_open()){
			bash_.filename_ = filename;
			bash_.file_<<"#!/bin/bash\n";
		} else { std::cerr<<__PRETTY_FUNCTION__<<" : the bash file can't be openend"<<std::endl; }
	}
}

void Linux::close(bool const& run_now){
	if(bash_.file_.is_open()){
		std::cerr<<__PRETTY_FUNCTION__<<" : close Bash file ("<<bash_.filename_<<")";
		bash_.file_.close();
		if(run_now){
			Linux command;
			command("chmod 755 "+Linux::bash_.filename_,false);
			command("./"+Linux::bash_.filename_+" &",false);
			bash_.filename_ = "";
			std::cerr<<" and run commands";
		}
		std::cerr<<std::endl;
	}
}

/*{methods that don't need system*/
void Linux::mkdir(const char *directory, mode_t mode){
	struct stat st;
	ev_ = 0;
	if (stat(directory, &st) != 0){
		if(::mkdir(directory, mode) != 0 && errno != EEXIST){
			ev_ = 1;
		}
	} else if(!S_ISDIR(st.st_mode)){
		errno = ENOTDIR;
		ev_ = 1;
	}
}

void Linux::mkpath(const char *path, mode_t mode){
	char *pp;
	char *sp;
	char *copypath = strdup(path);

	pp = copypath;
	while (ev_ == 0 && (sp = strchr(pp, '/')) != 0){
		if(sp != pp){
			*sp = '\0';
			mkdir(copypath, mode);
			*sp = '/';
		}
		pp = sp + 1;
	}
	if(ev_ == 0) { mkdir(path, mode); }
	free(copypath);
}
/*}*/

/*{methods returning a std::string that is then to by executed by system*/
std::string Linux::latex(std::string const& path, std::string const& filename){
	std::string cmd(MY_BIN_LATEX);
	cmd+= " -output-directory " +path + " ";
	cmd+= path+filename + ".tex";
	return cmd;
}

std::string Linux::pdflatex(std::string const& path, std::string const& filename){
	std::string cmd(MY_BIN_PDFLATEX);
	cmd+= " -shell-escape";
	cmd+= " -output-directory " +path + " ";
	cmd+= filename + ".tex";
	return cmd;
}

std::string Linux::dvipdf(std::string const& path, std::string const& filename){
	std::string cmd(MY_BIN_DVIPDF);
	cmd+= " " + path  + filename + ".dvi ";
	cmd+= path + filename + ".pdf ";
	return cmd;
}

std::string Linux::pdfcrop(std::string const& path, std::string const& filename){
	std::string cmd(MY_BIN_PDFCROP);
	cmd+= " "  + path + filename + ".pdf ";
	cmd+=       path + filename + ".pdf > /dev/null";
	return cmd;
}

std::string Linux::pdf2png(std::string const& infile, std::string const& outfile){
	std::string cmd(MY_BIN_PDF2PNG);
	cmd+= " -density 500 -resize 20% " + infile + ".pdf " + outfile + ".png";
	return cmd;
}

std::string Linux::gp2latex(std::string const& texfile, std::string const& path, std::string const& gpfile){
	std::string cmd(MY_BIN_GNUPLOT);
	cmd = "cd " + path + "; " + cmd;
	std::ifstream file(path+gpfile+".gp",std::ifstream::in);
	std::string size;
	if(file.is_open() && std::getline(file,size) && size.find("#latex_size") != std::string::npos){
		size = size.substr(12);
		std::cerr<<__PRETTY_FUNCTION__<<" : set size "<<size<<std::endl;
	} else { size = "12.15cm,7.54"; }
	cmd+= " -e \"set terminal epslatex color size "+size+" standalone lw 2 header \'\\\\usepackage{amsmath,amssymb}\'; set output \'" + texfile + ".tex\'\" ";
	cmd+= gpfile + ".gp";
	cmd+= "; cd - > /dev/null";
	return cmd;
}

std::string Linux::rst2latex(std::string const& texfile, std::string const& path, std::string const& filename){
	std::string cmd(MY_BIN_RST2LATEX);
	cmd+= " --latex-preamble='\\usepackage{grffile}\\usepackage[a4paper,total={13cm,27cm}]{geometry}'";
	cmd+= " " + path + filename + ".rst ";
	cmd+=       texfile + ".tex ";
	return cmd;
}

std::string Linux::rst2html(std::string const& path, std::string const& filename){
	std::string cmd(MY_BIN_RST2HTML);
	cmd+= " -st";
	cmd+= " --stylesheet=" + std::string(MY_RST2HTML_STYLESHEET);
	cmd+= " --field-name-limit=0 ";
	cmd+= " --link-stylesheet ";
	cmd+= path + filename + ".rst ";
	cmd+= path + filename + ".html ";
	return cmd;
}

std::string Linux::html_browser(std::string const& html){
	std::string cmd(MY_BIN_HTMLBROWSER);
	return cmd + " " + html;
}
/*}*/
