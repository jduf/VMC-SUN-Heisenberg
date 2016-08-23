/*! @file rst2pdf.cpp*/

#include "Linux.hpp"
#include "Miscellaneous.hpp"

int main(int argc,char* argv[]){
	const std::string ext(".rst");
	if(argc==2){
		std::string filename(argv[1]);
		if ( filename != ext && filename.size() > ext.size() && filename.substr(filename.size() - ext.size()) == ext ) {
			filename = filename.substr(0, filename.size() - ext.size());
			std::string path("");

			std::vector<std::string> tmp(my::string_split(path+filename,'/'));
			for(unsigned int i(0);i<tmp.size()-1;i++){ path += tmp[i] + "/"; }
			filename = tmp[tmp.size()-1];

			Linux command;
			if(path[0] != '/'){ path = command.pwd() + path; }
			std::string texfile("/tmp/rst2latex_tmp");
			command(Linux::rst2latex(texfile,path,filename),false);
			if(!command.status()){
				command("sed -i 's/scale=2\\.000000/width=1\\\\linewidth/g' "+texfile+".tex",false);
				command(Linux::pdflatex("/tmp/",texfile),true);
				command("mv " + texfile + ".pdf " + path + filename + ".pdf",false);
				command("rm " + texfile + "*",false);
			} else { std::cerr<<__PRETTY_FUNCTION__<<" : the script contains mistakes"<<std::endl; }
		} else { std::cerr<<__PRETTY_FUNCTION__<<" : the filename must have a '"<<ext<<"' extension"<<std::endl; }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : a reStructrueText file (*"<<ext<<") must be given as argument "<<std::endl; }
}
