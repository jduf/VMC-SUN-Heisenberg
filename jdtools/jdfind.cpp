#include "Directory.hpp"
#include "Read.hpp"
#include "Write.hpp"
#include "RST.hpp"
#include "Array2D.hpp"

#include <string>
#include <vector>
#include <cstdlib> // system(std::string commad)

void create_tag_file_list(std::string h, std::vector<std::string>& df){
	size_t t0(0);
	size_t t1(0);
	std::string s("");
	while(t0 != std::string::npos){
		t0 = h.find(':',t0);
		t1 = h.find(':',t0+1);
		if(h[t0-1] == '\n' ){
			//std::cout<<h.substr(t0 , t1-t0+1)<< f<<std::endl;
			//w<<h.substr(t0, t1-t0+1)<<f<<Write::endl;
			s+=h.substr(t0+1, t1-t0-1);
		}
		t0 = t1;
	}
	df.push_back(s);
}

int main(){
	Directory d;
	d.search_ext(".jdbin","/home/jdufour/travail/SU4HL");
	std::string save_in("/home/jdufour/travail/cpp-dev/jdtools/rst-output/");
	std::string rsthtml_file("");
	std::string data_file("");
	std::string file("");
	std::string command("");
	std::vector<std::string> df;

	command = "rm "+save_in +"*";
	system(command.c_str());
	for(unsigned int i(0); i<d.size();i++){
		data_file = d.get_path(i) + "/" + d.get_name(i) + d.get_ext(i);
		Read r(data_file);
		std::string h(r.header());

		rsthtml_file = save_in + d.get_name(i) ;
		Write w(rsthtml_file + ".rst");
		w<<h<<Write::endl;
		w<<"*"+data_file+"*"<<Write::endl;

		command = "rst2html " + rsthtml_file + ".rst " + rsthtml_file + ".html";  
		system(command.c_str());
		create_tag_file_list(h,df);
		df.push_back(rsthtml_file);
	}
	Array2D<std::string> DF(df.size()/2,2);
	for(unsigned int i(0);i<df.size()/2;i++){
		for(unsigned int j(0);j<2;j++){
			DF(i,j) = df[2*i+j];
		}
	}
	Write w_tag(save_in + "TAGS.bin");
	w_tag<<DF;

	std::string file_link;
	RST rst(save_in + "index");
	rst.title("List of all the files","-");
	for(unsigned int i(0);i<d.size();i++){
		file_link = save_in + d.get_name(i) + ".html";
		file = d.get_path(i) + "/" + d.get_name(i);
		rst.hyperlink(file,file_link);
		rst.np();
	}
	command = "firefox " +save_in + "index.html &";
	system(command.c_str());
}
