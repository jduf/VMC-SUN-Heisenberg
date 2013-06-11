/*!  @file jdfind.cpp */

#include "Directory.hpp"
#include "Parseur.hpp"
#include "Read.hpp"
#include "Write.hpp"
#include "RSTfile.hpp"
#include "Array2D.hpp"

#include <string>
#include <vector>
#include <cstdlib> // system(std::string commad)


void create_readme(std::string const& directory_name);
void list_all_simulation_files(Directory const& d, std::string const& save_in);
void create_all_simulation_files(Directory const& d, std::string const& save_in);
void create_tag_list(Directory const& d, std::string const& save_in);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	if(argc == 2){
		char buff[PATH_MAX];
		getcwd(buff,PATH_MAX);
		std::string directory_name("");
		std::string save_in("");
		std::string command("");

		P.set(directory_name);
		if(directory_name == "."){
			directory_name = buff;
			if(directory_name[directory_name.size()-1] != '/'){
				directory_name = directory_name + "/";
			}
			save_in = directory_name + "info/";
		} else {
			directory_name ="/"+directory_name;
			directory_name = buff+directory_name;
			if(directory_name[directory_name.size()-1] != '/'){
				directory_name = directory_name + "/";
			}
			save_in = directory_name + "info/";
		}
		Directory d;
		d.search_file_ext(".jdbin",directory_name);
		d.sort();
		command = "mkdir -p " + save_in;
		system(command.c_str());
		
		create_readme(directory_name);
		list_all_simulation_files(d,save_in);
		create_all_simulation_files(d,save_in);
		create_tag_list(d,save_in);

		command = "firefox " + save_in + "README.html &";
		system(command.c_str());
	} else {
		std::cerr<<"need to give a directory"<<std::endl;
	}
}

void create_readme(std::string const& directory_name){
	Read r(directory_name + "README");
	std::string h("");
	r>>h;

	RSTfile rst("info/README",directory_name);
	rst.text(h);
	rst.hyperlink("List of all simulations", directory_name + "info/index.html");
}

void list_all_simulation_files(Directory const& d, std::string const& save_in){
	RSTfile rst("index",save_in);
	rst.title("List of all the files","-");
	for(unsigned int i(0);i<d.size();i++){
		rst.hyperlink(d.get_name(i), save_in + d.get_name(i) + ".html");
		rst.np();
	}
}

void create_all_simulation_files(Directory const& d, std::string const& save_in){
	std::string h("");
	std::string data("");
	for(unsigned int i(0); i<d.size();i++){
		data = "";
		Read r(d[i]);
		h = r.get_header();
		Read r_data(d[i]+".dat");
		r_data>>data;

		RSTfile rst(d.get_name(i),save_in);
		rst.text(h);
		rst.textit(d[i]);
		rst.np();
		rst.title("Output","-");
		rst.lineblock(data);
	}
}

void create_tag_list(Directory const& d, std::string const& save_in){
	std::vector<std::string> df;
	std::string h("");
	for(unsigned int i(0);i<d.size();i++){
		Read r(d[i]);
		h = r.get_header();
		size_t t0(0);
		size_t t1(0);
		std::string s("");
		while(t0 != std::string::npos){
			t0 = h.find(':',t0);
			t1 = h.find(':',t0+1);
			if(h[t0-1] == '\n' ){
				s+=h.substr(t0+1, t1-t0-1);
			}
			t0 = t1;
		}
		df.push_back(s);
		df.push_back(save_in + d.get_name(i));
	}
	Array2D<std::string> DF(df.size()/2,2);
	for(unsigned int i(0);i<df.size()/2;i++){
		for(unsigned int j(0);j<2;j++){
			DF(i,j) = df[2*i+j];
		}
	}
	Write w(save_in + "TAGS.bin");
	w<<DF;
}

