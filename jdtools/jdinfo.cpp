/*!  @file jdinfo.cpp */

#include "Directory.hpp"
#include "Parseur.hpp"
#include "Read.hpp"
#include "RSTFile.hpp"

void update_readme(Directory const& d, std::string const& directory_name);
void create_rst_plot(unsigned int i, Directory const& d, std::string const& directory_name);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	std::string directory_name(P.get<std::string>("0"));
	if(!P.status()){
		Linux command;
		std::string save_in("");

		if(directory_name == "."){ directory_name = command.pwd(); }
		else { directory_name = command.pwd()+directory_name; }

		if(directory_name[directory_name.size()-1] != '/'){ directory_name += "/"; }
		save_in = directory_name + "info/";

		Directory d;
		d.search_file_ext(".png",directory_name,false,true);
		d.sort();

		update_readme(d,directory_name);

		command("firefox " + save_in + "README.html &");
	} else {
		std::cerr<<"need to give a directory"<<std::endl;
	}
}

void update_readme(Directory const& d, std::string const& directory_name){
	Read r(directory_name + "README");
	std::string h("");
	r>>h;

	RSTFile rst(directory_name + "info/","README");
	rst.text(h);
	rst.hyperlink("List of all simulations", directory_name + "info/index.html");
	rst.np();
	rst.title("PLOTS","=");
	if(d.size()<10){
		for(unsigned int i(0);i<d.size();i++){
			rst.figure(d[i],d[i],80);
		}
	} else {
		for(unsigned int i(0);i<d.size();i++){
			create_rst_plot(i,d,directory_name + "info/");
			rst.hyperlink(d.get_name(i), directory_name + "info/" + d.get_name(i) + ".html");
		}
	}
	rst.pdf();
}

void create_rst_plot(unsigned int i, Directory const& d, std::string const& directory_name){
	RSTFile rst(directory_name, d.get_name(i));
	rst.title(d.get_name(i),"=");
	if(i==0){ rst.hyperlink(d.get_name(i+1), directory_name + d.get_name(i+1) + ".html");}
	else if(i+1==d.size()){ rst.hyperlink(d.get_name(i-1), directory_name + d.get_name(i-1) + ".html"); }
	else { rst.hyperlink(d.get_name(i-1), directory_name + d.get_name(i-1) + ".html"); rst.hyperlink(d.get_name(i+1), directory_name + d.get_name(i+1) + ".html"); }
	rst.np();
	rst.hyperlink(directory_name,directory_name + "README.html");
	rst.figure(d[i],d[i],80);
}
