#include "Directory.hpp"

Directory::Directory():
	dir("."),
	fname(0),
	path(0),
	ext(0)
{ }


void Directory::search_ext(std::string extension, std::string curr_dir){
	DIR* dir_point = opendir(curr_dir.c_str());
	dirent* entry = readdir(dir_point);
	while (entry){
		if (entry->d_type == DT_DIR){
			std::string dir = entry->d_name;
			if (dir != "." && dir != ".."){
				search_ext(extension,curr_dir+"/"+dir);
			}
		}
		else if (entry->d_type == DT_REG){
			std::string f = entry->d_name;
			if (f.find(extension, (f.size() - extension.size())) != std::string::npos){
				path.push_back(curr_dir);
				split_fname(f);
			}
		}
		entry = readdir(dir_point);
	}
	closedir(dir_point);
}

void Directory::print(){
	for(unsigned int i(0);i<path.size();i++){
		std::cout<<path[i]<<"/"<<fname[i]<<ext[i]<<std::endl;
	}
}

void Directory::split_fname(std::string f){
	unsigned int pos(f.find("."));
	fname.push_back(f.substr(0,pos));
	ext.push_back(f.substr(pos));
}
