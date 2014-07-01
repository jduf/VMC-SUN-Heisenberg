#include "Directory.hpp"

Directory::Directory():
	dir("."),
	fname(0),
	path(0),
	ext(0)
{}

void Directory::search_file(std::string const& keyword, std::string curr_dir, bool follow_link, bool recursive){
	DIR* dir_point(opendir(curr_dir.c_str()));
	dirent* entry(readdir(dir_point));
	while (entry){
		if (entry->d_type == DT_DIR){
			std::string dir(entry->d_name);
			if (dir != "." && dir != ".." && recursive){
				search_file(keyword,curr_dir+dir+"/",follow_link,recursive);
			}
		}
		else if (entry->d_type == DT_REG || (follow_link && entry->d_type==DT_LNK) ){
			std::string f(entry->d_name);
			if (f.find(keyword) != std::string::npos){
				path.push_back(curr_dir);
				split_ext(f);
			}
		}
		entry = readdir(dir_point);
	}
	closedir(dir_point);
}

void Directory::search_file_ext(std::string const& extension, std::string curr_dir, bool follow_link, bool recursive){
	DIR* dir_point(opendir(curr_dir.c_str()));
	dirent* entry((readdir(dir_point)));
	while (entry){
		if (entry->d_type == DT_DIR){
			std::string dir(entry->d_name);
			if (dir != "." && dir != ".." && recursive){
				search_file_ext(extension,curr_dir+dir+"/",follow_link,recursive);
			}
		} else {
			if (entry->d_type == DT_REG || (follow_link && entry->d_type==DT_LNK) ){
				std::string f(entry->d_name);
				if (f.find(extension,f.size()-extension.size()) != std::string::npos){
					if(f.find(".") != std::string::npos){
						path.push_back(curr_dir);
						split_ext(f);
					}
				}
			}
		}
		entry = readdir(dir_point);
	}
	closedir(dir_point);
}

void Directory::list_dir(std::string curr_dir){
	DIR* dir_point = opendir(curr_dir.c_str());
	dirent* entry = readdir(dir_point);
	while (entry){
		if (entry->d_type == DT_DIR){
			std::string dir = entry->d_name;
			if (dir != "." && dir != ".."){
				path.push_back(curr_dir);
				fname.push_back(dir);
				ext.push_back("/");
			} 
		}
		entry = readdir(dir_point);
	}
	closedir(dir_point);
}

void Directory::split_ext(std::string f){
	if(f.find(".") != std::string::npos){
		size_t pos(f.find_last_of("."));
		fname.push_back(f.substr(0,pos));
		ext.push_back(f.substr(pos));
	} else {
		fname.push_back(f);
		ext.push_back("");
	}
}

void Directory::sort(){
	if(path.size()>0){
		bool sort(true);
		while(sort){
			sort=false;
			for(unsigned int i(0); i<path.size()-1;i++){
				if(path[i]+fname[i]+ext[i] > path[i+1]+fname[i+1]+ext[i+1]) { 
					sort= true;
					std::string tmp("");
					tmp = path[i];
					path[i] = path[i+1];
					path[i+1] = tmp;
					tmp =fname[i];
					fname[i] = fname[i+1];
					fname[i+1] = tmp;
					tmp = ext[i];
					ext[i] = ext[i+1];
					ext[i+1] = tmp;
				}
			}
		}
	} else {
		std::cerr<<"Directory::sort() : the file list is empty"<<std::endl;
	}
}

void Directory::clear(){
	fname.clear();
	path.clear();
	ext.clear();
}

void Directory::print(){
	for(unsigned int i(0);i<path.size();i++){
		std::cout<<path[i];
		if(i<fname.size()){
			std::cout<<fname[i]<<ext[i];
		}
		std::cout<<std::endl;
	}
}
