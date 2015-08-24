#include "Directory.hpp"

/*{public method*/
void Directory::set(){
	filename_.clear();
	path_.clear();
	ext_.clear();
}

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
				path_.push_back(curr_dir);
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
						path_.push_back(curr_dir);
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
			std::string dir(entry->d_name);
			if (dir != "." && dir != ".."){
				path_.push_back(curr_dir);
				filename_.push_back(dir);
				ext_.push_back("/");
			} 
		}
		entry = readdir(dir_point);
	}
	closedir(dir_point);
}

void Directory::sort(){
	if(path_.size()>0){
		bool sort(true);
		while(sort){
			sort=false;
			for(unsigned int i(0); i<path_.size()-1;i++){
				if(path_[i]+filename_[i]+ext_[i] > path_[i+1]+filename_[i+1]+ext_[i+1]) { 
					sort= true;
					std::swap(path_[i],path_[i+1]);
					std::swap(filename_[i],filename_[i+1]);
					std::swap(ext_[i],ext_[i+1]);
				}
			}
		}
	} else {
		std::cerr<<__PRETTY_FUNCTION__<<" : the file list is empty"<<std::endl;
	}
}

void Directory::print(){
	for(unsigned int i(0);i<path_.size();i++){
		std::cout<<path_[i];
		if(i<filename_.size()){
			std::cout<<filename_[i]<<ext_[i];
		}
		std::cout<<std::endl;
	}
}
/*}*/

/*{private method*/
void Directory::split_ext(std::string f){
	if(f.find(".") != std::string::npos){
		size_t pos(f.find_last_of("."));
		filename_.push_back(f.substr(0,pos));
		ext_.push_back(f.substr(pos));
	} else {
		filename_.push_back(f);
		ext_.push_back("");
	}
}
/*}*/
