#ifndef DEF_LINUX
#define DEF_LINUX

#include <cstdlib> 
#include <string>

class Linux {
	public:
		Linux(){};
		~Linux(){};

		void operator()(std::string c){ st=system(c.c_str()); }
		std::string pwd(){ return std::string(get_current_dir_name()) + '/'; }
		int status(){return st;}

	private:
		std::string command;
		int st;
};
#endif
