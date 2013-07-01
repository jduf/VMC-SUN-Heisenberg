#ifndef DEF_LINUX
#define DEF_LINUX

#include <cstdlib> 
#include <string>

class Linux {
	public:
		Linux(){};
		~Linux(){};

		void operator()(std::string c){ st=system(c.c_str()); }
		int status(){return st;}

	private:
		std::string command;
		int st;
};
#endif
