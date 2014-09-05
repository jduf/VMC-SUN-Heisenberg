#ifndef DEF_LINUX
#define DEF_LINUX

#include <cstdlib> 
#include <string>
#include <unistd.h>

class Linux {
	public:
		/*!Constructor*/
		Linux(){}
		/*!Destructor*/
		~Linux(){}

		/*!Execute a UNIX command and get its exit value*/
		void operator()(std::string c){ ev_=system(c.c_str()); }
		/*!Returns exit value of the last command*/
		int status(){return ev_;}
		/*!Returns a string containing the current path*/
		std::string pwd(){ return std::string(get_current_dir_name()) + '/'; }

	private:
		/*!Forbids copy*/
		Linux(Linux const& l);
		/*!Forbids assignment*/
		Linux& operator=(Linux const& l);

		int ev_;//!< exit value of the last UNIX command
};
#endif
