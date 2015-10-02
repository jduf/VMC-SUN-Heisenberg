#include "Linux.hpp" 
#include <iostream>
#include <vector>

//int mkdir(const char *path, mode_t mode){
//    struct stat st;
//
//    if (stat(path, &st) != 0 && mkdir(path, mode) != 0 && errno != EEXIST){
		//return -1;
//    } else if(!S_ISDIR(st.st_mode)){
//        errno = ENOTDIR;
		//return -1;
//    }
	//return 0;
//}

//int mkpath(const char *path, mode_t mode){
//    int  status;
//    char *pp;
//    char *sp;
//    char *copypath = strdup(path);
//
//    status = 0;
//    pp = copypath;
//    while (status == 0 && (sp = strchr(pp, '/')) != 0){
//        if(sp != pp){
//            *sp = '\0';
//            status = do_mkdir(copypath, mode);
//            *sp = '/';
//        }
//        pp = sp + 1;
//    }
//    if (status == 0)
//        status = do_mkdir(path, mode);
//    free(copypath);
//    return (status);
//}
//
int *errno_p = __errno_location();

int main(){
	std::vector<double> a(1e9,1);
	Linux command;
	command.mkpath("aa/bb/cc/kk/ff/ee/");
	command("gnuplot -e 'plot sin(x)'",false);
	std::cout<<command.status()<<std::endl;
}
