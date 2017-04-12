/*!@file min.cpp */

#include "Directory.hpp"
#include "VMCPSO.hpp"
#include "VMCInterpolation.hpp"
#include "VMCSystematic.hpp"
#include "VMCExtract.hpp"
#include "VMCACiD.hpp"

void error();
//class Input {
	//public:
		//Input(Parseur& P, unsigned int i=0): 
			//which_obs(P.find("obs",i) ?(P.get_type(i)?P.get<std::vector<unsigned int> >(i):Vector<unsigned int>(1,P.get<unsigned int>(i))):0),
			//what(P.find("what",i,true)?(P.get_type(i)?P.get<std::vector<unsigned int> >(i):Vector<unsigned int>(1,P.get<unsigned int>(i))):0),
			//dEoE(P.check_get("dEoE",1e-5)),
			//tmax(P.check_get("tmax",10)),
			//maxiter(P.check_get("maxiter",10)),
			//save_in(P.check_get("save_in",std::string("./")))
			//{}
//
		//Vector<unsigned int> const which_obs;
		//Vector<unsigned int> const what;
		//double const dEoE;
		//unsigned int const tmax;
		//unsigned int const maxiter;
		//std::string const save_in;
//};

int main(int argc, char* argv[]){
	Time running_time;
	Parseur P(argc,argv);
	//Input in(P);

	unsigned int i;
	Vector<unsigned int> what(P.find("what",i,true)?(P.get_type(i)?P.get<std::vector<unsigned int> >(i):Vector<unsigned int>(1,P.get<unsigned int>(i))):0);
	Vector<unsigned int> which_obs(P.find("obs",i) ?(P.get_type(i)?P.get<std::vector<unsigned int> >(i):Vector<unsigned int>(1,P.get<unsigned int>(i))):0);

	if(what(0)<100){
		if(!P.locked()){
			VMCMinimization m(P);
			if(P.find("reset_PS")){ m.set_phase_space(P); }
			if(m.ready() && !P.locked()){
				for(unsigned int w(0);w<what.size();w++){
					switch(what(w)){
						case 0:
							{
								if(P.find("dEoE",i)){
									double dEoE(P.get<double>(i));
									if(P.find("E",i)){ m.refine(P.get<double>(i),dEoE,P.check_get("save_in",std::string("./"))); }
									if(P.find("nmin",i)){ m.refine(which_obs,dEoE,P.check_get("maxiter",10),P.check_get("tmax",10),P.get<unsigned int>(i),P.check_get("save_in",std::string("./"))); }
								} else { m.refine(P.check_get("save_in",std::string("./"))); }
							}break;
						case 1:
							{ m.find_and_run_minima(10,which_obs,P.check_get("dEoE",1e-5),P.check_get("save_in",std::string("./"))); }break;
						case 2:
							{ m.improve_bad_samples(P.check_get("dEoE",1e-5),P.check_get("save_in",std::string("./"))); }break;
						case 3:
							{ m.save_parameters(P.check_get("nbest",10),P.check_get("save_in",std::string("./"))); }break;
						case 4:
							{ m.run_parameters(P.get<std::string>("parameters_file"),which_obs,P.check_get("dEoE",1e-5),P.check_get("maxiter",10),P.check_get("save_in",std::string("./"))); }break;
						case 10:
							{
								VMCPSO pso(P,m);
								unsigned int loop(P.get<unsigned int>("loop"));
								if(!P.locked()){
									for(unsigned int l(0);l<loop;l++){
										pso.init(true);
										pso.run(P.check_get("dEoE",1e-5),P.check_get("maxiter",5),P.check_get("save_in",std::string("./")));
									}
								} else { std::cerr<<__PRETTY_FUNCTION__<<" : some argument are not correctly set"<<std::endl; }
							}break;
						case 20:
							{ VMCSystematic(m).run(P.check_get("dEoE",1e-5),P.check_get("maxiter",1),P.check_get("save_in",std::string("./"))); }break;
						case 21:
							{ VMCSystematic(m).plot(); }break;
						case 30:
							{
								std::cerr<<"WARNING : this method has never been really tested"<<std::endl;
								VMCInterpolation interp(m);
								unsigned int loop(P.get<unsigned int>("loop"));
								if(!P.locked()){
									for(unsigned int l(0);l<loop;l++){
										interp.init();
										interp.run(true,P.check_get("save_in",std::string("./")));
									}
								}
							}break;
						case 40:
							{
								Vector<unsigned int> d(P.get<std::vector<unsigned int> >("d"));
								Vector<double> param;
								unsigned int i;
								if(P.find("param",i)){ param = P.get<std::vector<double> >(i); }
								VMCACiD min(m,d,param);
								if(!P.find("norun")){
									if(!P.locked()){
										for(unsigned int j(0);j<10;j++){
											min.run(P.check_get("dEoE",1e-6),P.check_get("maxiter",10),P.check_get("tmax",10),P.check_get("maxstep",100),P.check_get("save_in",std::string("./")));
											m.save(P.check_get("save_in",std::string("./")));
										}
									}
								} else { min.display_param_and_xmean(param); }
							}break;
						default:
							{ error(); }
					}
				}
			}
		}
	} else {
		Directory dir;
		dir.search_files(P.get<std::string>("load"),P.get<std::string>("pathload"),false,true);
		if(dir.size() == 1){
			IOFiles in(dir[0],false,false);
			VMCExtract extract(in,P.get<unsigned int>("min_sort"), P.get<unsigned int>("max_sort"));

			if(!P.locked()){
				for(unsigned int w(0);w<what.size();w++){
					switch(what(w)){
						case 100:
							{ extract.save(P.check_get("save_in",std::string("./"))); }break;
						case 101:
							{ extract.refine(which_obs,P.check_get("dEoE",1e-5),P.check_get("maxiter",10),P.check_get("tmax",60),P.check_get("nmin",1000),P.check_get("save_in",std::string("./"))); }break;
						case 102:
							{
								std::string fname("VMCExtract-tmp");
								RSTFile rst("/tmp/",fname);

								List<MCSim> kept_samples;
								extract.analyse(P.check_get("save_in",std::string("./")),"test",kept_samples);

								kept_samples.set_target();
								while(kept_samples.target_next()){
									kept_samples.get().display_results("","","","","","/tmp/",&rst);
								}
								rst.save(false,true);
								Linux()(Linux::html_browser("/tmp/"+fname+".html"),true);
							}break;
						case 103:
							{
								Vector<unsigned int> d(P.get<std::vector<unsigned int> >("d"));
								Vector<double> param;
								unsigned int i;
								if(P.find("param",i)){ param = P.get<std::vector<double> >(i); }
								VMCACiD min(extract,d,param);
								if(!P.find("norun")){
									if(!P.locked()){
										for(unsigned int j(0);j<10;j++){
											min.run(P.check_get("dEoE",1e-6),P.check_get("maxiter",10),P.check_get("tmax",10),P.check_get("maxstep",100),P.check_get("save_in",std::string("./")));
											extract.save(P.check_get("save_in",std::string("./")));
										}
									}
								} else { min.display_param_and_xmean(param); }
							}break;
						case 104:
							{
								IOFiles in_ACiD(P.get<std::string>("ACiD"),false,false);
								VMCACiD min(extract,in_ACiD);
								if(!P.find("norun") && !P.locked()){
									for(unsigned int j(0);j<10;j++){
										min.run(P.check_get("dEoE",1e-6),P.check_get("maxiter",10),P.check_get("tmax",10),P.check_get("maxstep",100),P.check_get("save_in",std::string("./")));
										extract.save(P.check_get("save_in",std::string("./")));
									}
								}
							}break;
						default:
							{ error(); }
					}
				}
			}
		} else {
			if(dir.size()){
				std::cerr<<__PRETTY_FUNCTION__<<" : multiples files found :"<<std::endl;
				dir.print(std::cerr);
			} else {
				std::cerr<<__PRETTY_FUNCTION__<<" : no file found"<<std::endl;
			}
		}
	}
	std::cout<<"code has run for "<<my::convert_seconds(running_time.elapsed())<<" ("<<running_time.elapsed()<<"s)"<<std::endl;
}

void error(){
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown option 'what', options are :"<<std::endl;
	std::cerr<<"    - refine                      :  0"<<std::endl;
	std::cerr<<"    - find and run minima         :  1"<<std::endl;
	std::cerr<<"    - improve bad samples         :  2"<<std::endl;
	std::cerr<<"    - minimize with PSO           : 10"<<std::endl;
	std::cerr<<"    - measure all phase space     : 20"<<std::endl;
	std::cerr<<"    - plot all phase space        : 21"<<std::endl;
	std::cerr<<"    - d-dimensional interpolation : 30"<<std::endl;
	std::cerr<<"    - adaptive coordinate descent : 40"<<std::endl;
	std::cerr<<"    - save relevant samples       :100"<<std::endl;
	std::cerr<<"    - refine relevant samples     :101"<<std::endl;
	std::cerr<<"    - display relevant samples    :102"<<std::endl;
	std::cerr<<"    - perform a free minimization :103"<<std::endl;
}
