/*!@file min.cpp */

#include "VMCPSO.hpp"
#include "VMCInterpolation.hpp"
#include "VMCSystematic.hpp"
#include "VMCExtract.hpp"

void error();

int main(int argc, char* argv[]){
	Parseur P(argc,argv);

	unsigned int i;
	P.find("what",i,true);
	Vector<unsigned int> what(P.get_type(i)?P.get<std::vector<unsigned int> >(i):Vector<unsigned int>(1,P.get<unsigned int>(i)));
	double dEoE(P.find("dEoE",i,false)?P.get<double>(i):0.01);
	unsigned int maxiter(P.find("maxiter",i,false)?P.get<unsigned int>(i):1);
	Vector<unsigned int> obs(P.find("obs",i,false)?(P.get_type(i)?P.get<std::vector<unsigned int> >(i):Vector<unsigned int>(1,P.get<unsigned int>(i))):0);

	if(what(0)<100){
		if(!P.locked()){
			VMCMinimization m(P);
			if(P.find("reset_PS",i,false)){ m.set_phase_space(P); }
			if(m.ready() && !P.locked()){
				for(unsigned int w(0);w<what.size();w++){
					switch(what(w)){
						case 0:
							{
								if(P.find("dEoE",i,false)){
									if(P.find("E",i,false)){ m.refine(P.get<double>(i),dEoE); }
									if(P.find("nmin",i,false)){ m.refine(P.get<unsigned int>(i),obs,dEoE,maxiter); }
								} else { m.refine(); }
							}break;
						case 1:
							{ m.find_and_run_minima(10,obs,1e-5); }break;
						case 2:
							{ m.improve_bad_samples(dEoE); }break;
						case 3:
							{ m.save_parameters(P); }break;
						case 4:
							{ m.run_parameters(P); }break;
						case 10:
							{
								VMCPSO pso(P,m);
								unsigned int loop(P.get<unsigned int>("loop"));
								if(!P.locked()){
									for(unsigned int l(0);l<loop;l++){
										pso.init(true);
										pso.run();
										pso.complete_analysis(1e-5);
										pso.refine(30,false,1e-5,5);
										pso.save();
									}
								} else { std::cerr<<__PRETTY_FUNCTION__<<" : some argument are not correctly set"<<std::endl; }
							}break;
						case 20:
							{
								VMCSystematic sys(m);
								sys.run(obs,dEoE,maxiter);
								sys.save();
							}break;
						case 21:
							{ VMCSystematic(m).plot(); }break;
						case 30:
							{
								std::cerr<<"WARNING : this method has never really optimized anything"<<std::endl;
								VMCInterpolation interp(m);
								unsigned int loop(P.get<unsigned int>("loop"));
								if(!P.locked()){
									for(unsigned int l(0);l<loop;l++){
										interp.init();
										interp.run(true);
										interp.save();
									}
								}
							}break;
						case 40:{ /*kim*/ }break;
						default:
							{ error(); }
					}
				}
			}
		}
	} else {
		IOFiles in(P.get<std::string>("load"),false,false);
		unsigned int ttotal(P.get<unsigned int>("ttotal"));
		std::string dirname(P.get<std::string>("dirname"));

		VMCExtract extract(in,P.get<unsigned int>("min_sort"),P.get<unsigned int>("max_sort"));
		if(!P.locked()){
			for(unsigned int w(0);w<what.size();w++){
				switch(what(w)){
					case 100:
						{ extract.save(dirname); }break;
					case 101:
						{
							extract.refine(obs,dEoE,ttotal);
							extract.save(dirname);
						}break;
					case 102:
						{
							std::string fname("VMCExtract-tmp");
							RSTFile rst("/tmp/",fname);

							List<MCSim> kept_samples;
							extract.analyse(dirname,"test",kept_samples);

							kept_samples.set_target();
							while(kept_samples.target_next()){
								kept_samples.get().display_results("","","","","","/tmp/",&rst);
							}
							rst.save(false,true);
							Linux()(Linux::html_browser("/tmp/"+fname+".html"),true);
						}break;
					default:
						{ error(); }break;
				}
			}
		}
	}
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
	std::cerr<<"    - save relevant samples       :100"<<std::endl;
	std::cerr<<"    - refine relevant samples     :101"<<std::endl;
	std::cerr<<"    - display relevant samples    :102"<<std::endl;
}
