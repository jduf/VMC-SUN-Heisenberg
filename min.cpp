/*!@file min.cpp */

#include "VMCPSO.hpp"
#include "VMCInterpolation.hpp"
#include "VMCSystematic.hpp"
#include "VMCExtract.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);

	unsigned int i;
	P.find("what",i,true);
	Vector<unsigned int> what(P.get_type(i)?P.get<std::vector<unsigned int> >(i):Vector<unsigned int>(1,P.get<unsigned int>(i)));
	if(!P.locked()){
		double dEoE(P.find("dEoE",i,false)?P.get<double>(i):0.01);
		unsigned int maxiter(P.find("maxiter",i,false)?P.get<unsigned int>(i):1);
		Vector<unsigned int> obs(P.find("obs",i,false)?(P.get_type(i)?P.get<std::vector<unsigned int> >(i):Vector<unsigned int>(1,P.get<unsigned int>(i))):0);
		if(what(0)<20){
			VMCMinimization m(P);
			if(P.find("reset_PS",i,false)){ m.set_phase_space(P); }
			if(m.ready() && !P.locked()){
				for(unsigned int w(0);w<what.size();w++){
					switch(what(w)){
						case 0:
							{
								VMCPSO m1(P,m);
								VMCInterpolation m2(m);
								unsigned int loop(P.get<unsigned int>("loop"));
								if(!P.locked()){
									for(unsigned int l(0);l<loop;l++){
										m1.init(true);
										m1.run();
										m1.refine(30,false,1e-5,5);
										m1.save();

										//m.refine(10,0,1e-5,5);
										//m.explore_around_minima(10,0,1e-5,0.05);
										//m.complete_analysis(1e-5);
										//m.save();
										//m2.init();
										//m2.run(true);
										//m2.save();
										//
									}

									m.complete_analysis(1e-5);
									m.refine();
								} else { std::cerr<<__PRETTY_FUNCTION__<<" : some argument are not corretly set"<<std::endl; }
							}break;
						case 1:
							{
								m.complete_analysis(1e-5);
								unsigned int j(0);
								if(P.find("dEoE",j,false)){
									if(P.find("E",i,false)){ m.refine(P.get<double>(i),P.get<double>(j)); }
									if(P.find("nmin",i,false)){ m.refine(P.get<unsigned int>(i),obs,P.get<double>(j),maxiter); }
								} else { m.refine(); }
							}break;
						case 2:
							{
								m.find_and_run_minima(10,obs,1e-5);
								m.save();
							}break;
						case 3:
							{
								m.improve_bad_samples(dEoE);
								m.save();
								m.refine();
							}break;
						case 5:
							{
								VMCSystematic m3(m);
								m3.run(obs,dEoE,maxiter);
								m3.save();
								m3.plot();
							}break;
						case 6:
							{
								m.save_parameters(P.get<unsigned int>("nparam"));
							}break;
						case 7:
							{
								m.run_parameters(P);
								m.save();
							}break;
						case 9:
							{
								VMCSystematic m3(m);
								m3.plot();
							}break;
						case 10:
							{
								IOFiles out("out.dat",true,false);
								m.find_save_and_plot_minima(10,out,"./","bla");
							}break;
						case 11:
							{ m.check(1000); }break;
						default:
							{
								std::cerr<<__PRETTY_FUNCTION__<<" : unknown option 'what', options are :"<<std::endl;
								std::cerr<<"    - complete run (PSO)           : 0"<<std::endl;
								std::cerr<<"    - refine                       : 1"<<std::endl;
								std::cerr<<"    - find and run minima + save   : 2"<<std::endl;
								std::cerr<<"    - redefine phase space + save  : 3"<<std::endl;
								std::cerr<<"    - PSO run with symmetries      : 4"<<std::endl;
								std::cerr<<"    - systematic run + save + plot : 5"<<std::endl;
								std::cerr<<"    - save parameters              : 6"<<std::endl;
								std::cerr<<"    - run parameters               : 7"<<std::endl;
							}
					}
				}
			}
		} else {
			if(P.find("min_sort",i,true) && P.find("max_sort",i,true) && P.find("load",i,true) && P.find("dirname",i,true)){
				std::string dirname(P.get<std::string>("dirname"));
				IOFiles in(P.get<std::string>("load"),false,false);
				VMCExtract m4(in,P.get<unsigned int>("min_sort"),P.get<unsigned int>("max_sort"));
				for(unsigned int w(0);w<what.size();w++){
					switch(what(w)){
						case 20:
							{ m4.save(dirname); }break;
						case 21:
							{
								m4.refine(obs,dEoE,P.get<unsigned int>("ttotal"));
								for(unsigned int i(0);i<obs.size();i++){
									if(obs(i)==1){ dirname += "_bond-energy"; }
									if(obs(i)==2){ dirname += "_long-range-correlation"; }
									if(obs(i)==3){ dirname += "_color-occupation"; }
								}
								m4.save(dirname);
							}break;
						case 22:
							{
								List<MCSim> kept_samples;
								m4.analyse("./","test-m4",kept_samples);
								kept_samples.set_target();

								std::string fname("bla");
								RSTFile rst("/tmp/",fname);
								while(kept_samples.target_next()){
									kept_samples.get().display_results("","","","","","/tmp/",&rst);
								}
								//rst.text(iof.get_header());
								rst.save(false,true);
								Linux()(Linux::html_browser("/tmp/"+fname+".html"),true);
							}break;
					}
				}
			} else { std::cerr<<__PRETTY_FUNCTION__<<" : missing arguments"<<std::endl; }
		}
	}
}
