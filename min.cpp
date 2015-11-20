/*! @file min.cpp */

#include "VMCPSO.hpp"
#include "VMCInterpolation.hpp"
#include "VMCSystematic.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	
	unsigned int i;
	VMCMinimization m(P);
	if(m.ready()){
		switch(P.find("what",i,true)?P.get<unsigned int>(i):666){
			case 0:
				{
					VMCPSO m1(P,m,P.get<bool>("symmetry"));
					VMCInterpolation m2(m);
					unsigned int loop(P.get<unsigned int>("loop"));
					if(!P.locked()){
						for(unsigned int i(0);i<loop;i++){
							m1.init(true, 3*i<loop?false:true);
							m1.run();
							m1.save();

							m2.init();
							m2.run(true);
							m2.save();
						}

						m.get_header().title("Minimization",'>');
						m.complete_analysis(1e-5);
						m.save();

						m.refine();
						m.save();

						m.set_tmax(20);
						m.find_and_run_minima(10,-1,1e-4);
						m.save();
					} else { std::cerr<<__PRETTY_FUNCTION__<<" : some argument are not corretly set"<<std::endl; }
				}break;
			case 1:
				{
					unsigned int i(0);
					unsigned int j(0);
					if(P.find("E",i,false),P.find("dE",j,false)){ m.refine(P.get<double>(i),P.get<double>(j)); }
					else{ m.refine(); }
					m.save();
				}break;
			case 2:
				{
					m.set_tmax(30);
					m.find_and_run_minima(10,-1,1e-4);
					m.save();
				}break;
			case 3:
				{
					m.set_phase_space(P);
					m.save();
				}break;
			case 4:
				{
					IOFiles out("test.jdbin",true);
					m.find_save_and_plot_minima(10,out);
					VMCInterpolation m2(m);
					m2.plot();
				}break;
			case 5:
				{
					Vector<double> param(11,1.0);
					Matrix<int> sym(9,3);
					Vector<double> J(P.get<std::vector<double> >("Jp"));
					int A(J(0)>J(1)?0:-1); //link between sites 0-1 (rung) 
					int B(J(0)>J(1)?-1:0); //link between sites 0-2 (leg)

					sym(0,0) = 1;
					sym(0,1) = B;
					sym(0,2) = -1;

					sym(1,0) = 2;
					sym(1,1) = B;
					sym(1,2) = 1;

					sym(2,0) = 3;
					sym(2,1) = 0;
					sym(2,2) = 1;

					sym(3,0) = 4;
					sym(3,1) = B;
					sym(3,2) = -1;

					sym(4,0) = 5;
					sym(4,1) = B;
					sym(4,2) = 1;

					sym(5,0) = 6;
					sym(5,1) = 0;
					sym(5,2) = 1;

					sym(6,0) = 7;
					sym(6,1) = B;
					sym(6,2) = -1;

					sym(7,0) = 9;
					sym(7,1) = A;
					sym(7,2) = 1;

					sym(8,0) = 10;
					sym(8,1) = 8;
					sym(8,2) = -1;

					VMCSystematic m3(m,param,sym,0,8);
					m.set_tmax(30);
					m3.run(0,1e-6,20);
					m3.save();
					m3.plot();
					//m3.test();
				}break;
			case 6:
				{
					Vector<double> param(11,1.0);
					Vector<double> J(P.get<std::vector<double> >("Jp"));
					int A(J(0)>J(1)?0:-1); //link between sites 0-1 (rung) 
					int B(J(0)>J(1)?-1:0); //link between sites 0-2 (leg)

					Matrix<int> sym(9,3);
					switch(P.get<unsigned int>("s")){
						case 0:
							{
								/*facing tetramerization*/
								sym(0,0) = 1;
								sym(0,1) = B;
								sym(0,2) = 1;

								sym(1,0) = 2;
								sym(1,1) = B;
								sym(1,2) = 1;

								sym(2,0) = 3;
								sym(2,1) = 0;
								sym(2,2) = 1;

								sym(3,0) = 4;
								sym(3,1) = B;
								sym(3,2) = 1;

								sym(4,0) = 5;
								sym(4,1) = B;
								sym(4,2) = 1;

								sym(5,0) = 6;
								sym(5,1) = 0;
								sym(5,2) = 1;

								sym(6,0) = 7;
								sym(6,1) = B;
								sym(6,2) = 1;

								sym(7,0) = 9;
								sym(7,1) = A;
								sym(7,2) = 1;

								sym(8,0) = 10;
								sym(8,1) = 8;
								sym(8,2) = 1;
							}break;
						case 1:
							{
								/*facing tetramerization with pi*/
								sym(0,0) = 1;
								sym(0,1) = B;
								sym(0,2) = -1;

								sym(1,0) = 2;
								sym(1,1) = B;
								sym(1,2) = 1;

								sym(2,0) = 3;
								sym(2,1) = 0;
								sym(2,2) = 1;

								sym(3,0) = 4;
								sym(3,1) = B;
								sym(3,2) = -1;

								sym(4,0) = 5;
								sym(4,1) = B;
								sym(4,2) = 1;

								sym(5,0) = 6;
								sym(5,1) = 0;
								sym(5,2) = 1;

								sym(6,0) = 7;
								sym(6,1) = B;
								sym(6,2) = -1;

								sym(7,0) = 9;
								sym(7,1) = A;
								sym(7,2) = 1;

								sym(8,0) = 10;
								sym(8,1) = 8;
								sym(8,2) = -1;
							}break;
						case 2:
							{
								/*shifted by 2 tetramerization*/
								tmp(0,0) = 1;
								tmp(0,1) = B;
								tmp(0,2) = 1;

								tmp(1,0) = 2;
								tmp(1,1) = B;
								tmp(1,2) = 1;

								tmp(2,0) = 3;
								tmp(2,1) = A;
								tmp(2,2) = 1;

								tmp(3,0) = 4;
								tmp(3,1) = 8;
								tmp(3,2) = 1;

								tmp(4,0) = 5;
								tmp(4,1) = B;
								tmp(4,2) = 1;

								tmp(5,0) = 6;
								tmp(5,1) = A;
								tmp(5,2) = 1;

								tmp(6,0) = 7;
								tmp(6,1) = B;
								tmp(6,2) = 1;

								tmp(7,0) = 9;
								tmp(7,1) = A;
								tmp(7,2) = 1;

								tmp(8,0) = 10;
								tmp(8,1) = B;
								tmp(8,2) = 1;
							}break;
						case 3:
							{
								/*shifted by 2 tetramerization with pi*/
								tmp(0,0) = 1;
								tmp(0,1) = B;
								tmp(0,2) = -1;

								tmp(1,0) = 2;
								tmp(1,1) = B;
								tmp(1,2) = 1;

								tmp(2,0) = 3;
								tmp(2,1) = A;
								tmp(2,2) = 1;

								tmp(3,0) = 4;
								tmp(3,1) = 8;
								tmp(3,2) = -1;

								tmp(4,0) = 5;
								tmp(4,1) = B;
								tmp(4,2) = 1;

								tmp(5,0) = 6;
								tmp(5,1) = A;
								tmp(5,2) = 1;

								tmp(6,0) = 7;
								tmp(6,1) = B;
								tmp(6,2) = -1;

								tmp(7,0) = 9;
								tmp(7,1) = A;
								tmp(7,2) = 1;

								tmp(8,0) = 10;
								tmp(8,1) = B;
								tmp(8,2) = -1;
							}break;
					}
					VMCSystematic m3(m,param,sym,0,8);
					m.set_tmax(150);
					m3.rerun(10,0,1e-6,20);
					m3.save();
					m3.plot();
					//m3.test();
				}break;
			default:
				{
					std::cerr<<__PRETTY_FUNCTION__<<" : unknown option 'what', options are :"<<std::endl;
					std::cerr<<"    - complete run                 : 0"<<std::endl;
					std::cerr<<"    - refine + save                : 1"<<std::endl;
					std::cerr<<"    - find and run minima + save   : 2"<<std::endl;
					std::cerr<<"    - redefine phase space + save  : 3"<<std::endl;
					std::cerr<<"    - plot                         : 4"<<std::endl;
					std::cerr<<"    - systematic run + save + plot : 5"<<std::endl;
				}
		}
	}
}
