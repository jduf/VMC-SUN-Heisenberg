#include "ExtractSystem.hpp"

ExtractSystem::ExtractSystem(std::string filename):
	file(filename),
	data(false)
{
	file.extract<Vector<unsigned int> >(ref);
}

void ExtractSystem::extract(){
	extract(data,data);
}

bool ExtractSystem::use_complex(){
	if(ref(1) == 1){ return false; }
	else { return true; }
}

void ExtractSystem::print(){
	std::cout<<"sts=" <<std::endl;
	std::cout<<data.get<Matrix<unsigned int> >("sts")<<std::endl;

	switch(ref(0)){
		case 2:
			{
				switch(ref(1)){
					case 1:
						{ 
							std::cout<<"EVec="<<std::endl;
							std::cout<<data.get<Matrix<double> >("EVec")<<std::endl;
							switch(ref(2)){
								case 0:{}break;
								case 1:{ std::cout<<"delta="<<data.get<double>("delta")<<std::endl; }break;
							}
						}break;
					default:{ std::cerr<<"ref = [02.] unknown"<<std::endl; }break;
				}
				std::cout<<"bc="<<data.get<double>("bc")<<std::endl;
			}break;
		case 3:
			{
				switch(ref(1)){
					case 0:
						{
							std::cout<<"nn="<<std::endl;
							std::cout<<data.get<Matrix<unsigned int> >("nn")<<std::endl;
							std::cout<<"cc="<<std::endl;
							std::cout<<data.get<Matrix<unsigned int> >("cc")<<std::endl;
							std::cout<<"sl="<<std::endl;
							std::cout<<data.get<Vector<unsigned int> >("sl")<<std::endl;
							std::cout<<"omega="<<std::endl;
							std::cout<<data.get<Matrix<std::complex<double> > >("omega")<<std::endl;
						}break;
					case 1:
						{ 
							std::cout<<"EVec="<<std::endl;
							std::cout<<data.get<Matrix<double> >("EVec")<<std::endl; 
							switch(ref(2)){ case 1:{ std::cout<<"mu="<<data.get<double>("mu")<<std::endl; }break;} //TriangleMu 
						}break;
					case 2:
						{ 
							std::cout<<"EVec="<<std::endl;
							std::cout<<data.get<Matrix<std::complex<double> > >("EVec")<<std::endl;
							switch(ref(2)) { case 4: { std::cout<<"phi="<<data.get<double>("phi")<<std::endl; }break;}//TrianglePhi 
						}break;
				}
			}break;
			std::cout<<"bc="<<data.get<double>("bc")<<std::endl;
			std::cout<<"Lx="<<data.get<unsigned int>("Lx")<<std::endl;
			std::cout<<"Ly="<<data.get<unsigned int>("Ly")<<std::endl;
		case 4:
			{
				switch(ref(1)){
					case 0:
						{
							std::cout<<"nn="<<std::endl;
							std::cout<<data.get<Matrix<unsigned int> >("nn")<<std::endl;
							std::cout<<"cc="<<std::endl;
							std::cout<<data.get<Matrix<unsigned int> >("cc")<<std::endl;
							std::cout<<"sl="<<std::endl;
							std::cout<<data.get<Vector<unsigned int> >("sl")<<std::endl;
							std::cout<<"omega="<<std::endl;
							std::cout<<data.get<Matrix<std::complex<double> > >("omega")<<std::endl;
						}break;
					case 1:
						{ 
							std::cout<<"EVec="<<std::endl;
							std::cout<<data.get<Matrix<double> >("EVec")<<std::endl; 
							switch(ref(2)){ case 1:{ std::cout<<"mu="<<data.get<double>("mu")<<std::endl; }break;} //SquareMu 
						}break;
					case 2:
						{ 
							std::cout<<"EVec="<<std::endl;
							std::cout<<data.get<Matrix<std::complex<double> > >("EVec")<<std::endl;
							switch(ref(2)){ case 4:{ std::cout<<"phi="<<data.get<double>("phi")<<std::endl; }break;}//SquarePhiFlux
						}break;
				}
				std::cout<<"bc="<<data.get<double>("bc")<<std::endl;
				std::cout<<"Lx="<<data.get<unsigned int>("Lx")<<std::endl;
				std::cout<<"Ly="<<data.get<unsigned int>("Ly")<<std::endl;
			}break;
		case 6:
			{
				switch(ref(1)){ case 1:{ std::cout<<data.get<Matrix<double> >("EVec")<<std::endl; }break; } //Honeycomb
				std::cout<<"bc="<<data.get<double>("bc")<<std::endl;
				std::cout<<"Lx="<<data.get<unsigned int>("Lx")<<std::endl;
				std::cout<<"Ly="<<data.get<unsigned int>("Ly")<<std::endl;
			}break;
	}

	std::cout<<"N ="<<data.get<unsigned int>("N")<<std::endl;
	std::cout<<"m ="<<data.get<unsigned int>("m")<<std::endl;
}

void ExtractSystem::extract(Container& input, Container& param){
	input.set("ref",ref);
	unsigned int m(0),N(0);
	file.extract<unsigned int>(N);
	file.extract<unsigned int>(m);
	input.set("N",N);
	input.set("m",m);
	input.set("n",N*m);
	param.set("N",N);
	param.set("m",m);
	param.set("n",N*m);
	file.extract<Matrix<unsigned int> >("sts",input);

	switch(ref(0)){
		case 2:
			{
				switch(ref(1)){
					case 1:
						{ 
							switch(ref(2)) {
								case 0:{ }break;
								case 1:{ file.extract<double>("delta",param); }break;
								default:{ std::cerr<<"ref = ["<<ref(0)<<ref(1)<<ref(2)<<"] unknown"<<std::endl; }break;
							}
							file.extract<Matrix<double> >("EVec",input);
						}break;
					default:{ std::cerr<<"ref = ["<<ref(0)<<ref(1)<<ref(2)<<"] unknown"<<std::endl; }break;
				}
				file.extract<double>("bc",param);
			}break;
		case 3:
			{
				switch(ref(1)){
					case 0:
						{//TriangleJastrow
							file.extract<Matrix<unsigned int> >("nn",input);
							file.extract<Matrix<unsigned int> >("cc",input);
							file.extract<Vector<unsigned int> >("sl",input);
							file.extract<Matrix<std::complex<double> > >("omega",input);
						}break;
					case 1:
						{ 
							switch(ref(2)){
								case 0:{ }break;//TriangleFermi
								case 1:{ file.extract<double>("mu",param); }break;//TriangleMu
								default:{ std::cerr<<"ref = ["<<ref(0)<<ref(1)<<ref(2)<<"] unknown"<<std::endl; }break;
							}
							file.extract<Matrix<double> >("EVec",input);
						}break;
					case 2:
						{ 
							switch(ref(2)){
								case 4:{ file.extract<double>("phi",param); }break;//TrianglePhi
								default:{ std::cerr<<"ref = ["<<ref(0)<<ref(1)<<ref(2)<<"] unknown"<<std::endl; }break;
							}
							file.extract<Matrix<std::complex<double> > >("EVec",input);
						}break;
					default:{ std::cerr<<"ref = ["<<ref(0)<<ref(1)<<ref(2)<<"] unknown"<<std::endl; }break;
				}
				file.extract<double>("bc",param);
				file.extract<unsigned int>("Lx",param);
				file.extract<unsigned int>("Ly",param);
			}break;
		case 4:
			{
				switch(ref(1)){
					case 0:
						{//SquareJastrow
							file.extract<Matrix<unsigned int> >("nn",input);
							file.extract<Matrix<unsigned int> >("cc",input);
							file.extract<Vector<unsigned int> >("sl",input);
							file.extract<Matrix<std::complex<double> > >("omega",input);
						}break;
					case 1:
						{ 
							switch(ref(2)){
								case 0:{ }break; //SquareFermi
								case 1:{ file.extract<double>("mu",param); }break; //SquareMu
								default:{ std::cerr<<"ref = ["<<ref(0)<<ref(1)<<ref(2)<<"] unknown"<<std::endl; }break;
							}
							file.extract<Matrix<double> >("EVec",input);
						}break;
					case 2:
						{ 
							switch(ref(2)){
								case 3:{ }break;//SquarePiFlux
								case 4:{ file.extract<double>("phi",param); }break;//SquarePhiFlux
								default:{ std::cerr<<"ref = ["<<ref(0)<<ref(1)<<ref(2)<<"] unknown"<<std::endl; }break;
							}
							file.extract<Matrix<std::complex<double> > >("EVec",input);
						}break;
					default:{ std::cerr<<"ref = ["<<ref(0)<<ref(1)<<ref(2)<<"] unknown"<<std::endl; }break;
				}
				file.extract<double>("bc",param);
				file.extract<unsigned int>("Lx",param);
				file.extract<unsigned int>("Ly",param);
			}break;
		case 6:
			{
				switch(ref(1)){
					case 1:
						{ 
							switch(ref(2)){
								case 0:{ }break; //HoneycombSU3
								case 1:{ }break; //HoneycombSU4
								default:{ std::cerr<<"ref = ["<<ref(0)<<ref(1)<<ref(2)<<"] unknown"<<std::endl; }break;
							}
							file.extract<Matrix<double> >("EVec",input);
						}break;
					default:{ std::cerr<<"ref = ["<<ref(0)<<ref(1)<<ref(2)<<"] unknown"<<std::endl; }break;
				}
				file.extract<double>("bc",param);
				file.extract<unsigned int>("Lx",param);
				file.extract<unsigned int>("Ly",param);
			}break;
		default:{ std::cerr<<"ref = ["<<ref(0)<<ref(1)<<ref(2)<<"] unknown"<<std::endl; }break;
	}
}
