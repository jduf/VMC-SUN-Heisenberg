#include "Rand.hpp"
#include "Matrix.hpp"
#include<omp.h>

void check_basic();
void check_time();
void check_openmp();

int main(){
	check_basic();
	//check_time();
	//check_openmp();
}

void check_basic(){
	Rand rdm(1e4);

	std::cout<<rdm.get()<<std::endl;
	std::cout<<rdm.get(10)<<std::endl;
	std::cout<<rdm.get(2)<<std::endl;

	Rand rnd(10,1802,9373);
	for(unsigned int i(0);i<20000;i++){
		rnd.get();
	}
	std::cout<<"if the code works properly, the two following lines should be identical (cf Rand.cpp)"<<std::endl;
	for(unsigned int i(0);i<6;i++){
		std::cout<<int(rnd.get()*4096*4096)<<" ";
	}
	std::cout<<std::endl<<6533892<<" "<<14220222<<" "<<7275067 <<" "<<6172232<<" "<<8354498 <<" "<<10633180<<std::endl;
}

void check_time(){
	Time t;
	Rand rnd1(10);
	for(unsigned int i(0);i<1e9;i++){ rnd1.get(); }
	std::cout<<t.elapsed()<<std::endl;

	t.set();
	Rand rnd2(1e6);
	for(unsigned int i(0);i<1e9;i++){ rnd2.get(); }
	std::cout<<t.elapsed()<<std::endl;
}

void check_openmp(){
	std::cout<<"Rand declared inside openmp"<<std::endl;
	Matrix<double> m(20,4);
#pragma omp parallel for num_threads(m.col())
	for(unsigned int i=0;i<m.col();i++){ 
		unsigned int thread(omp_get_thread_num());
		Rand rnd0(1e4,thread);
		for(unsigned int j(0);j<m.row();j++){ 
			m(j,thread)=rnd0.get(); 
		}
	}
	std::cout<<m<<std::endl;

	std::cout<<"two Rand's declared inside openmp with the same/different seed"<<std::endl;
	std::cout<<"(is problematic if same seed)"<<std::endl;
#pragma omp parallel for num_threads(m.col())
	for(unsigned int i=0;i<m.col();i++){ 
		unsigned int thread(omp_get_thread_num());
		Rand rnd1(1e4,thread);
		Rand rnd2(1e4,thread);
		//Rand rnd2(1e4,omp_get_thread_num()+rnd1.get(10)+1);
		for(unsigned int j(0);j<m.row();j++){ m(j,thread)=rnd1.get()-rnd2.get(); }
	}
	std::cout<<m<<std::endl;

	std::cout<<"Rand declared outside openmp (different threads try to access"<<std::endl;
	std::cout<<"the same rnd is problematic because the generation of random numbers"<<std::endl;
	std::cout<<"takes time, more visible if rnd3(>1e7))"<<std::endl;
	Rand rnd3(1e7);
#pragma omp parallel for num_threads(m.col())
	for(unsigned int i=0;i<m.col();i++){ 
		unsigned int thread(omp_get_thread_num());
		for(unsigned int j(0);j<m.row();j++){ m(j,thread)=rnd3.get(); }
	}

	std::cout<<m<<std::endl;

}
