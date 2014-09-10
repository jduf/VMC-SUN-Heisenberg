#include "Rand.hpp"
#include "Vector.hpp"
#include<omp.h>

void check_basic();
void check_openmp_mt();
void check_minimal_number_mt();

int main(){
	std::random_device rd;
	std::mt19937_64 mt(rd());

	unsigned int const n(16);
	unsigned int const m(2);
	unsigned int const N(4);
	unsigned int const M(n*m/N);
	Matrix<unsigned int> v(n,m);
	for(unsigned int i(0);i<N;i++){ 
		for(unsigned int j(0);j<M;j++){ 
			v.ptr()[i*M+j] = i;
		}
	}

	std::cout<<v.transpose()<<std::endl;
	std::cout<<std::endl;
	std::shuffle(v.ptr(),v.ptr()+n,mt);
	std::shuffle(v.ptr()+n,v.ptr()+v.size(),mt);
	std::cout<<v.transpose()<<std::endl;

	//check_basic();
	//check_openmp_mt();
	//check_minimal_number_mt();
}

void check_basic(){
	RandUnsignedInt rndui(0,10);
	for(unsigned int i(0);i<20;i++){
		std::cout<<rndui.get()<<std::endl;
	}
}

void check_openmp_mt(){
	std::cout<<"mt declared inside openmp"<<std::endl;
	Matrix<double> m(20,4);
#pragma omp parallel for num_threads(m.col())
	for(unsigned int i=0;i<m.col();i++){ 
		unsigned int thread(omp_get_thread_num());
		RandDouble rnd(0,1);
		for(unsigned int j(0);j<m.row();j++){ 
			m(j,thread)=rnd.get(); 
		}
	}
	std::cout<<m<<std::endl;

	std::cout<<"two different mt declared inside openmp"<<std::endl;
#pragma omp parallel for num_threads(m.col())
	for(unsigned int i=0;i<m.col();i++){ 
		unsigned int thread(omp_get_thread_num());
		RandDouble rnd0(0,1);
		RandDouble rnd1(0,1);
		for(unsigned int j(0);j<m.row();j++){ m(j,thread)=rnd1.get()-rnd0.get(); }
	}
	std::cout<<m<<std::endl;
}

void check_minimal_number_mt(){
	Time t;
	double min(10);
	double tmp;
	unsigned int i(0);
	std::random_device rd;
	std::mt19937_64 mt(rd());
	std::uniform_real_distribution<double> dist(0,1);

	while(!t.limit_reached(300)){ 
		tmp = dist(mt); 
		i++;
		if(tmp<min){ min = tmp; }
		if(t.progress(1)){ std::cout<<i<<" "<<min<<std::endl; }
	}
}
