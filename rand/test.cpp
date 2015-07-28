#include "Rand.hpp"
#include "Vector.hpp"
#include <chrono>
#include <omp.h>

void check_basic();
void check_openmp_mt();
void check_openmp_mt_time();

void check_minimal_number_mt();
void check_shuffle();

void check_array();

class TMP{
	public:
		TMP():rnd_(0.0,1.0){}
		double get() { return rnd_.get(); }

	private:
		Rand<double> rnd_;
};

int main(){
	//check_basic();
	//check_openmp_mt();
	check_openmp_mt_time();
	
	//check_minimal_number_mt();
	//check_shuffle();
	
	//check_array();
}

void check_basic(){
	Rand<unsigned int> rndui(0,10);
	for(unsigned int i(0);i<20;i++){
		std::cout<<rndui.get()<<std::endl;
	}
	std::cout<<std::endl;
}

void check_openmp_mt(){
	std::cout<<"mt declared inside openmp"<<std::endl;
	Matrix<double> m(20,omp_get_max_threads());
#pragma omp parallel for num_threads(m.col())
	for(unsigned int i=0;i<m.col();i++){ 
		unsigned int thread(omp_get_thread_num());
		Rand<double> rnd(0,1);
		for(unsigned int j(0);j<m.row();j++){ 
			m(j,thread)=rnd.get(); 
		}
	}
	std::cout<<m<<std::endl<<std::endl;

	std::cout<<omp_get_max_threads()<<" different mt declared inside openmp"<<std::endl;
#pragma omp parallel for num_threads(m.col())
	for(unsigned int i=0;i<m.col();i++){ 
		unsigned int thread(omp_get_thread_num());
		Rand<double> rnd0(0,1);
		Rand<double> rnd1(0,1);
		for(unsigned int j(0);j<m.row();j++){ m(j,thread)=rnd1.get()-rnd0.get(); }
	}
	std::cout<<m<<std::endl<<std::endl;

	std::cout<<"one mt declared outside openmp (note the same value if #pragma critical is removed) : takes too long"<<std::endl;
	Rand<double> rnd_out(0,1);
#pragma omp parallel for num_threads(m.col())
	for(unsigned int i=0;i<m.col();i++){ 
		unsigned int thread(omp_get_thread_num());
#pragma omp critical
		{
			for(unsigned int j(0);j<m.row();j++){ m(j,thread)=rnd_out.get(); }
		}
	}
	std::cout<<m<<std::endl<<std::endl;

	std::cout<<"manual array of mt declared outside ompen, one mt for each thread"<<std::endl;
	unsigned int N(omp_get_max_threads());
	Rand<double>** rnd_array(new Rand<double>*[N]);
	for(unsigned int i(0);i<N;i++){
		rnd_array[i] = new Rand<double>(i,i+1);
	}
#pragma omp parallel for num_threads(m.col())
	for(unsigned int i=0;i<m.col();i++){ 
		unsigned int thread(omp_get_thread_num());
		for(unsigned int j(0);j<m.row();j++){ m(j,thread)=rnd_array[thread]->get(); }
	}
	std::cout<<m<<std::endl<<std::endl;
	for(unsigned int i(0);i<N;i++){ delete rnd_array[i]; }
	delete[] rnd_array;
	rnd_array = NULL;

	std::cout<<"auto array of mt declared outside ompen, one mt for each thread"<<std::endl;
	RandArray<double> rnd_array_class(N);
	for(unsigned int i(0);i<N;i++){
		rnd_array_class.set(i,i,i+1);
	}
#pragma omp parallel for num_threads(m.col())
	for(unsigned int i=0;i<m.col();i++){ 
		unsigned int thread(omp_get_thread_num());
		for(unsigned int j(0);j<m.row();j++){ m(j,thread)=rnd_array_class.get(thread); }
	}
	std::cout<<m<<std::endl;
}

void check_openmp_mt_time(){
	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;
	unsigned int duration;

	unsigned iter(1e7);
	unsigned int N(omp_get_max_threads());
	std::cout<<"will run with "<<N<<" threads"<<std::endl;

	RandArray<double> rnd_array_class(N);
	for(unsigned int i(0);i<N;i++){ rnd_array_class.set(i,i,i+1); }
	Rand<double>** rnd_array(new Rand<double>*[N]);
	for(unsigned int i(0);i<N;i++){ rnd_array[i] = new Rand<double>(i,i+1); }

	unsigned int d1(0);
	unsigned int d2(0);
	for(unsigned int k(0);k<80;k++){
		t1 = std::chrono::high_resolution_clock::now();
#pragma omp parallel for
		for(unsigned int i=0;i<N;i++){ 
			for(unsigned int j(0);j<iter;j++){ rnd_array[i]->get(); }
		}
		t2 = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
		d1+= duration;
		std::cout<<"Rand**    : "<<duration<<std::endl;

		t1 = std::chrono::high_resolution_clock::now();
#pragma omp parallel for 
		for(unsigned int i=0;i<N;i++){ 
			for(unsigned int j(0);j<iter;j++){ rnd_array_class.get(i); }
		}
		t2 = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
		d2+= duration;
		std::cout<<"RandArray : "<<duration<<std::endl;
	}
	std::cout<<"total time "<<d1<<" "<<d2<<std::endl;

	for(unsigned int i(0);i<N;i++){ delete rnd_array[i]; }
	delete[] rnd_array;
	rnd_array = NULL;
}

void check_shuffle(){
	std::random_device rd;
	std::mt19937_64 mt(rd());

	unsigned int const n(16);
	unsigned int const m(2);
	Matrix<unsigned int> v(n,m);
	for(unsigned int i(0);i<n;i++){ 
		for(unsigned int j(0);j<m;j++){ 
			v(i,j) = j;
		}
	}

	std::cout<<v.transpose()<<std::endl;
	std::cout<<std::endl;
	std::shuffle(v.ptr(),v.ptr()+n*m,mt);
	std::cout<<v.transpose()<<std::endl;
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
	}
}

void check_array(){
	unsigned int N(10);
	TMP* tmp(new TMP[N]);
	TMP* test(tmp);
	for(unsigned int i(0);i<N;i++){
		std::cout<<tmp[i].get()<<" "<<test[i].get()<<std::endl;
	}
	delete[] tmp;
}
