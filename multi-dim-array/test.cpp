/* @file test.cpp */

#include "Vector.hpp"
#include "omp.h"

unsigned int iter(0);
bool init(Vector<double>* A, unsigned int& i, unsigned int& j, unsigned int const& imax);
bool eval(Vector<double>* A, Vector<unsigned int>& idx);
bool eval_border(Vector<double>* A, Vector<unsigned int>& idx, bool min0=true, bool max0=true);
void normal();
void parallel();

int main(){
	//normal();
	parallel();
}

bool init(Vector<double>* A, unsigned int& i, unsigned int& j, unsigned int const& imax){
	if(j==A[i].size()){ 
		if(++i==imax){ return false; }
		j=0;
	} else { 
		A[i](j) = j*pow(10,i);
		j += 1;
	}
	return init(A,i,j,imax);
}

bool eval(Vector<double>* A, Vector<unsigned int>& idx){
	unsigned int tmp(0);
	for(unsigned int i(0);i<idx.size();i++){ tmp += A[i](idx(i)); }
	std::cout<<idx<<" ->"<<tmp<<std::endl;

	idx(0)++;
	for(unsigned int i(0);i<idx.size();i++){
		if(idx(i) == A[i].size()){ 
			if(i+1 == idx.size()){ return false; }
			idx(i) = 0;
			idx(i+1)++;
		}
	}
	return eval(A,idx);
}

bool eval_border(Vector<double>* A, Vector<unsigned int>& idx, bool min0, bool max0){
	if( (min0 && idx(0) == 0) || (max0 && idx(0) == A[0].size()-1) ) { 
		std::cout<<idx<<std::endl;
	} else {
		for(unsigned int l(1);l<idx.size();l++){
			if(l<idx.size() && (idx(l) == 0 || idx(l) == A[l].size()-1) ) { 
				std::cout<<idx<<std::endl;
				l = idx.size();
			}
		}
	}

	idx(0)++;
	for(unsigned int i(0);i<idx.size();i++){
		if(idx(i) == A[i].size()){ 
			if(i+1 == idx.size()){ return false; }
			idx(i) = 0;
			idx(i+1)++;
		}
	}
	return eval_border(A,idx);
}

void normal(){
	unsigned int N(3);
	Vector<double>* A(new Vector<double>[N]);
	A[0].set(10,0);
	A[1].set(10,0);
	A[2].set(10,0);
	//A[3].set(10,0);
	//for(unsigned int i(0);i<N;i++){ std::cout<<A[i].size()<<std::endl; }

	unsigned int i(0);
	unsigned int j(0);
	while(init(A,i,j,N));
	//for(unsigned int i(0);i<N;i++){ std::cout<<A[i]<<std::endl; }
	Vector<unsigned int> idx;
	//idx.set(N,0);
	//while(eval(A,idx));
	idx.set(N,0);
	while(eval_border(A,idx));
}

void parallel(){
	unsigned int N(3);
	Vector<double>* A(new Vector<double>[N]);
	A[0].set(80,0);
	A[1].set(10,0);
	A[2].set(10,0);
	//A[3].set(10,0);
	//for(unsigned int i(0);i<N;i++){ std::cout<<A[i].size()<<std::endl; }

	unsigned int i(0);
	unsigned int j(0);
	while(init(A,i,j,N));
	//for(unsigned int i(0);i<N;i++){ std::cout<<A[i]<<std::endl; }
	//idx.set(N,0);
	//while(eval(A,idx));

	Vector<double> backup_A0(A[0]);
#pragma omp parallel
	{
		Vector<unsigned int> idx;
		idx.set(N,0);
		Vector<double>* B(new Vector<double>[N]);
		B[0].set(A[0].size()/8);
		for(unsigned int i(0);i<B[0].size();i++){
			B[0](i) = A[0](i+omp_get_thread_num()*B[0].size());
		}
		for(unsigned int i(1);i<N;i++){ B[i] = A[i]; }
#pragma omp critical
		eval_border(B,idx,(omp_get_thread_num()==0?true:false),(omp_get_thread_num()==7?true:false));
		delete[] B;
	}
}
