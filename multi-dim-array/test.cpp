/* @file test.cpp */

#include "Vector.hpp"
#include "Rand.hpp"
#include "omp.h"

unsigned int hypercube_element(int m, int n);
bool init(Vector<double>* A, unsigned int& i, unsigned int& j, unsigned int const& imax);
bool eval(Vector<double>* A, Vector<unsigned int>& idx);
bool eval_border(Vector<double>* A, Vector<unsigned int>& idx, unsigned int const& min0 = 0, unsigned int const& max0 = 0);
void normal(Vector<unsigned int> v);
void parallel(Vector<unsigned int> v);
unsigned int number_of_border_element(Vector<unsigned int> edges_length);

int main(){
	Rand<unsigned int> rnd_N(0,5);
	unsigned int N(rnd_N.get());
	Vector<unsigned int> v(N);
	Rand<unsigned int> rnd_l(0,500);
	for(unsigned int i(0);i<N;i++){ v(i) = rnd_l.get(); }
	v.set(4);
	v(0) = 80;
	v(1) = 10;
	v(2) = 10;
	v(3) = 10;
	std::cout<<v<<std::endl;
	std::cout<<number_of_border_element(v)<<std::endl;
	parallel(v);
	//normal(v);
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

bool eval_border(Vector<double>* A, Vector<unsigned int>& idx, unsigned int const& min0, unsigned int const& max0){
	for(unsigned int l(0);l<idx.size();l++){
		if(l<idx.size() && (idx(l) == 0 || idx(l) == A[l].size()-1) ) { 
#pragma omp critical
			std::cerr<<idx<<std::endl;
			l = idx.size();
		}
	}

	idx(0)++;
	if(min0==0 && max0==0){
		for(unsigned int i(0);i<idx.size();i++){
			if(idx(i) == A[i].size()){ 
				if(i+1 == idx.size()){ return false; }
				idx(i) = 0;
				idx(i+1)++;
			}
		}
	} else {
		if(idx(0) == max0){ 
			if(1 == idx.size()){ return false; }
			idx(0) = min0;
			idx(1)++;
		}
		for(unsigned int i(1);i<idx.size();i++){
			if(idx(i) == A[i].size()){ 
				if(i+1 == idx.size()){ return false; }
				idx(i) = 0;
				idx(i+1)++;
			}
		}
	}

	return eval_border(A,idx,min0,max0);
}

void normal(Vector<unsigned int> v){
	unsigned int N(v.size());
	Vector<double>* A(new Vector<double>[N]);
	for(unsigned int i(0);i<N;i++){ A[i].set(v(i),0); }

	unsigned int i(0);
	unsigned int j(0);
	while(init(A,i,j,N));

	Vector<unsigned int> idx;
	//idx.set(N,0);
	//while(eval(A,idx));
	idx.set(N,0);
	while(eval_border(A,idx));
}

void parallel(Vector<unsigned int> v){
	unsigned int N(v.size());
	Vector<double>* A(new Vector<double>[N]);
	for(unsigned int i(0);i<N;i++){ A[i].set(v(i),0); }

	unsigned int i(0);
	unsigned int j(0);
	while(init(A,i,j,N));

#pragma omp parallel
	{
		unsigned int min0(omp_get_thread_num()*A[0].size()/omp_get_num_threads());
		unsigned int max0((omp_get_thread_num()+1)*A[0].size()/omp_get_num_threads());
		Vector<unsigned int> idx;
		idx.set(N,0);
		idx(0) = min0;
		eval_border(A,idx,min0,max0);
	}
}

unsigned int hypercube_element(int m, int n){
	if(m==0 && n==0){ return 1; }
	if(m<0 || n<0 || n<m){ return 0; }
	return 2*hypercube_element(m,n-1)+hypercube_element(m-1,n-1);
}

unsigned int number_of_border_element(Vector<unsigned int> edges_length){
	unsigned int N(edges_length.size());
	bool hypercube(true);
	for(unsigned int i(1);i<N;i++){
		if(edges_length(0) != edges_length(i)){ 
			i = N;
			hypercube = false;
		}
	}
	unsigned int s(0);
	if(hypercube){
		for(unsigned int i(1);i<N+1;i++){
			s += hypercube_element(N-i,N)*(i%2?1:-1)*pow(10,N-i);
		}
	} else {
		for(unsigned int i(1);i<N+1;i++){
			Vector<unsigned int> coeff(my::comb(N,N-i,edges_length));
			if(i%2){ s += coeff.sum()*hypercube_element(N-i,N)/coeff.size(); }
			else{ s -= coeff.sum()*hypercube_element(N-i,N)/coeff.size(); }
		}
	}
	return s;
}
