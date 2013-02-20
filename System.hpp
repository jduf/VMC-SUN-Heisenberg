#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "Vecteur.hpp"
#include "Matrice.hpp"
#include "Lapack.hpp"
#include "Array2D.hpp"

#include <cstdlib>
#include <cmath>

template<typename T>
class System{
	public:
		System(unsigned int N_spin, unsigned int N_m, Array2D<unsigned int> const& sts, Matrice<T> const& H, Matrice<T> const& EVec);
		~System();

		unsigned int const N_spin, N_m, N_site;
		Array2D<unsigned int> const sts;
		
		void update();
		T ratio(bool update=false);
		void swap();
		void swap(unsigned int a, unsigned int b);

		unsigned int const& operator[](unsigned int const& i) const { return wis[i];};

		void print();

	private:
		System();
		System& operator=(System const & S);
		Matrice<T> *A;
		Matrice<T> *Ainv;
		Matrice<T> tmp_mat;
		Matrice<T> const H;
		unsigned int *s;
		unsigned int *wis;
		unsigned int cc[2]; // column changed
		unsigned int mc[2]; // matrix changed
		T w[2];

		void init_state(Matrice<T> const& EVec);
};

/*Constructors and destructor*/
/*{*/
template<typename T>
System<T>::System(unsigned int N_spin, unsigned int N_m, Array2D<unsigned int> const& sts, Matrice<T> const& H, Matrice<T> const& EVec):
	N_spin(N_spin),
	N_m(N_m),
	N_site(N_m*N_spin),
	sts(sts),
	A(new Matrice<T>[N_spin]),
	Ainv(new Matrice<T>[N_spin]),
	tmp_mat(N_m,0.0),
	H(H),
	s(new unsigned int[N_spin*N_m]),
	wis(new unsigned int[N_spin*N_m])
{
	srand(time(NULL)^(getpid()<<16));
	init_state(EVec);
	for(unsigned int i(0);i<2;i++){
		w[i] = 0;
		cc[i] = 0;
		mc[i] = 0;
	}
}

template<typename T>
System<T>::~System(){
	delete[] s;
	delete[] wis;
	delete[] A;
	delete[] Ainv;
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename T> 
T System<T>::ratio(bool update){
	if(update){
		if(mc[0] == mc[1]){
			return -1;
		} else { 
			w[0] = 0.0;
			w[1] = 0.0;
			for(unsigned int i(0);i<N_m;i++){
				w[0] += Ainv[mc[0]](cc[0],i)*A[mc[1]](i,cc[1]);
				w[1] += Ainv[mc[1]](cc[1],i)*A[mc[0]](i,cc[0]);
			}
			//std::cout<<alpha<<" "<<beta<<" "<<H(alpha,beta)<<std::endl;
			//std::cout<<H(s[mc[0]*N_m + cc[0]],s[mc[1]*N_m + cc[1]])<<std::endl;
			//std::cout<<w[0]<<" "<<w[1]<<" "<<w[0]*w[1]<<std::endl;
			//std::cout<<H(s[mc[0]*N_m + cc[0]],s[mc[1]*N_m + cc[1]]) * w[0]*w[1]<<std::endl;
			return H(s[mc[0]*N_m + cc[0]],s[mc[1]*N_m + cc[1]]) * w[0]*w[1];
		}
	} else {
		w[0] = 0.0;
		w[1] = 0.0;
		for(unsigned int i(0);i<N_m;i++){
			w[0] += Ainv[mc[0]](cc[0],i)*A[mc[1]](i,cc[1]);
			w[1] += Ainv[mc[1]](cc[1],i)*A[mc[0]](i,cc[0]);
		}
		return w[0]*w[1];
	}
}

template<typename T>
void System<T>::print(){
	std::cout<<"{";
	for(unsigned int i(0); i<N_spin; i++){
		std::cout<<"{ ";
		for(unsigned int j(0); j<N_m; j++){
			std::cout<<s[i*N_m+j]<<" ";
		}
		std::cout<<"}";
	}
	std::cout<<"}"<<std::endl;
	//for(unsigned int i(0); i<N_site; i++){
		//std::cout<<wis[i]<<" ";
	//}
	//std::cout<<std::endl;
	for(unsigned int i(0);i<2;i++){
		std::cout<<mc[i]<<" "<< cc[i]<<std::endl;
	}
	for(unsigned int i(0);i<N_spin;i++){
		std::cout<<std::endl;
		A[i].print();
	}
	//for(unsigned int i(0);i<N_spin;i++){
		//std::cout<<std::endl;
		//(A[i]*Ainv[i]).print();
	//}
}
/*}*/

/*methods that modify the class*/
/*{*/
template<typename T>
void System<T>::init_state(Matrice<T> const& EVec){
	unsigned int site(0);
	unsigned int N_as(N_site);
	unsigned int available_sites[N_site];
	for(unsigned int i(0); i < N_as; i++){
		available_sites[i]  = i;	
	}

	for(unsigned int i(0); i < N_spin; i++){
		A[i] = Matrice<T> (N_m);
		Ainv[i] = Matrice<T> (N_m);
		for(unsigned int j(0); j < N_m; j++){
			site = rand() % N_as;
			s[j+i*N_m] = available_sites[site];
			wis[available_sites[site]] = i*N_m+j; 	
			for(unsigned int k(0); k < N_m; k++){
				A[i](k,j) = EVec(available_sites[site],k);
			}
			for(unsigned int k(site); k < N_as-1; k++){
				available_sites[k] = available_sites[k+1];
			}
			N_as--;
		}
		Ainv[i] = A[i];
		Lapack<T> A_(Ainv[i].ptr(),Ainv[i].size(),'G');
		A_.inv();
	}
	//for(unsigned int i(0);i<N_spin;i++){
		//(Ainv[i]*A[i]).print();
	//}
}

template<typename T>
void System<T>::swap() {
	mc[0] = rand() % N_spin;
	mc[1] = rand() % N_spin;
	while(mc[0]==mc[1]){
		mc[1] = rand() % N_spin;
	}
	cc[0] = rand() % N_m;
	cc[1] = rand() % N_m;
}

template<typename T>
void System<T>::swap(unsigned int a, unsigned int b) {
	mc[0] = wis[a] / N_m;
	mc[1] = wis[b] / N_m;
	cc[0] = wis[a] % N_m;
	cc[1] = wis[b] % N_m;
}

template<typename T>
void System<T>::update(){
	unsigned int a(0),b(0);
	a = mc[0]*N_m+cc[0];
	b = mc[1]*N_m+cc[1];

	unsigned int s_tmp(0);
	s_tmp = s[b];
	s[b] = s[a];
	s[a] = s_tmp;

	a = s[a];
	b = s[b];

	s_tmp = wis[b];
	wis[b] = wis[a];
	wis[a] = s_tmp;

	T tmp(0.0);
	for(unsigned int i(0); i<N_m; i++){
		tmp = A[mc[0]](i,cc[0]);
		A[mc[0]](i,cc[0]) = A[mc[1]](i,cc[1]);
		A[mc[1]](i,cc[1]) = tmp;
	}
	
	T tmp_start(0.0);
	for(unsigned int m(0);m<2;m++){
		for(unsigned int i(0);i<N_m;i++){
			if(cc[m] == i){ tmp_start = -1.0; }
			else { tmp_start = 0.0;}
			for(unsigned int j(0);j<N_m;j++){
				tmp = tmp_start;
				for(unsigned int k(0);k<N_m;k++){
					tmp += Ainv[mc[m]](i,k)*A[mc[m]](k,cc[m]);
				}
				tmp_mat(i,j) = tmp*Ainv[mc[m]](cc[m],j)/w[m];
			}
		}
		Ainv[mc[m]] -= tmp_mat;
	}
}
/*}*/
#endif
