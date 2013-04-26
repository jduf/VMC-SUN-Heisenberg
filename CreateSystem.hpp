#ifndef DEF_CREATESYSTEM
#define DEF_CREATESYSTEM

#include "Write.hpp"
#include "Array2D.hpp"
#include "Matrice.hpp"
#include "Lapack.hpp"
#include "RST.hpp"

#include <complex>
#include <sstream>
#include <cmath>

/*!Class that creates a file containing all the necessary information to run a
 * Monte-Carlo simulation.
 *
 *  
 *
*/
template<typename M>
class CreateSystem{
	public:
		CreateSystem(unsigned int N_m, unsigned int N_spin, unsigned int N_n, double parity);
		~CreateSystem();

	private:
		unsigned int const N_m, N_n, N_spin, N_site;//!<
		Array2D<unsigned int> sts;//!< 
		Matrice<double> H;//!< hopping matrix
		Matrice<M> T;//!< eigenvectors matrix (transfer matrix)
		std::string filename; //!<
		bool is_complex;//!<
		bool successful;//!<
		char mat_type;//!<
		RST rst;//!< will be added before the values in the header
		Write* w;//!<

		/*!Compute the hopping matrix for a chain*/
		void compute_H();
		/*!Compute the hopping matrix for a square lattice*/
		void compute_H(unsigned int N_row, unsigned int N_col, double parity);
		/*!Compute the eigenvectors from the mean field hamiltonian*/
		void compute_EVec();
		void compute_sts();

		void save();
		void set_size(unsigned int& N_row, unsigned int& N_col);
};

template<typename M>
CreateSystem<M>::CreateSystem(unsigned int N_m, unsigned int N_spin, unsigned int N_n, double parity):
	N_m(N_m),
	N_n(N_n),
	N_spin(N_spin), 
	N_site(N_m*N_spin),
	sts(N_spin*N_m*N_n/2,2),
	H(N_spin*N_m),
	T(N_spin*N_m),
	filename(""),
	mat_type('U'),
	is_complex(false),
	successful(false),
	rst(),
	w(NULL)
{
	std::stringstream ss1;
	std::stringstream ss2;
	ss1<<N_spin;
	ss2<<N_spin*N_m;
	filename = "-N"+ss1.str() + "-S" + ss2.str();
	switch(N_n){
		case 2:
			{
				filename="chain"+filename;
				compute_H();
				compute_EVec();
				successful = true;
				rst.text("Chain");
				save();
				break;
			}
		case 4:
			{
				filename="square"+filename;
				unsigned int N_row(floor(sqrt(N_site)));
				unsigned int N_col(floor(sqrt(N_site)));
				if(N_site==N_row*N_col){
					compute_H(N_row,N_col,parity);
				} else {
					std::cerr<<"CreateSystem : the cluster is not a square"<<std::endl;
				}
				compute_EVec();
				if(successful){
					if(is_complex){ rst.text("Chiral spin liquid"); }
					else{ rst.text("Fermi sea"); }
					std::stringstream ss3;
					std::stringstream ss4;
					ss3<<N_row;
					ss4<<N_col;
					if(parity == -1){ filename += "-" + ss3.str()+"x"+ ss4.str() + "-P" ;}
					else { filename += "-" + ss3.str()+"x"+ ss4.str() + "-A" ;}
					save();
					(*w)("N_row",N_row);
					(*w)("N_col",N_col);
					(*w)("parity",parity);
				} else {
					if(parity == -1){
						std::cerr<<"CreateSystem : degeneate for PBC"<<std::endl;
					} else {
						std::cerr<<"CreateSystem : degeneate for APBC"<<std::endl;
					}
				}
				break;
			}
		default:
			{
				std::cerr<<"CreateSystem : lattice type undefined"<<std::endl;
			}
	}
}

template<typename M>
CreateSystem<M>::~CreateSystem(){
	if(w){ delete w; }
}

template<typename M>
void CreateSystem<M>::compute_sts(){
	unsigned int k(0);
	for(unsigned int i(0); i<N_site;i++){
		for(unsigned int j(i+1); j<N_site;j++){
			if ( std::abs(H(i,j)) > 1e-4){
				sts(k,0) = i;
				sts(k,1) = j;
				k++;
			}
		}
	}
}

template<typename M>
void CreateSystem<M>::compute_EVec(){
	successful = false;
	Lapack<M> ES(T.ptr(),N_site, mat_type);
	Vecteur<double> EVal(N_site);
	ES.eigensystem(EVal);
	//std::cout<<EVal<<std::endl;
	if(std::abs(EVal(N_m) - EVal(N_m-1))>1e-10){ successful = true; }
}

template<typename M>
void CreateSystem<M>::save(){
	if(successful){
		if(is_complex){ filename += "-c"; }
		else{filename += "-r"; }
		compute_sts();

		w = new Write(filename+".jdbin");
		rst.title("Input values","~");
		w->header(rst.get());
		(*w)("complex",is_complex);
		(*w)("N_spin",N_spin);
		(*w)("N_m",N_m);
		(*w)("N_n",N_n);
		(*w)("sts",sts);
		(*w)("H",H);
		(*w)("T",T);
	}
}

template<typename M>
void CreateSystem<M>::set_size(unsigned int& N_row, unsigned int& N_col){
	N_col = ceil(sqrt(N_site))+N_spin;
	bool test(true);
	while( test ){
		test = false;
		N_col--; 
		if(N_col % N_spin != 0){ test = true; }
		if(N_site % N_col != 0){ test = true; }
		if(N_col < 3) { test = true; }
	}
	N_row = N_site/N_col;
	if(N_col % N_spin != 0){
		unsigned int tmp(N_col);
		N_col=N_row;
		N_row=tmp;
	}
}

#endif
