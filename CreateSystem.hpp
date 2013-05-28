#ifndef DEF_CREATESYSTEM
#define DEF_CREATESYSTEM

#include "Write.hpp"
#include "Read.hpp"
#include "Array2D.hpp"
#include "Matrice.hpp"
#include "Lapack.hpp"
#include "RST.hpp"

#include <complex>
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
		CreateSystem(unsigned int N_m, unsigned int N_spin, unsigned int N_n, double bc);
		CreateSystem(unsigned int N_m, unsigned int N_spin, unsigned int N_n, double td, std::string hopfile);
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
		Write* w;//!<output file that contains all the informations of the system

		/*!Compute the hopping matrix for a chain*/
		void compute_T();
		/*!Compute the hopping matrix for a square or honeycomb lattice*/
		void compute_T(unsigned int N_row, unsigned int N_col, double bc);
		/*!Compute the eigenvectors from the mean field hamiltonian*/
		void compute_EVec();
		/*!Compute the array of pairs of swapping sites*/
		void compute_sts();

		void save();
};

template<typename M>
CreateSystem<M>::CreateSystem(unsigned int N_m, unsigned int N_spin, unsigned int N_n, double bc):
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
	filename = "-N"+tostring(N_spin) + "-S" + tostring(N_spin*N_m);
	switch(N_n){
		case 2:
			{
				filename="chain"+filename;
				if(N_m % 2 == 0){ 
					filename += "-A";
					bc = -1;
				} else {
					filename += "-P";
					bc = 1;
				}
				H(0, N_site -1 ) = -1.0;
				T(0, N_site -1 ) = -bc;
				compute_T();
				compute_sts();
				compute_EVec();
				successful = true;
				rst.text("Spin chain");
				save();
				(*w)("bc",bc);
				break;
			}
		case 3:
			{
				filename="honeycomb"+filename;
				unsigned int N_row(6);
				unsigned int N_col(3);
				assert(N_row*N_col*4==N_site);
				compute_T(N_row,N_col,bc);
				compute_sts();
				compute_EVec();
				if(successful){
					if(bc == 1){
						std::cout<<"system with PBC created"<<std::endl;
						filename += "-" + tostring(N_row) +"x"+ tostring(N_col) + "-P" ;
					} else {
						std::cout<<"system with APBC created"<<std::endl;
						filename += "-" + tostring(N_row) +"x"+ tostring(N_col) + "-A" ;
					}
					save();
					(*w)("N_col",N_col);
					(*w)("N_row",N_row);
					(*w)("bc",bc);
				} else {
					if(bc == 1){
						std::cerr<<"CreateSystem : degeneate for PBC"<<std::endl;
					} else {
						std::cerr<<"CreateSystem : degeneate for APBC"<<std::endl;
					}
				}
				break;
			}
		case 4:
			{
				filename="square"+filename;
				unsigned int N_row(floor(sqrt(N_site)));
				unsigned int N_col(floor(sqrt(N_site)));
				if(N_site==N_row*N_col){
					compute_T(N_row,N_col,bc);
					compute_sts();
					compute_EVec();
					if(successful){
						if(is_complex){ rst.text("Chiral spin liquid"); }
						else{ rst.text("Fermi sea"); }
						if(bc == 1){
							std::cout<<"system with PBC created"<<std::endl;
							filename += "-" + tostring(N_row) +"x"+ tostring(N_col) + "-P" ;
						} else {
							std::cout<<"system with APBC created"<<std::endl;
							filename += "-" + tostring(N_row) +"x"+ tostring(N_col) + "-A" ;
						}
						save();
						(*w)("N_row",N_row);
						(*w)("N_col",N_col);
						(*w)("bc",bc);
					} else {
						if(bc == 1){
							std::cerr<<"CreateSystem : degeneate for PBC"<<std::endl;
						} else {
							std::cerr<<"CreateSystem : degeneate for APBC"<<std::endl;
						}
					}
				} else {
					std::cerr<<"CreateSystem : the cluster is not a square"<<std::endl;
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
CreateSystem<M>::CreateSystem(unsigned int N_m, unsigned int N_spin, unsigned int N_n, double td, std::string hopfile):
	N_m(N_m),
	N_n(N_n),
	N_spin(N_spin), 
	N_site(N_m*N_spin),
	sts(N_spin*N_m*N_n/2,2),
	H(N_spin*N_m,0.0),
	T(N_spin*N_m,0.0),
	filename(""),
	mat_type('U'),
	is_complex(false),
	successful(false),
	rst(),
	w(NULL)
{
	std::cerr<<"System created with "<<hopfile<<std::endl;
	filename = "-N" + tostring(N_spin) + "-S" + tostring(N_spin*N_m) + "-td" + tostring(td);
	switch(N_n){
		case 3:
			{
				filename="honeycomb"+filename;
				Read r(hopfile);
				Matrice<double> tmp(N_spin*N_m);
				r>>tmp;
				double t(-1.0);
				double th(-1.0);
				double bc(-1.0);
				for(unsigned int i(0);i<N_site;i++){
					for(unsigned int j(0);j<N_site;j++){
						if(std::abs(tmp(i,j)+1.0) < 1e-4){//original file th=-1
							H(i,j) = t;
							T(i,j) = th;
						}
						if(std::abs(tmp(i,j)-2.0) < 1e-4){ //original file td=2
							H(i,j) = t;
							T(i,j) = td;
						}
						if(std::abs(tmp(i,j)-1.0) < 1e-4){
							H(i,j) = t;
							T(i,j) = bc*th;
						}
						if(std::abs(tmp(i,j)+2.0) < 1e-4){
							H(i,j) = t;
							T(i,j) = bc*td;
						}
					}
				}
				mat_type = 'S';
				compute_sts();
				//for(unsigned int i(0);i<sts.row();i++){
					//std::cout<<sts(i,0)+1<<" "<<sts(i,1)+1<<" "<<H(sts(i,0),sts(i,1))<<std::endl;
				//}
				compute_EVec();
				if(successful){
					save();
					(*w)("th",th);
					(*w)("td",td);
					(*w)("bc",bc);
				} else {
					std::cerr<<"CreateSystem : degeneate"<<std::endl;
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
	//Matrice<M> A(T);
	Lapack<M> ES(T.ptr(),N_site, mat_type);
	Vecteur<double> EVal(N_site);
	ES.eigensystem(EVal);
	//std::cout<<EVal<<std::endl;
	//Matrice<M> Tinv(T);
	//Lapack<M> Tinv_(Tinv.ptr(),N_site, 'G');
	//Tinv_.inv();
	//Matrice<M> vp(Tinv*A*T);
	//std::cout<<vp.diag()<<std::endl;
	//A.print_mathematica();
	
	if(std::abs(EVal(N_m) - EVal(N_m-1))>1e-10){ successful = true; }
}

template<typename M>
void CreateSystem<M>::save(){
	if(successful){
		if(is_complex){ filename += "-c"; }
		else{filename += "-r"; }

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
#endif
