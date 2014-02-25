#ifndef DEF_GENERICSYSTEM
#define DEF_GENERICSYSTEM

#include "Lapack.hpp"
#include "Gnuplot.hpp"
#include "PSTricks.hpp"
#include "Container.hpp"

/*!Class that creates a file containing all the necessary information to run a
 * Monte-Carlo simulation.
*/
template<typename Type>
class GenericSystem{
	public:
		/*!Parseur needs N and m, z is the coordination number*/
		GenericSystem(Container const& param, unsigned int z, std::string filename); 
		/*Simple destructor*/
		virtual ~GenericSystem();

		virtual void save()=0;
		virtual void create()=0;

		virtual void get_input(Container& input);
		virtual void get_param(Container& param);

	protected:
		Vector<unsigned int> const ref_;//!< type of wavefunction
		unsigned int const n_;		//!< sites' number
		unsigned int const N_;		//!< colors' number
		unsigned int const m_;		//!< particles per site's number
		unsigned int const M_;		//!< particles' number of each color
		unsigned int const z_;		//!< coordination number
		double bc_;					//!< boundary condition
		Matrix<unsigned int> sts_;	//!< list of connected sites
		Matrix<Type> T_;			//!< Gutzwiller Hamiltonian
		Matrix<Type> EVec_;			//!< eigenvectors Matrix (transfer Matrix)
		bool degenerate_;			//!< no degeneracy at the fermi level
		std::string filename_;		//!< filename of the System

		/*!return the neighbours of site i*/
		virtual Vector<unsigned int> get_neighbourg(unsigned int i)=0;
		/*!compute the array of pairs of swapping sites*/
		void compute_sts();

		/*!compute the eigenvectors from the mean field Hamiltonian*/
		void diagonalize_T(char mat_type);
};

template<typename Type>
GenericSystem<Type>::GenericSystem(Container const& param, unsigned int z, std::string filename): 
	ref_(param.get<Vector<unsigned int> >("ref")),
	n_(param.get<unsigned int>("n")),
	N_(param.get<unsigned int>("N")), 
	m_(param.get<unsigned int>("m")),
	M_((m_*n_)/N_), 
	z_(z),
	bc_(1),
	sts_(n_*z/2,2),
	T_(n_,n_,0.0),
	EVec_(N_*n_,M_),
	degenerate_(false),
	filename_(filename)
{ 
	if(M_*N_ != m_*n_){ std::cerr<<"GenericSystem::GenericSystem(param,z,filename) : There is not an equal number of color"<<std::endl; }
	if(m_ > N_){ std::cerr<<"GenericSystem::GenericSystem(param,z,filename) : m>N is impossible"<<std::endl;} 
	filename_ += "-N" + tostring(N_);
	filename_ += "-m" + tostring(m_);
	filename_ += "-S" + tostring(n_);
}

template<typename Type>
GenericSystem<Type>::~GenericSystem(){}

template<typename Type>
void GenericSystem<Type>::compute_sts(){
	unsigned int k(0);
	Vector<unsigned int> neighbourg;
	for(unsigned int i(0); i<n_;i++){
		neighbourg = get_neighbourg(i);
		for(unsigned int j(0); j<z_/2; j++){
			sts_(k,0) = i;
			sts_(k,1) = neighbourg(j);
			k++;
		}
	}
}

template<typename Type>
void GenericSystem<Type>::diagonalize_T(char mat_type){
	Lapack<Type> ES(&T_,false, mat_type);
	Vector<double> EVal;
	ES.eigensystem(&EVal,true);
	if(std::abs(EVal(M_) - EVal(M_-1))<1e-12){ degenerate_ = true; }
}

template<typename Type>
void GenericSystem<Type>::get_param(Container& param){
	param.set("n",n_);
	param.set("N",N_);
	param.set("m",m_);
	param.set("bc",bc_);
}

template<typename Type>
void GenericSystem<Type>::get_input(Container& input){
	input.set("n",n_);
	input.set("N",N_);
	input.set("m",m_);
	input.set("sts",sts_);
	input.set("EVec",EVec_);
}
#endif
