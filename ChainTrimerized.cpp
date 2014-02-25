#include"ChainTrimerized.hpp"

ChainTrimerized::ChainTrimerized(Container const& param):
	Chain<double>(param,"chain-trimerized")
{
	std::cerr<<"Ensure the correct number of site in function of the unit cell"<<std::endl;
}

ChainTrimerized::~ChainTrimerized(){}

void ChainTrimerized::create(Vector<double> const& t){
	a_ = t.size();
	t_= t;
	compute_T();
	diagonalize_T('S');
	if(degenerate_){
		bc_ = -1;
		T_.set(n_,n_,0.0);
		compute_T();
		diagonalize_T('S');
	}
	if(degenerate_){
		bc_ = 1;
		T_.set(n_,n_,0.0);
		compute_T();
		diagonalize_T('S');
		filename_ += "-degenerate";
	}
	for(unsigned int spin(0);spin<N_;spin++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_;j++){
				EVec_(i+spin*n_,j) = T_(i,j);
			}
		}
	}
}

void ChainTrimerized::compute_T(){
	Vector<unsigned int> neighbourg;
	for(unsigned int i(0); i < n_; i += a_){
		for(unsigned int j(0); j<a_; j++){
			neighbourg = get_neighbourg(i+j);
			T_(i+j,neighbourg(0)) = t_(j);
		}
	}
	T_(n_-1,0) = bc_*t_(a_-1);
	T_ += T_.transpose();
}

void ChainTrimerized::compute_P(Matrix<double>& P){
	P.set(n_,n_);
	P(n_ -1,0) = bc_;
	for(unsigned int i(0); i< n_-1; i++){
		P(i,i+1) = 1.0;
	}
}

void ChainTrimerized::save(){
	Write w(filename_+".jdbin");
	RST rst;
	rst.text("Spin chain, with different hopping term on the odd and even sites");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("ref (wave function)",ref_);
	w("n (particles' number)",n_);
	w("N (N of SU(N))",N_);
	w("m (particles per site' number)",m_);
	w("M (particles' number of each color)",M_);
	w("sts (connected sites)",sts_);
	w("t (hopping terms)",t_);
	w("EVec (unitary matrix)",EVec_);
	w("bc (boundary condition)",bc_);
}

void ChainTrimerized::get_param(Container& param){
	GenericSystem<double>get_param(param);
	param.set("t",t_);
}
