#include"ChainPolymerized.hpp"

ChainPolymerized::ChainPolymerized(Container const& param):
	Chain<double>(param,"chain-Polymerized")
{}

ChainPolymerized::~ChainPolymerized(){}

void ChainPolymerized::create(double delta){
	delta_=delta;
	compute_T();
	diagonalize_T('S');
	for(unsigned int spin(0);spin<N_;spin++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_;j++){
				EVec_(i+spin*n_,j) = T_(i,j);
			}
		}
	}
}

void ChainPolymerized::compute_T(){
	Vector<unsigned int> neighbourg;
	double t(1.0);
	for(unsigned int i(0); i < n_; i += a_){
		for(unsigned int j(0); j<a_-1; j++){
			neighbourg = get_neighbourg(i+j);
			T_(i+j,neighbourg(0)) = t+delta_;
		}
		neighbourg = get_neighbourg(i+a_-1);
		T_(i+a_-1,neighbourg(0)) = t-delta_;
	}
	T_(n_-1,0) = bc_*(t-delta_);
	T_ += T_.transpose();
}

void ChainPolymerized::compute_P(Matrix<double>& P){
	P.set(n_,n_);
	P(n_ -1,0) = bc_;
	for(unsigned int i(0); i< n_-1; i++){
		P(i,i+1) = 1.0;
	}
}

void ChainPolymerized::save(){
	filename_ += "-delta" + tostring(delta_);
	filename_ += "-a" + tostring(a_);

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
	w("delta (t+-delta)",delta_);
	w("a (unit vector)",a_);
	w("EVec (unitary matrix)",EVec_);
	w("bc (boundary condition)",bc_);
}

void ChainPolymerized::get_param(Container& param){
	GenericSystem<double>::get_param(param);
	param.set("delta",delta_);
}
