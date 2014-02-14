#include"ChainTrimerized.hpp"

ChainTrimerized::ChainTrimerized(Container const& param):
	Chain<double>(param,"chain-trimerized")
{
	std::cerr<<"Ensure the correct number of site in function of the unit cell"<<std::endl;
}

ChainTrimerized::~ChainTrimerized(){}

void ChainTrimerized::compute_T(Vector<double> const& t){
	unsigned int npuc(t.size()+1);
	t_.set(npuc);
	t_(0) = 1.0;
	for(unsigned int i(1); i<npuc;i++){
		t_(i) = t(i-1);
	}
	T_.set(n_,n_,0);
	for(unsigned int i(0); i < n_-npuc; i += npuc){
		for(unsigned int j(0); j<npuc; j++){
			T_(i+j,i+j+1) = t_(j);
		}
	}
	for(unsigned int i(0); i< npuc-1; i++){
		T_(n_-npuc+i,n_-npuc+i+1) = bc_*t_(i);
	}
	T_(n_-1,0) = t_(npuc-1);
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

void ChainTrimerized::properties(Container& c){
	c.set("ref",ref_);
	c.set("n",n_);
	c.set("N",N_);
	c.set("m",m_);
	c.set("M",M_);
	c.set("sts",sts_);
	c.set("t",t_);
	c.set("EVec",EVec_);
	c.set("bc",bc_);
}

void ChainTrimerized::create(Vector<double> const& t){
	compute_T(t);
	diagonalize_T('S');
	for(unsigned int spin(0);spin<N_;spin++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_;j++){
				EVec_(i+spin*n_,j) = T_(i,j);
			}
		}
	}
}
