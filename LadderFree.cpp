#include "LadderFree.hpp"

LadderFree::LadderFree(System const& s, Vector<double> const& t, Vector<double> const& mu):
	System(s),
	Ladder<double>(set_spuc(t,mu,N_/m_),"ladder-free"),
	t_(t),
	mu_(mu)
{
	if(status_==2 && t_.ptr()){
		init_fermionic();

		system_info_.text("LadderFree :");
		system_info_.text(" Each color has the same Hamiltonian.");

		filename_ += "-t";
		for(unsigned int i(0);i<t_.size();i++){
			filename_ += ((t_(i)>=0)?"+":"")+my::tostring(t_(i));
		}
		filename_ += "-mu";
		for(unsigned int i(0);i<mu_.size();i++){
			filename_ += ((mu_(i)>=0)?"+":"")+my::tostring(mu_(i));
		}
	}
}

/*{method needed for running*/
void LadderFree::compute_H(){
	H_.set(n_,n_,0);

	unsigned int s0(0);
	unsigned int s1(0);
	unsigned int k(0);
	unsigned int l(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		if(obs_[0](i,4)){ H_(s0,s1) = mu_(0)*t_(k); }
		else { H_(s0,s1) = t_(k); }
		if(!obs_[0](i,3)){
			H_(s0,s0) = mu_(l)/2.0; 
			l = (l+1)%spuc_;
		}
		k = (k+1)%spuc_;
	}
	H_ += H_.transpose();
}

void LadderFree::create(){
	compute_H();
	diagonalize(true);

	if(status_==1){
		for(unsigned int c(0);c<N_;c++){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
			}
		}
	}
}

void LadderFree::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("t=(");
		Vector<double> param(t_.size()+mu_.size());

		for(unsigned int i(0);i<t_.size()-1;i++){ 
			param(i) = t_(i); 
			s += my::tostring(t_(i))+",";
		}
		param(t_.size()-1) = t_.back(); 
		s += my::tostring(t_.back())+") "+RST::math("\\mu")+"=(";

		for(unsigned int i(0);i<mu_.size()-1;i++){
			param(i+t_.size()) = mu_(i); 
			s += my::tostring(mu_(i))+",";
		}
		param.back() = mu_.back(); 
		s += my::tostring(mu_.back())+")";

		w.add_header()->title(s,'<');
		w<<param;
		w.add_header()->add(system_info_.get());
	} else { w<<t_<<" "<<mu_<<" "; }
}

unsigned int LadderFree::set_spuc(Vector<double> const& t, Vector<double> const& mu, unsigned int const& spuc){
	if((t.size()*2/3)%spuc == 0 && mu.size()%spuc==0 && mu.size()<9){ return mu.size(); }
	else { 
		std::cerr<<__PRETTY_FUNCTION__<<" : invalid or incoherent t and mu sizes : t:="<<t.size()<<", mu:="<<mu.size()<<std::endl;
		return 1; 
	}
}

void LadderFree::get_wf_symmetries(std::vector<Matrix<int> >& sym) const {
	int A(J_(0)>J_(1)?0:-1); //link between sites 0-1 (rung)
	int B(J_(0)>J_(1)?-1:0); //link between sites 0-2 (leg)
	switch(spuc_){
		case 2:
			{
				Matrix<int> tmp(1,3,1);
				tmp(0,1) = B;
				sym.push_back(tmp);
			}break;
		case 4:
			{
				Matrix<int> tmp;
				/*{no symmetry breaking*/
				/*0,0*/
				sym.push_back(tmp);
				/*{ 1 pi-flux*/
				tmp.set(1,3);
				/*pi,0*/
				tmp(0,0) = 1;
				tmp(0,1) = 1;
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*0,pi*/
				tmp(0,0) = 4;
				tmp(0,1) = 4;
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*}*/

				/*{ 2 pi-flux*/
				tmp.set(2,3);
				/*pi,pi*/
				tmp(0,0) = 1;
				tmp(0,1) = 1;
				tmp(0,2) = -1;
				tmp(1,0) = 4;
				tmp(1,1) = 4;
				tmp(1,2) = -1;
				sym.push_back(tmp);
				/*}*/

				/*{facing dimerization*/
				tmp.set(3,3);
				tmp(0,0) = 1;
				tmp(0,1) = B;
				tmp(0,2) = 1;

				tmp(1,0) = 3;
				tmp(1,1) = A;
				tmp(1,2) = 1;

				tmp(2,0) = 4;
				tmp(2,1) = 2;
				tmp(2,2) = 1;
				/*0,0*/
				sym.push_back(tmp);
				/*pi,0*/
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*0,pi*/
				tmp(0,2) = 1;
				tmp(2,2) = -1;
				sym.push_back(tmp);
				/*pi,pi*/
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*}*/
				/*}*/

				/*{shifted dimerization*/
				tmp(0,0) = 1;
				tmp(0,1) = 2;
				tmp(0,2) = 1;

				tmp(1,0) = 4;
				tmp(1,1) = B;
				tmp(1,2) = 1;

				tmp(2,0) = 3;
				tmp(2,1) = A;
				tmp(2,2) = 1;
				/*0,0*/
				sym.push_back(tmp);
				/*pi,0*/
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*0,pi*/
				tmp(0,2) = 1;
				tmp(1,2) = -1;
				sym.push_back(tmp);
				/*0,pi*/
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*}*/
			}break;
		case 6:
			{
				Matrix<int> tmp;
				/*{no symmetry breaking*/
				/*0,0,0*/
				sym.push_back(tmp);
				/*{ 1 pi-flux*/
				tmp.set(1,3);
				/*pi,0,0*/
				tmp(0,0) = 1;
				tmp(0,1) = 1;
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*0,pi,0*/
				tmp(0,0) = 4;
				tmp(0,1) = 4;
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*0,0,pi*/
				tmp(0,0) = 7;
				tmp(0,1) = 7;
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*}*/

				/*{ 2 pi-flux*/
				tmp.set(2,3);
				/*pi,pi,0*/
				tmp(0,0) = 1;
				tmp(0,1) = 1;
				tmp(0,2) = -1;
				tmp(1,0) = 4;
				tmp(1,1) = 4;
				tmp(1,2) = -1;
				sym.push_back(tmp);
				/*pi,0,pi*/
				tmp(0,0) = 1;
				tmp(0,1) = 1;
				tmp(0,2) = -1;
				tmp(1,0) = 7;
				tmp(1,1) = 7;
				tmp(1,2) = -1;
				sym.push_back(tmp);
				/*0,pi,pi*/
				tmp(0,0) = 4;
				tmp(0,1) = 4;
				tmp(0,2) = -1;
				tmp(1,0) = 7;
				tmp(1,1) = 7;
				tmp(1,2) = -1;
				sym.push_back(tmp);
				/*}*/

				/*{ 3 pi-flux*/
				tmp.set(3,3);
				/*pi,pi,pi*/
				tmp(0,0) = 1;
				tmp(0,1) = 1;
				tmp(0,2) = -1;
				tmp(1,0) = 4;
				tmp(1,1) = 4;
				tmp(1,2) = -1;
				tmp(2,0) = 7;
				tmp(2,1) = 7;
				tmp(2,2) = -1;
				sym.push_back(tmp);
				/*}*/
				/*}*/

				/*{facing trimerization*/
				tmp.set(6,3);
				tmp(0,0) = 1;
				tmp(0,1) = B;
				tmp(0,2) = 1;

				tmp(1,0) = 2;
				tmp(1,1) = B;
				tmp(1,2) = 1;

				tmp(2,0) = 3;
				tmp(2,1) = A;
				tmp(2,2) = 1;

				tmp(3,0) = 4;
				tmp(3,1) = B;
				tmp(3,2) = 1;

				tmp(4,0) = 6;
				tmp(4,1) = A;
				tmp(4,2) = 1;

				tmp(5,0) = 7;
				tmp(5,1) = 5;
				tmp(5,2) = 1;
				/*0,0,0*/
				sym.push_back(tmp);
				/*pi,0,0*/
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*pi,pi,0*/
				tmp(3,2) = -1;
				sym.push_back(tmp);
				/*0,pi,0*/
				tmp(0,2) = 1;
				sym.push_back(tmp);
				/*0,pi,pi*/
				tmp(5,2) = -1;
				sym.push_back(tmp);
				/*0,0,pi*/
				tmp(3,2) = 1;
				sym.push_back(tmp);
				/*pi,0,pi*/
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*pi,pi,pi*/
				tmp(3,2) = -1;
				sym.push_back(tmp);
				/*}*/

				/*{shifted trimerization*/
				tmp(0,0) = 1;
				tmp(0,1) = 2;
				tmp(0,2) = 1;

				tmp(1,0) = 3;
				tmp(1,1) = A;
				tmp(1,2) = 1;

				tmp(2,0) = 4;
				tmp(2,1) = 2;
				tmp(2,2) = 1;

				tmp(3,0) = 5;
				tmp(3,1) = 2;
				tmp(3,2) = 1;

				tmp(4,0) = 6;
				tmp(4,1) = A;
				tmp(4,2) = 1;

				tmp(5,0) = 7;
				tmp(5,1) = B;
				tmp(5,2) = 1;
				/*0,0,0*/
				sym.push_back(tmp);
				/*pi,0,0*/
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*pi,pi,0*/
				tmp(2,2) = -1;
				sym.push_back(tmp);
				/*0,pi,0*/
				tmp(0,2) = 1;
				sym.push_back(tmp);
				/*0,pi,pi*/
				tmp(5,2) = -1;
				sym.push_back(tmp);
				/*0,0,pi*/
				tmp(2,2) = 1;
				sym.push_back(tmp);
				/*pi,0,pi*/
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*pi,pi,pi*/
				tmp(2,2) = -1;
				sym.push_back(tmp);
				/*}*/

				/*{dimer and box*/
				tmp.set(4,3);
				tmp(0,0) = 1;
				tmp(0,1) = B;
				tmp(0,2) = 1;

				tmp(1,0) = 3;
				tmp(1,1) = 6;
				tmp(1,2) = 1;

				tmp(2,0) = 4;
				tmp(2,1) = 2;
				tmp(2,2) = 1;

				tmp(3,0) = 7;
				tmp(3,1) = 5;
				tmp(3,2) = 1;
				/*0,0,0*/
				sym.push_back(tmp);
				/*pi,0,0*/
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*pi,pi,0*/
				tmp(2,2) = -1;
				sym.push_back(tmp);
				/*0,pi,0*/
				tmp(0,2) = 1;
				sym.push_back(tmp);
				/*0,pi,pi*/
				tmp(3,2) = -1;
				sym.push_back(tmp);
				/*0,0,pi*/
				tmp(2,2) = 1;
				sym.push_back(tmp);
				/*pi,0,pi*/
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*pi,pi,pi*/
				tmp(2,2) = -1;
				sym.push_back(tmp);
				/*}*/
			}break;
			//case 8:
			//{
			//Matrix<int> tmp;
			///*{no symmetry breaking*/
			///*0,0,0,0*/
			//sym.push_back(tmp);
			///*{ 1 pi-flux*/
			//tmp.set(1,3);
			///*pi,0,0,0*/
			//tmp(0,0) = 1;
			//tmp(0,1) = 1;
			//tmp(0,2) = -1;
			//sym.push_back(tmp);
			///*0,pi,0,0*/
			//tmp(0,0) = 4;
			//tmp(0,1) = 4;
			//tmp(0,2) = -1;
			//sym.push_back(tmp);
			///*0,0,pi,0*/
			//tmp(0,0) = 7;
			//tmp(0,1) = 7;
			//tmp(0,2) = -1;
			//sym.push_back(tmp);
			///*0,0,0,pi*/
			//tmp(0,0) = 10;
			//tmp(0,1) = 10;
			//tmp(0,2) = -1;
			//sym.push_back(tmp);
			///*}*/
			//
			///*{ 2 pi-flux*/
			//tmp.set(2,3);
			///*pi,pi,0,0*/
			//tmp(0,0) = 1;
			//tmp(0,1) = 1;
			//tmp(0,2) = -1;
			//tmp(1,0) = 4;
			//tmp(1,1) = 4;
			//tmp(1,2) = -1;
			//sym.push_back(tmp);
			///*pi,0,pi,0*/
			//tmp(0,0) = 1;
			//tmp(0,1) = 1;
			//tmp(0,2) = -1;
			//tmp(1,0) = 7;
			//tmp(1,1) = 7;
			//tmp(1,2) = -1;
			//sym.push_back(tmp);
			///*0,pi,0,pi*/
			//tmp(0,0) = 4;
			//tmp(0,1) = 4;
			//tmp(0,2) = -1;
			//tmp(1,0) = 10;
			//tmp(1,1) = 10;
			//tmp(1,2) = -1;
			//sym.push_back(tmp);
			///*0,pi,pi,0*/
			//tmp(0,0) = 4;
			//tmp(0,1) = 4;
			//tmp(0,2) = -1;
			//tmp(1,0) = 7;
			//tmp(1,1) = 7;
			//tmp(1,2) = -1;
			//sym.push_back(tmp);
			///*0,0,pi,pi*/
			//tmp(0,0) = 7;
			//tmp(0,1) = 7;
			//tmp(0,2) = -1;
			//tmp(1,0) = 10;
			//tmp(1,1) = 10;
			//tmp(1,2) = -1;
			//sym.push_back(tmp);
			///*pi,0,0,pi*/
			//tmp(0,0) = 1;
			//tmp(0,1) = 1;
			//tmp(0,2) = -1;
			//tmp(1,0) = 10;
			//tmp(1,1) = 10;
			//tmp(1,2) = -1;
			//sym.push_back(tmp);
			///*}*/
			//
			///*{ 3 pi-flux*/
			//tmp.set(3,3);
			///*pi,pi,pi,0*/
			//tmp(0,0) = 1;
			//tmp(0,1) = 1;
			//tmp(0,2) = -1;
			//tmp(1,0) = 4;
			//tmp(1,1) = 4;
			//tmp(1,2) = -1;
			//tmp(2,0) = 7;
			//tmp(2,1) = 7;
			//tmp(2,2) = -1;
			//sym.push_back(tmp);
			///*pi,pi,0,pi*/
			//tmp(0,0) = 1;
			//tmp(0,1) = 1;
			//tmp(0,2) = -1;
			//tmp(1,0) = 4;
			//tmp(1,1) = 4;
			//tmp(1,2) = -1;
			//tmp(2,0) = 10;
			//tmp(2,1) = 10;
			//tmp(2,2) = -1;
			//sym.push_back(tmp);
			///*pi,0,pi,pi*/
			//tmp(0,0) = 1;
			//tmp(0,1) = 1;
			//tmp(0,2) = -1;
			//tmp(1,0) = 7;
			//tmp(1,1) = 7;
			//tmp(1,2) = -1;
			//tmp(2,0) = 10;
			//tmp(2,1) = 10;
			//tmp(2,2) = -1;
			//sym.push_back(tmp);
			///*0,pi,pi,pi*/
			//tmp(0,0) = 4;
			//tmp(0,1) = 4;
			//tmp(0,2) = -1;
			//tmp(1,0) = 7;
			//tmp(1,1) = 7;
			//tmp(1,2) = -1;
			//tmp(2,0) = 10;
			//tmp(2,1) = 10;
			//tmp(2,2) = -1;
			//sym.push_back(tmp);
			///*}*/
			//
			///*{ 4 pi-flux*/
			//tmp.set(4,3);
			///*pi,pi,pi,pi*/
			//tmp(0,0) = 1;
			//tmp(0,1) = 1;
			//tmp(0,2) = -1;
			//tmp(1,0) = 4;
			//tmp(1,1) = 4;
			//tmp(1,2) = -1;
			//tmp(2,0) = 7;
			//tmp(2,1) = 7;
			//tmp(2,2) = -1;
			//tmp(3,0) = 10;
			//tmp(3,1) = 10;
			//tmp(3,2) = -1;
			//sym.push_back(tmp);
			///*}*/
			///*}*/
			//
			///*{facing tetramerization*/
			//tmp.set(7,3);
			//tmp(0,0) = 1;
			//tmp(0,1) = B;
			//tmp(0,2) = 1;
			//
			//tmp(1,0) = 4;
			//tmp(1,1) = 2;
			//tmp(1,2) = 1;
			//
			//tmp(2,0) = 5;
			//tmp(2,1) = B;
			//tmp(2,2) = 1;
			//
			//tmp(3,0) = 6;
			//tmp(3,1) = 3;
			//tmp(3,2) = 1;
			//
			//tmp(4,0) = 7;
			//tmp(4,1) = B;
			//tmp(4,2) = 1;
			//
			//tmp(5,0) = 9;
			//tmp(5,1) = A;
			//tmp(5,2) = 1;
			//
			//tmp(6,0) = 10;
			//tmp(6,1) = 8;
			//tmp(6,2) = 1;
			///*0,0,0,0*/
			//sym.push_back(tmp);
			///*pi,0,0,0*/
			//tmp(0,2) = -1;
			//sym.push_back(tmp);
			///*pi,pi,0,0*/
			//tmp(1,2) = -1;
			//sym.push_back(tmp);
			///*0,pi,0,0*/
			//tmp(0,2) = 1;
			//sym.push_back(tmp);
			///*0,pi,pi,0*/
			//tmp(4,2) = -1;
			//sym.push_back(tmp);
			///*0,0,pi,0*/
			//tmp(1,2) = 1;
			//sym.push_back(tmp);
			///*0,0,pi,pi*/
			//tmp(6,2) = -1;
			//sym.push_back(tmp);
			///*0,0,0,pi*/
			//tmp(4,2) = 1;
			//sym.push_back(tmp);
			///*0,pi,0,pi*/
			//tmp(1,2) = -1;
			//sym.push_back(tmp);
			///*pi,0,pi,0*/
			//tmp(0,2) = -1;
			//tmp(1,2) = 1;
			//tmp(4,2) = -1;
			//tmp(6,2) = 1;
			//sym.push_back(tmp);
			///*pi,0,0,pi*/
			//tmp(4,2) = 1;
			//tmp(6,2) = -1;
			//sym.push_back(tmp);
			///*pi,pi,pi,0*/
			//tmp(1,2) = -1;
			//tmp(4,2) = -1;
			//tmp(6,2) = 1;
			//sym.push_back(tmp);
			///*pi,pi,0,pi*/
			//tmp(4,2) = 1;
			//tmp(6,2) = -1;
			//sym.push_back(tmp);
			///*pi,0,pi,pi*/
			//tmp(1,2) = 1;
			//tmp(4,2) = -1;
			//sym.push_back(tmp);
			///*0,pi,pi,pi*/
			//tmp(0,2) = 1;
			//tmp(1,2) = -1;
			//sym.push_back(tmp);
			///*pi,pi,pi,pi*/
			//tmp(0,2) = -1;
			//sym.push_back(tmp);
			///*}*/
			//
			///*{shifted by one tetramerization*/
			//tmp.set(6,3);
			//tmp(0,0) = 1;
			//tmp(0,1) = 2;
			//tmp(0,2) = 1;
			//
			//tmp(1,0) = 4;
			//tmp(1,1) = 5;
			//tmp(1,2) = 1;
			//
			//tmp(2,0) = 7;
			//tmp(2,1) = 2;
			//tmp(2,2) = 1;
			//
			//tmp(3,0) = 8;
			//tmp(3,1) = 2;
			//tmp(3,2) = 1;
			//
			//tmp(4,0) = 9;
			//tmp(4,1) = 3;
			//tmp(4,2) = 1;
			//
			//tmp(5,0) = 10;
			//tmp(5,1) = B;
			//tmp(5,2) = 1;
			///*0,0,0,0*/
			//sym.push_back(tmp);
			///*pi,0,0,0*/
			//tmp(0,2) = -1;
			//sym.push_back(tmp);
			///*pi,pi,0,0*/
			//tmp(1,2) = -1;
			//sym.push_back(tmp);
			///*0,pi,0,0*/
			//tmp(0,2) = 1;
			//sym.push_back(tmp);
			///*0,pi,pi,0*/
			//tmp(2,2) = -1;
			//sym.push_back(tmp);
			///*0,0,pi,0*/
			//tmp(1,2) = 1;
			//sym.push_back(tmp);
			///*0,0,pi,pi*/
			//tmp(5,2) = -1;
			//sym.push_back(tmp);
			///*0,0,0,pi*/
			//tmp(2,2) = 1;
			//sym.push_back(tmp);
			///*0,pi,0,pi*/
			//tmp(1,2) = -1;
			//sym.push_back(tmp);
			///*pi,0,pi,0*/
			//tmp(0,2) = -1;
			//tmp(1,2) = 1;
			//tmp(2,2) = -1;
			//tmp(5,2) = 1;
			//sym.push_back(tmp);
			///*pi,pi,pi,0*/
			//tmp(1,2) = -1;
			//sym.push_back(tmp);
			///*pi,pi,0,pi*/
			//tmp(2,2) = 1;
			//tmp(5,2) = -1;
			//sym.push_back(tmp);
			///*pi,0,pi,pi*/
			//tmp(1,2) = 1;
			//tmp(2,2) = -1;
			//sym.push_back(tmp);
			///*0,pi,pi,pi*/
			//tmp(0,2) = 1;
			//tmp(1,2) = -1;
			//sym.push_back(tmp);
			///*0,pi,pi,pi*/
			//tmp(0,2) = -1;
			//sym.push_back(tmp);
			///*}*/
			//
			///*{shifted by two tetramerization*/
			//tmp.set(8,3);
			//tmp(0,0) = 1;
			//tmp(0,1) = B;
			//tmp(0,2) = 1;
			//
			//tmp(1,0) = 3;
			//tmp(1,1) = A;
			//tmp(1,2) = 1;
			//
			//tmp(2,0) = 4;
			//tmp(2,1) = 8;
			//tmp(2,2) = 1;
			//
			//tmp(3,0) = 5;
			//tmp(3,1) = B;
			//tmp(3,2) = 1;
			//
			//tmp(4,0) = 6;
			//tmp(4,1) = A;
			//tmp(4,2) = 1;
			//
			//tmp(5,0) = 7;
			//tmp(5,1) = B;
			//tmp(5,2) = 1;
			//
			//tmp(6,0) = 9;
			//tmp(6,1) = A;
			//tmp(6,2) = 1;
			//
			//tmp(7,0) = 10;
			//tmp(7,1) = 2;
			//tmp(7,2) = 1;
			///*0,0,0,0*/
			//sym.push_back(tmp);
			///*pi,0,0,0*/
			//tmp(0,2) = -1;
			//sym.push_back(tmp);
			///*pi,pi,0,0*/
			//tmp(2,2) = -1;
			//sym.push_back(tmp);
			///*0,pi,0,0*/
			//tmp(0,2) = 1;
			//sym.push_back(tmp);
			///*0,pi,pi,0*/
			//tmp(5,2) = -1;
			//sym.push_back(tmp);
			///*0,0,pi,0*/
			//tmp(2,2) = 1;
			//sym.push_back(tmp);
			///*0,0,pi,pi*/
			//tmp(7,2) = -1;
			//sym.push_back(tmp);
			///*0,0,0,pi*/
			//tmp(5,2) = 1;
			//sym.push_back(tmp);
			///*0,pi,0,pi*/
			//tmp(2,2) = -1;
			//sym.push_back(tmp);
			///*pi,0,pi,0*/
			//tmp(0,2) = -1;
			//tmp(2,2) = 1;
			//tmp(5,2) = -1;
			//tmp(7,2) = 1;
			//sym.push_back(tmp);
			///*pi,0,0,pi*/
			//tmp(5,2) = 1;
			//tmp(7,2) = -1;
			//sym.push_back(tmp);
			///*pi,pi,pi,0*/
			//tmp(2,2) = -1;
			//tmp(5,2) = -1;
			//tmp(7,2) = 1;
			//sym.push_back(tmp);
			///*pi,pi,0,pi*/
			//tmp(5,2) = 1;
			//tmp(7,2) = -1;
			//sym.push_back(tmp);
			///*pi,0,pi,pi*/
			//tmp(2,2) = 1;
			//tmp(5,2) = -1;
			//sym.push_back(tmp);
			///*0,pi,pi,pi*/
			//tmp(0,2) = 1;
			//tmp(2,2) = -1;
			//sym.push_back(tmp);
			///*pi,pi,pi,pi*/
			//tmp(0,2) = -1;
			//sym.push_back(tmp);
			///*}*/
			//
			///*{double dimerization*/
			//tmp.set(7,3);
			//tmp(0,0) = 1;
			//tmp(0,1) = B;
			//tmp(0,2) = 1;
			//
			//tmp(1,0) = 2;
			//tmp(1,1) = B;
			//tmp(1,2) = 1;
			//
			//tmp(2,0) = 4;
			//tmp(2,1) = B;
			//tmp(2,2) = 1;
			//
			//tmp(3,0) = 6;
			//tmp(3,1) = A;
			//tmp(3,2) = 1;
			//
			//tmp(4,0) = 7;
			//tmp(4,1) = 5;
			//tmp(4,2) = 1;
			//
			//tmp(5,0) = 8;
			//tmp(5,1) = 5;
			//tmp(5,2) = 1;
			//
			//tmp(6,0) = 10;
			//tmp(6,1) = 5;
			//tmp(6,2) = 1;
			///*0,0,0,0*/
			//sym.push_back(tmp);
			///*pi,0,0,0*/
			//tmp(0,2) = -1;
			//sym.push_back(tmp);
			///*pi,pi,0,0*/
			//tmp(2,2) = -1;
			//sym.push_back(tmp);
			///*0,pi,0,0*/
			//tmp(0,2) = 1;
			//sym.push_back(tmp);
			///*0,pi,pi,0*/
			//tmp(4,2) = -1;
			//sym.push_back(tmp);
			///*0,0,pi,0*/
			//tmp(2,2) = 1;
			//sym.push_back(tmp);
			///*0,0,pi,pi*/
			//tmp(6,2) = -1;
			//sym.push_back(tmp);
			///*0,0,0,pi*/
			//tmp(4,2) = 1;
			//sym.push_back(tmp);
			///*0,pi,0,pi*/
			//tmp(2,2) = -1;
			//sym.push_back(tmp);
			///*pi,0,pi,0*/
			//tmp(0,2) = -1;
			//tmp(2,2) = 1;
			//tmp(4,2) = -1;
			//tmp(6,2) = 1;
			//sym.push_back(tmp);
			///*pi,pi,pi,0*/
			//tmp(2,2) = -1;
			//sym.push_back(tmp);
			///*pi,pi,0,pi*/
			//tmp(4,2) = 1;
			//tmp(6,2) = -1;
			//sym.push_back(tmp);
			///*pi,0,pi,pi*/
			//tmp(2,2) = 1;
			//tmp(4,2) = -1;
			//sym.push_back(tmp);
			///*0,pi,pi,pi*/
			//tmp(0,2) = 1;
			//tmp(2,2) = -1;
			//sym.push_back(tmp);
			///*pi,pi,pi,pi*/
			//tmp(0,2) = -1;
			//sym.push_back(tmp);
			///*}*/
			//}break;
		case 8:
			{
				Matrix<int> tmp;
				/*{no symmetry breaking*/
				/*0,0,0,0*/
				sym.push_back(tmp);
				/*{ 1 pi-flux*/
				tmp.set(1,3);
				/*pi,0,0,0*/
				tmp(0,0) = 1;
				tmp(0,1) = 1;
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*0,pi,0,0*/
				tmp(0,0) = 4;
				tmp(0,1) = 4;
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*0,0,pi,0*/
				tmp(0,0) = 7;
				tmp(0,1) = 7;
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*0,0,0,pi*/
				tmp(0,0) = 10;
				tmp(0,1) = 10;
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*}*/

				/*{ 2 pi-flux*/
				tmp.set(2,3);
				/*pi,pi,0,0*/
				tmp(0,0) = 1;
				tmp(0,1) = 1;
				tmp(0,2) = -1;
				tmp(1,0) = 4;
				tmp(1,1) = 4;
				tmp(1,2) = -1;
				sym.push_back(tmp);
				/*pi,0,pi,0*/
				tmp(0,0) = 1;
				tmp(0,1) = 1;
				tmp(0,2) = -1;
				tmp(1,0) = 7;
				tmp(1,1) = 7;
				tmp(1,2) = -1;
				sym.push_back(tmp);
				/*0,pi,0,pi*/
				tmp(0,0) = 4;
				tmp(0,1) = 4;
				tmp(0,2) = -1;
				tmp(1,0) = 10;
				tmp(1,1) = 10;
				tmp(1,2) = -1;
				sym.push_back(tmp);
				/*0,pi,pi,0*/
				tmp(0,0) = 4;
				tmp(0,1) = 4;
				tmp(0,2) = -1;
				tmp(1,0) = 7;
				tmp(1,1) = 7;
				tmp(1,2) = -1;
				sym.push_back(tmp);
				/*0,0,pi,pi*/
				tmp(0,0) = 7;
				tmp(0,1) = 7;
				tmp(0,2) = -1;
				tmp(1,0) = 10;
				tmp(1,1) = 10;
				tmp(1,2) = -1;
				sym.push_back(tmp);
				/*pi,0,0,pi*/
				tmp(0,0) = 1;
				tmp(0,1) = 1;
				tmp(0,2) = -1;
				tmp(1,0) = 10;
				tmp(1,1) = 10;
				tmp(1,2) = -1;
				sym.push_back(tmp);
				/*}*/

				/*{ 3 pi-flux*/
				tmp.set(3,3);
				/*pi,pi,pi,0*/
				tmp(0,0) = 1;
				tmp(0,1) = 1;
				tmp(0,2) = -1;
				tmp(1,0) = 4;
				tmp(1,1) = 4;
				tmp(1,2) = -1;
				tmp(2,0) = 7;
				tmp(2,1) = 7;
				tmp(2,2) = -1;
				sym.push_back(tmp);
				/*pi,pi,0,pi*/
				tmp(0,0) = 1;
				tmp(0,1) = 1;
				tmp(0,2) = -1;
				tmp(1,0) = 4;
				tmp(1,1) = 4;
				tmp(1,2) = -1;
				tmp(2,0) = 10;
				tmp(2,1) = 10;
				tmp(2,2) = -1;
				sym.push_back(tmp);
				/*pi,0,pi,pi*/
				tmp(0,0) = 1;
				tmp(0,1) = 1;
				tmp(0,2) = -1;
				tmp(1,0) = 7;
				tmp(1,1) = 7;
				tmp(1,2) = -1;
				tmp(2,0) = 10;
				tmp(2,1) = 10;
				tmp(2,2) = -1;
				sym.push_back(tmp);
				/*0,pi,pi,pi*/
				tmp(0,0) = 4;
				tmp(0,1) = 4;
				tmp(0,2) = -1;
				tmp(1,0) = 7;
				tmp(1,1) = 7;
				tmp(1,2) = -1;
				tmp(2,0) = 10;
				tmp(2,1) = 10;
				tmp(2,2) = -1;
				sym.push_back(tmp);
				/*}*/

				/*{ 4 pi-flux*/
				tmp.set(4,3);
				/*pi,pi,pi,pi*/
				tmp(0,0) = 1;
				tmp(0,1) = 1;
				tmp(0,2) = -1;
				tmp(1,0) = 4;
				tmp(1,1) = 4;
				tmp(1,2) = -1;
				tmp(2,0) = 7;
				tmp(2,1) = 7;
				tmp(2,2) = -1;
				tmp(3,0) = 10;
				tmp(3,1) = 10;
				tmp(3,2) = -1;
				sym.push_back(tmp);
				/*}*/
				/*}*/

				/*{facing tetramerization*/
				tmp.set(13,3);
				tmp(0,0) = 1;
				tmp(0,1) = B;
				tmp(0,2) = 1;

				tmp(1,0) = 4;
				tmp(1,1) = 2;
				tmp(1,2) = 1;

				tmp(2,0) = 5;
				tmp(2,1) = B;
				tmp(2,2) = 1;

				tmp(3,0) = 6;
				tmp(3,1) = 3;
				tmp(3,2) = 1;

				tmp(4,0) = 7;
				tmp(4,1) = B;
				tmp(4,2) = 1;

				tmp(5,0) = 9;
				tmp(5,1) = A;
				tmp(5,2) = 1;

				tmp(6,0) = 10;
				tmp(6,1) = 8;
				tmp(6,2) = 1;

				tmp(7,0) = 11;
				tmp(7,1) = -1;
				tmp(7,2) = 0;

				tmp(8,0) = 13;
				tmp(8,1) = 12;
				tmp(8,2) = 1;

				tmp(9,0) = 14;
				tmp(9,1) = 12;
				tmp(9,2) = 1;

				tmp(10,0) = 15;
				tmp(10,1) = 12;
				tmp(10,2) = 1;

				tmp(11,0) = 16;
				tmp(11,1) = -1;
				tmp(11,2) = 0;

				tmp(12,0) = 17;
				tmp(12,1) = -1;
				tmp(12,2) = 0;
				/*0,0,0,0*/
				sym.push_back(tmp);
				/*pi,0,0,0*/
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*pi,pi,0,0*/
				tmp(1,2) = -1;
				sym.push_back(tmp);
				/*0,pi,0,0*/
				tmp(0,2) = 1;
				sym.push_back(tmp);
				/*0,pi,pi,0*/
				tmp(4,2) = -1;
				sym.push_back(tmp);
				/*0,0,pi,0*/
				tmp(1,2) = 1;
				sym.push_back(tmp);
				/*0,0,pi,pi*/
				tmp(6,2) = -1;
				sym.push_back(tmp);
				/*0,0,0,pi*/
				tmp(4,2) = 1;
				sym.push_back(tmp);
				/*0,pi,0,pi*/
				tmp(1,2) = -1;
				sym.push_back(tmp);
				/*pi,0,pi,0*/
				tmp(0,2) = -1;
				tmp(1,2) = 1;
				tmp(4,2) = -1;
				tmp(6,2) = 1;
				sym.push_back(tmp);
				/*pi,0,0,pi*/
				tmp(4,2) = 1;
				tmp(6,2) = -1;
				sym.push_back(tmp);
				/*pi,pi,pi,0*/
				tmp(1,2) = -1;
				tmp(4,2) = -1;
				tmp(6,2) = 1;
				sym.push_back(tmp);
				/*pi,pi,0,pi*/
				tmp(4,2) = 1;
				tmp(6,2) = -1;
				sym.push_back(tmp);
				/*pi,0,pi,pi*/
				tmp(1,2) = 1;
				tmp(4,2) = -1;
				sym.push_back(tmp);
				/*0,pi,pi,pi*/
				tmp(0,2) = 1;
				tmp(1,2) = -1;
				sym.push_back(tmp);
				/*pi,pi,pi,pi*/
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*}*/

				/*{shifted by two tetramerization*/
				tmp.set(14,3);
				tmp(0,0) = 1;
				tmp(0,1) = B;
				tmp(0,2) = 1;

				tmp(1,0) = 3;
				tmp(1,1) = A;
				tmp(1,2) = 1;

				tmp(2,0) = 4;
				tmp(2,1) = 8;
				tmp(2,2) = 1;

				tmp(3,0) = 5;
				tmp(3,1) = B;
				tmp(3,2) = 1;

				tmp(4,0) = 6;
				tmp(4,1) = A;
				tmp(4,2) = 1;

				tmp(5,0) = 7;
				tmp(5,1) = B;
				tmp(5,2) = 1;

				tmp(6,0) = 9;
				tmp(6,1) = A;
				tmp(6,2) = 1;

				tmp(7,0) = 10;
				tmp(7,1) = 2;
				tmp(7,2) = 1;

				tmp(8,0) = 11;
				tmp(8,1) = 12;
				tmp(8,2) = 1;

				tmp(9,0) = 13;
				tmp(9,1) = -1;
				tmp(9,2) = 0;

				tmp(10,0) = 14;
				tmp(10,1) = 12;
				tmp(10,2) = 1;

				tmp(11,0) = 15;
				tmp(11,1) = -1;
				tmp(11,2) = 0;

				tmp(12,0) = 16;
				tmp(12,1) = -1;
				tmp(12,2) = 0;

				tmp(13,0) = 17;
				tmp(13,1) = 12;
				tmp(13,2) = 1;
				/*0,0,0,0*/
				sym.push_back(tmp);
				/*pi,0,0,0*/
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*pi,pi,0,0*/
				tmp(2,2) = -1;
				sym.push_back(tmp);
				/*0,pi,0,0*/
				tmp(0,2) = 1;
				sym.push_back(tmp);
				/*0,pi,pi,0*/
				tmp(5,2) = -1;
				sym.push_back(tmp);
				/*0,0,pi,0*/
				tmp(2,2) = 1;
				sym.push_back(tmp);
				/*0,0,pi,pi*/
				tmp(7,2) = -1;
				sym.push_back(tmp);
				/*0,0,0,pi*/
				tmp(5,2) = 1;
				sym.push_back(tmp);
				/*0,pi,0,pi*/
				tmp(2,2) = -1;
				sym.push_back(tmp);
				/*pi,0,pi,0*/
				tmp(0,2) = -1;
				tmp(2,2) = 1;
				tmp(5,2) = -1;
				tmp(7,2) = 1;
				sym.push_back(tmp);
				/*pi,0,0,pi*/
				tmp(5,2) = 1;
				tmp(7,2) = -1;
				sym.push_back(tmp);
				/*pi,pi,pi,0*/
				tmp(2,2) = -1;
				tmp(5,2) = -1;
				tmp(7,2) = 1;
				sym.push_back(tmp);
				/*pi,pi,0,pi*/
				tmp(5,2) = 1;
				tmp(7,2) = -1;
				sym.push_back(tmp);
				/*pi,0,pi,pi*/
				tmp(2,2) = 1;
				tmp(5,2) = -1;
				sym.push_back(tmp);
				/*0,pi,pi,pi*/
				tmp(0,2) = 1;
				tmp(2,2) = -1;
				sym.push_back(tmp);
				/*pi,pi,pi,pi*/
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*}*/
			}break;
		default:{ std::cerr<<__PRETTY_FUNCTION__<<"unknown spuc_"<<std::endl; }
	}
}
/*}*/

/*{method needed for checking*/
void LadderFree::check(){
	//check_lattice();
	//compute_H();
	//plot_band_structure();
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="ladder-free";
	display_results();
}

void LadderFree::plot(bool const& create_image){
	if(obs_.size()>3){
		/*!long range correlations*/
		/*{*/
		unsigned int llr(obs_[1].nval());
		Vector<std::complex<double> > Ck_intra(llr,0.0);
		Vector<std::complex<double> > Ck_inter(llr,0.0);
		std::complex<double> normalize_intra(0.0);
		std::complex<double> normalize_inter(0.0);
		double dk(2.0*M_PI/llr);

		for(unsigned int k(0);k<llr;k++){
			for(unsigned int i(0);i<llr;i++){
				Ck_intra(k) += std::polar(obs_[1][i].get_x(),dk*k*i);
				Ck_inter(k) += std::polar(obs_[2][i].get_x(),dk*k*i);
			}
			normalize_intra += Ck_intra(k);
			normalize_inter += Ck_inter(k);
		}
		Ck_intra /= dk*normalize_intra;
		Ck_inter /= dk*normalize_inter;

		IOFiles file_c(analyse_+path_+dir_+filename_+"-lr-c.dat",true);
		IOFiles file_sf(analyse_+path_+dir_+filename_+"-lr-sf.dat",true);
		for(unsigned int l(0);l<llr;l++){
			file_c<<l<<" "<<obs_[1][l]<<" "<<obs_[2][l]<<IOFiles::endl;
			file_sf<<dk*l<<" "<<Ck_intra(l).real()<<" "<<Ck_intra(l).imag()<<" "<<Ck_inter(l).real()<<" "<<Ck_inter(l).imag()<<IOFiles::endl;
		}

		Gnuplot gp(analyse_+path_+dir_,filename_+"-lr");
		gp.multiplot();
		/*{correlations*/
		gp.range("x","0",llr/2);

		gp.tics("x");
		gp.margin("0.1","0.5","0.9","0.5");
		gp+="plot '"+filename_+"-lr-c.dat' u 1:2:3 w errorbars lt 1 lc 6 t 'intra'";

		gp.margin("0.1","0.5","0.5","0.1");
		gp.tics("x","");
		gp+="plot '"+filename_+"-lr-c.dat' u 1:6:7  w errorbars lt 1 lc 7 t 'inter'";
		/*}*/
		/*{structure factor*/
		gp.range("x","0","pi");
		//gp.range("y2","0","");

		gp.key("left");
		gp.tics("x");
		gp.tics("y");
		gp.tics("x2","('' pi/3, '' pi/2, '' 2*pi/3, '' pi) mirror");
		gp.tics("y2","mirror");
		gp.margin("0.5","0.9","0.9","0.5");
		gp+="plot '"+filename_+"-lr-sf.dat' u 1:2 axes x1y2 lt 1 lc 6 notitle,\\";
		gp+="     '"+filename_+"-lr-sf.dat' u 1:3 axes x1y2 lt 2 lc 6 notitle";

		gp.margin("0.5","0.9","0.5","0.1");
		gp.tics("x","('$\\pi/3$' pi/3, '$\\pi/2$' pi/2, '$2\\pi/3$' 2*pi/3, '$\\pi$' pi)");
		gp+="plot '"+filename_+"-lr-sf.dat' u 1:4 axes x1y2 lt 1 lc 7 notitle,\\";
		gp+="     '"+filename_+"-lr-sf.dat' u 1:5 axes x1y2 lt 2 lc 7 notitle";
		/*}*/
		gp.save_file();
		gp.create_image(create_image,true);
		/*}*/
	}
	if(obs_.size()==5){
		/*!(anti)symmetric correlations and structure factors*/
		/*{*/
		unsigned int llr(obs_[2].nval());
		Vector<std::complex<double> > Ckm(llr,0.0);
		Vector<std::complex<double> > Ckp(llr,0.0);
		std::complex<double> normalize_m(0.0);
		std::complex<double> normalize_p(0.0);
		double dk(2.0*M_PI/llr);

		for(unsigned int k(0);k<llr;k++){
			for(unsigned int i(0);i<llr;i++){
				Ckm(k) += std::polar(obs_[1][i].get_x()-obs_[2][i].get_x()-obs_[3][i].get_x()+obs_[4][i].get_x(),dk*k*i);
				Ckp(k) += std::polar(obs_[1][i].get_x()+obs_[2][i].get_x()+obs_[3][i].get_x()+obs_[4][i].get_x(),dk*k*i);
			}
			normalize_m += Ckm(k);
			normalize_p += Ckp(k);
		}
		Ckm /= dk*normalize_m;
		Ckp /= dk*normalize_p;

		IOFiles file_c(analyse_+path_+dir_+filename_+"-as-c.dat",true);
		IOFiles file_sf(analyse_+path_+dir_+filename_+"-as-sf.dat",true);
		for(unsigned int l(0);l<llr;l++){
			file_c<<l<<" "<<obs_[1][l]<<" "<<obs_[2][l]<<" "<<obs_[3][l]<<" "<<obs_[4][l]<<IOFiles::endl;
			file_sf<<dk*l<<" "<<Ckm(l).real()<<" "<<Ckm(l).imag()<<" "<<Ckp(l).real()<<" "<<Ckp(l).imag()<<IOFiles::endl;
		}

		Gnuplot gp(analyse_+path_+dir_,filename_+"-as");
		gp.multiplot();
		/*{correlations*/
		gp.range("x","0",llr/2);

		gp.tics("x");
		gp.margin("0.1","0.5","0.9","0.5");
		gp+="plot '"+filename_+"-as-c.dat' u 1:($2-$6-$10+$14):($3+$7+$18+$18) w errorbars lt 1 lc 6 t '$-$'";

		gp.margin("0.1","0.5","0.5","0.1");
		gp.tics("x","");
		gp+="plot '"+filename_+"-as-c.dat' u 1:($2+$6+$10+$14):($3+$7+$18+$18)  w errorbars lt 1 lc 7 t '$+$'";
		/*}*/
		/*{structure factor*/
		gp.range("x","0","pi");
		//gp.range("y2","0","");

		gp.key("left");
		gp.tics("x");
		gp.tics("y");
		gp.tics("x2","('' pi/3, '' pi/2, '' 2*pi/3, '' pi)");
		gp.tics("y2","");
		gp.margin("0.5","0.9","0.9","0.5");
		gp+="plot '"+filename_+"-as-sf.dat' u 1:2 axes x1y2 lt 1 lc 6 notitle,\\";
		gp+="     '"+filename_+"-as-sf.dat' u 1:3 axes x1y2 lt 2 lc 6 notitle";

		gp.margin("0.5","0.9","0.5","0.1");
		gp.tics("x","('$\\pi/3$' pi/3, '$\\pi/2$' pi/2, '$2\\pi/3$' 2*pi/3, '$\\pi$' pi)");
		gp+="plot '"+filename_+"-as-sf.dat' u 1:4 axes x1y2 lt 1 lc 7 notitle,\\";
		gp+="     '"+filename_+"-as-sf.dat' u 1:5 axes x1y2 lt 2 lc 7 notitle";
		/*}*/
		gp.save_file();
		gp.create_image(create_image,true);
		/*}*/
	}
}

void LadderFree::lattice(){
	compute_H();

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	double x_shift(spuc_/2+2);

	PSTricks ps(info_+path_+dir_,filename_+"-pstricks");
	ps.begin(-1,-5,n_/1.5,2,filename_+"-pstricks");
	double t;
	double mu;
	double corr;
	unsigned int s0;
	unsigned int s1;
	std::string str;
	for(unsigned int i(0);i<3*spuc_/2;i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		xy0(0) = s0/2;
		xy0(1) = s0%2;
		xy1(0) = s1/2;
		xy1(1) = s1%2;

		t=H_(s0,s1);
		if(std::abs(t)>1e-4){
			if(xy1(0)<xy0(0)){
				xy1(0) = xy0(0)+1;
				linestyle="dashed";
			} else { linestyle="solid"; }
			if(t>0){ color = "blue"; }
			else   { color = "red"; }
			linewidth = my::tostring(std::abs(t))+"mm";
			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1),"linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}

		mu = H_(s0,s0);
		if(std::abs(mu)>1e-4){
			if(mu<0){ color = "magenta"; }
			else { color = "cyan"; }
			ps.circle(xy0,std::abs(mu),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
		}

		ps.put(xy0(0),xy0(1),"\\tiny{"+my::tostring(mu)+"}");
		switch(i%3){
			case 0: { ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+0.2,"\\tiny{"+my::tostring(t)+"}"); }break;
			case 1: { ps.put(xy0(0)+0.2,(xy0(1)+xy1(1))/2.0,"\\tiny{"+my::tostring(t)+"}"); }break;
			case 2: { ps.put((xy0(0)+xy1(0))/2.0,xy0(1)-0.2,"\\tiny{"+my::tostring(t)+"}"); }break;
		}

		if(obs_.size()>1){/*bound energy*/
			corr = obs_[1][i].get_x();
			if(std::abs(corr)>1e-4){
				if(corr>0){ color = "blue"; }
				else      { color = "red"; }
				linewidth = my::tostring(std::abs(corr))+"mm";

				ps.line("-",xy0(0)+x_shift,xy0(1),xy1(0)+x_shift,xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
			}

			if(i%3!=1){
				if(i%3==0){
					ps.put(xy0(0)+x_shift,xy0(1)-0.2,"\\tiny{"+my::tostring(s0)+"}");
				}
				str = my::tostring(corr);
				ps.put((xy0(0)+xy1(0))/2.0+2*x_shift,xy0(1),"\\tiny{"+str.substr(0,8)+"}");
				str = my::tostring(obs_[1][i].get_dx());
				if(obs_[1][i].get_dx()<1e-4){
					ps.put((xy0(0)+xy1(0))/2.0+2*x_shift,xy0(1)-0.2,"\\tiny{"+str.substr(0,4)+"e-"+str.substr(str.size()-2,2)+"}");
				} else {
					ps.put((xy0(0)+xy1(0))/2.0+2*x_shift,xy0(1)-0.2,"\\tiny{"+str.substr(0,8)+"}");
				}
			} else {
				ps.put(xy1(0)+x_shift,xy1(1)+0.2,"\\tiny{"+my::tostring(s1)+"}");
				str = my::tostring(corr);
				ps.put((xy0(0)+xy1(0))/2.0+2*x_shift,(xy0(1)+xy1(1))/2.0,"\\tiny{"+str.substr(0,8)+"}");
				str = my::tostring(obs_[1][i].get_dx());
				if(obs_[1][i].get_dx()<1e-4){
					ps.put((xy0(0)+xy1(0))/2.0+2*x_shift,(xy0(1)+xy1(1))/2.0-0.2,"\\tiny{"+str.substr(0,4)+"e-"+str.substr(str.size()-2,2)+"}");
				} else {
					ps.put((xy0(0)+xy1(0))/2.0+2*x_shift,(xy0(1)+xy1(1))/2.0-0.2,"\\tiny{"+str.substr(0,8)+"}");
				}
			}
		}
	}
	if(obs_.size()==5){/*long range correlations*/
		double rescale(0.75/obs_[2][0].get_x());
		unsigned int n(std::min(4*N_/m_+1,obs_[2].nval()));
		unsigned int idx;
		for(unsigned int i(0);i<n;i++){
			idx = (obs_[2].nval()-n/2+i)%obs_[1].nval();

			corr = obs_[2][idx].get_x()*rescale;
			xy0(0) = i;
			xy0(1) = -2;
			if(std::abs(corr)>1e-4){
				if(i!=idx){
					if(corr>0){ color = "blue"; }
					else      { color = "red"; }
				} else { color = "black"; }
				ps.circle(xy0,std::abs(corr),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
			}

			corr = obs_[3][idx].get_x()*rescale;
			xy0(1) = -1;
			if(std::abs(corr)>1e-4){
				if(corr>0){ color = "blue"; }
				else      { color = "red"; }
				ps.circle(xy0,std::abs(corr),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
			}
		}
	}
	ps.end(true,true,true);
}

void LadderFree::display_results(){
	lattice();
	plot(true);
	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title(RST::math("\\theta=")+my::tostring(acos(J_(0))) + " : t=(");
		std::string run_cmd("./mc -s:wf ladder-free");
		run_cmd += " -u:N " + my::tostring(N_);
		run_cmd += " -u:m " + my::tostring(m_);
		run_cmd += " -u:n " + my::tostring(n_);
		run_cmd += " -i:bc "+ my::tostring(bc_);
		run_cmd += " -d:theta " + my::tostring(acos(J_(0)));
		run_cmd += " -d:t ";
		for(unsigned int i(0);i<t_.size()-1;i++){
			title   += my::tostring(t_(i)) + ","; 
			run_cmd += my::tostring(t_(i)) + ","; 
		}
		title   += my::tostring(t_.back()) + "), "+RST::math("\\mu")+"=(";
		run_cmd += my::tostring(t_.back()) + " -d:mu ";
		for(unsigned int i(0);i<mu_.size()-1;i++){
			title   += my::tostring(mu_(i)) + ","; 
			run_cmd += my::tostring(mu_(i)) + ","; 
		}
		title   += my::tostring(mu_.back()) + ")";
		run_cmd += my::tostring(mu_.back()) + " -d -u:tmax 10";
		if(dir_ == "P/" || dir_ == "O/" || dir_ == "A/"){
			rst_file_->title("|theta"+my::tostring(acos(J_(0)))+"|_",'-');
			rst_file_->replace("theta"+my::tostring(acos(J_(0))),title);
		} else { rst_file_->title(title,'-'); }

		rst_file_->change_text_onclick("run command",run_cmd);

		rst_file_->figure(dir_+filename_+"-pstricks.png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+"-pstricks.pdf")+RST::scale("200"));
		if(obs_.size()>1){
			rst_file_->figure(relative_path+filename_+"-lr.png","long range correlations",RST::target(relative_path+filename_+"-lr.gp")+RST::scale("200"));
		} 
		if(obs_.size()==5){
			rst_file_->figure(relative_path+filename_+"-as.png","(anti)symmetric correlations",RST::target(relative_path+filename_+"-as.gp")+RST::scale("200"));
		}
	}
}
/*}*/

/*{method needed for analysing*/
std::string LadderFree::extract_level_6(){
	(*data_write_)<<N_<<" "<<m_<<" "<<n_<<" "<<bc_<<" "<<asin(J_(1))<<" "<<obs_[0][0]<<IOFiles::endl;

	display_results();

	save_param(*jd_write_);
	save(*jd_write_);

	return filename_;
}
/*}*/
