#include "LadderFree.hpp"

LadderFree::LadderFree(System const& s, Vector<double> const& t):
	System(s),
	Ladder<double>(set_spuc(t),"ladderfree"),
	t_(t)
{
	if(status_==2 && t_.ptr()){
		init_fermionic();
		system_info_.text("LadderFree : all colors experience the same Hamiltonian");
		filename_ += "-t";
		for(unsigned int j(0);j<t_.size();j++){
			filename_ += ((t_(j)>0)?"+":"")+my::tostring(t_(j));
		}
	}
}

/*{method needed for running*/
void LadderFree::compute_H(){
	H_.set(n_,n_,0);
	Matrix<int> nb;
	unsigned int k(0);

	for(unsigned int i(0);i<n_;i++){
		nb = get_neighbourg(i);
		if(i%spuc_){
			if(i%2){
				H_(i,nb(0,0)) = nb(0,1)*t_(k++);
			} else {
				H_(i,nb(0,0)) = nb(0,1)*t_(k++);
				H_(i,nb(1,0)) = nb(1,1)*t_(k++);
			}
		} else {
			/*!decoupled chains in the limit J⊥(0)=J_(1)->0, therefore, in
			 * that case the variational parameter is t⊥ (t‖ is set to 1),
			 * otherwise the inverse is done*/
			if(J_(0)>J_(1)){
				H_(i,nb(0,0)) = nb(0,1);
				H_(i,nb(1,0)) = nb(1,1)*t_(k++);
			} else {
				H_(i,nb(0,0)) = nb(0,1)*t_(k++);
				H_(i,nb(1,0)) = nb(1,1);
			}
		}
		k = k%t_.size();
	}
	H_ += H_.transpose();
}

void LadderFree::create(unsigned int const& which_observables){
	compute_H();
	diagonalize(true);
	(void)(which_observables);

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
	std::string t_string("");
	for(unsigned int i(0);i<t_.size()-1;i++){
		t_string += my::tostring(t_(i))+",";
	}
	t_string += my::tostring(t_.back());
	w.add_header()->title("t=("+t_string+")",'<');
	w<<t_;
	GenericSystem<double>::save_param(w);
}

unsigned int LadderFree::set_spuc(Vector<double> const& t){
	switch(t.size()){
		case 0: { return 1; } break;//!to allow a silent construction in case of undefined t
		case 2: { return 2; } break;
		case 5: { return 4; } break;
		case 8: { return 6; } break;
		case 11:{ return 8; } break;
		default:{
					std::cerr<<__PRETTY_FUNCTION__<<" : invalid t size : "<<t.size()<<std::endl;
					return 1;
				}
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
				tmp.set(7,3);
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

				/*{shifted by one tetramerization*/
				tmp.set(6,3);
				tmp(0,0) = 1;
				tmp(0,1) = 2;
				tmp(0,2) = 1;

				tmp(1,0) = 4;
				tmp(1,1) = 5;
				tmp(1,2) = 1;

				tmp(2,0) = 7;
				tmp(2,1) = 2;
				tmp(2,2) = 1;

				tmp(3,0) = 8;
				tmp(3,1) = 2;
				tmp(3,2) = 1;

				tmp(4,0) = 9;
				tmp(4,1) = 3;
				tmp(4,2) = 1;

				tmp(5,0) = 10;
				tmp(5,1) = B;
				tmp(5,2) = 1;
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
				tmp(2,2) = -1;
				sym.push_back(tmp);
				/*0,0,pi,0*/
				tmp(1,2) = 1;
				sym.push_back(tmp);
				/*0,0,pi,pi*/
				tmp(5,2) = -1;
				sym.push_back(tmp);
				/*0,0,0,pi*/
				tmp(2,2) = 1;
				sym.push_back(tmp);
				/*0,pi,0,pi*/
				tmp(1,2) = -1;
				sym.push_back(tmp);
				/*pi,0,pi,0*/
				tmp(0,2) = -1;
				tmp(1,2) = 1;
				tmp(2,2) = -1;
				tmp(5,2) = 1;
				sym.push_back(tmp);
				/*pi,pi,pi,0*/
				tmp(1,2) = -1;
				sym.push_back(tmp);
				/*pi,pi,0,pi*/
				tmp(2,2) = 1;
				tmp(5,2) = -1;
				sym.push_back(tmp);
				/*pi,0,pi,pi*/
				tmp(1,2) = 1;
				tmp(2,2) = -1;
				sym.push_back(tmp);
				/*0,pi,pi,pi*/
				tmp(0,2) = 1;
				tmp(1,2) = -1;
				sym.push_back(tmp);
				/*0,pi,pi,pi*/
				tmp(0,2) = -1;
				sym.push_back(tmp);
				/*}*/

				/*{shifted by two tetramerization*/
				tmp.set(8,3);
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

				/*{double dimerization*/
				tmp.set(7,3);
				tmp(0,0) = 1;
				tmp(0,1) = B;
				tmp(0,2) = 1;

				tmp(1,0) = 2;
				tmp(1,1) = B;
				tmp(1,2) = 1;

				tmp(2,0) = 4;
				tmp(2,1) = B;
				tmp(2,2) = 1;

				tmp(3,0) = 6;
				tmp(3,1) = A;
				tmp(3,2) = 1;

				tmp(4,0) = 7;
				tmp(4,1) = 5;
				tmp(4,2) = 1;

				tmp(5,0) = 8;
				tmp(5,1) = 5;
				tmp(5,2) = 1;

				tmp(6,0) = 10;
				tmp(6,1) = 5;
				tmp(6,2) = 1;
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
				tmp(4,2) = -1;
				sym.push_back(tmp);
				/*0,0,pi,0*/
				tmp(2,2) = 1;
				sym.push_back(tmp);
				/*0,0,pi,pi*/
				tmp(6,2) = -1;
				sym.push_back(tmp);
				/*0,0,0,pi*/
				tmp(4,2) = 1;
				sym.push_back(tmp);
				/*0,pi,0,pi*/
				tmp(2,2) = -1;
				sym.push_back(tmp);
				/*pi,0,pi,0*/
				tmp(0,2) = -1;
				tmp(2,2) = 1;
				tmp(4,2) = -1;
				tmp(6,2) = 1;
				sym.push_back(tmp);
				/*pi,pi,pi,0*/
				tmp(2,2) = -1;
				sym.push_back(tmp);
				/*pi,pi,0,pi*/
				tmp(4,2) = 1;
				tmp(6,2) = -1;
				sym.push_back(tmp);
				/*pi,0,pi,pi*/
				tmp(2,2) = 1;
				tmp(4,2) = -1;
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
		default:{  std::cerr<<__PRETTY_FUNCTION__<<" bla"<<std::endl; }
	}
}
/*}*/

/*{method needed for checking*/
void LadderFree::check(){
	lattice("./","lattice");
}

void LadderFree::lattice(std::string const& path, std::string const& filename){
	compute_H();
	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	double x_shift(spuc_/2+2);
	double y_shift(2);

	PSTricks ps(path,filename);
	ps.begin(-1,-5,n_/2,2,filename);
	double t;
	double corr;
	unsigned int s0;
	unsigned int s1;
	for(unsigned int i(0);i<spuc_+2;i++){
		s0 = link_types_[0](i,0);
		s1 = link_types_[0](i,1);
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

			if(t<0){ color = "red"; }
			else { color = "blue"; }
			linewidth = my::tostring(std::abs(t))+"mm";

			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1),"linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
			switch(i%3){
				case 0: { ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+0.2,"\\tiny{"+my::tostring(t)+"}"); }break;
				case 1: { ps.put(xy0(0)+0.2,(xy0(1)+xy1(1))/2.0,"\\tiny{"+my::tostring(t)+"}"); }break;
				case 2: { ps.put((xy0(0)+xy1(0))/2.0,xy0(1)-0.2,"\\tiny{"+my::tostring(t)+"}"); }break;
			}
		}

		if(i<corr_types_[0].size()){
			corr = corr_types_[0][i].get_x();
			if(std::abs(corr)>1e-4){
				if(corr<0){ color = "red"; }
				else { color = "blue"; }
				linewidth = my::tostring(std::abs(corr))+"mm";

				ps.line("-",xy0(0)+x_shift,xy0(1),xy1(0)+x_shift,xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
			}
		}

		if(i%3==0){ 
			ps.put(xy0(0),xy0(1)-0.2,"\\tiny{"+my::tostring(s0)+"}"); 
			ps.put(xy0(0)+x_shift,xy0(1)-0.2,"\\tiny{"+my::tostring(s0)+"}"); 
		}
		if(i%3==1){ 
			ps.put(xy1(0),xy1(1)+0.2,"\\tiny{"+my::tostring(s1)+"}"); 
			ps.put(xy1(0)+x_shift,xy1(1)+0.2,"\\tiny{"+my::tostring(s1)+"}"); 
		}
	}
	double rescale(corr_types_[1].size()?0.75/corr_types_[1][0].get_x():0);
	for(unsigned int i(0);i<corr_types_[0].size();i++){
		corr = corr_types_[1][i].get_x()*rescale;
		if(std::abs(corr)>1e-4){
			s0 = link_types_[1](i,0);
			s1 = link_types_[1](i,1);
			xy0(0) = s0/2;
			xy0(1) = s0%2-y_shift;
			xy1(0) = s1&2;
			xy1(1) = s1%2-y_shift;

			if(i){
				if(corr<0){ color = "red"; }
				else { color = "blue"; }
			} else { color = "black"; }

			ps.circle(xy0,std::abs(corr),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
		}
	}

	ps.end(true,true,true);
}
/*}*/

/*{method needed for analysing*/
std::string LadderFree::extract_level_6(){
	(*data_write_)<<N_<<" "<<m_<<" "<<n_<<" "<<bc_<<" "<<asin(J_(1))<<" "<<E_<<IOFiles::endl;

	lattice(info_+path_+dir_,filename_);
	rst_file_->figure(dir_+filename_+".png",RST::math("\\theta="+my::tostring(asin(J_(1))))+" : "+RST::math("E="+my::tostring(E_.get_x())+"\\pm"+my::tostring(E_.get_dx())),RST::target(dir_+filename_+".pdf")+RST::scale("200")); 

	save_param(*jd_write_);
	save_input(*jd_write_);
	save_output(*jd_write_);

	return filename_;
}
/*}*/
