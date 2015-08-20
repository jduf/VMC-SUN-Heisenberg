#include "LadderFree.hpp"

LadderFree::LadderFree(System const& s, Vector<double> const& t):
	System(s,3),
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
			if(J_(0)>J_(1)){
				/*!decoupled chains in the limit J⊥=J_(1)->0, therefore, in
				 * that case the variational parameter is t⊥ (t‖ is set to 1)*/
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
		case 0: { return 1; } break;
		case 2: { return 2; } break;
		case 5: { return 4; } break;
		case 8: { return 6; } break;
		case 11:{ return 8; } break;
		default:{
					std::cerr<<"unsigned int LadderFree::set_spuc(Vector<double> const& t) : invalid t size : "<<t.size()<<std::endl;
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
				/*no symmetry breaking*/
				/*0,0*/
				sym.push_back(tmp);
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

				tmp.set(3,3);
				/*facing dimerization*/
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

				/*shifted dimerization*/
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
			}break;
		case 6:
			{
				Matrix<int> tmp;
				/*no symmetry breaking*/
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

				tmp.set(6,3);
				/*{facing trimerization*/
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
				/*}*/
			}break;
		case 8:
			{
				Matrix<int> tmp;
				/*no symmetry breaking*/
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
				tmp(1,0) = 7;
				tmp(1,1) = 7;
				tmp(1,2) = -1;
				/*pi,pi,0,pi*/
				tmp(0,0) = 1;
				tmp(0,1) = 1;
				tmp(0,2) = -1;
				tmp(1,0) = 4;
				tmp(1,1) = 4;
				tmp(1,2) = -1;
				tmp(1,0) = 10;
				tmp(1,1) = 10;
				tmp(1,2) = -1;
				/*pi,0,pi,pi*/
				tmp(0,0) = 1;
				tmp(0,1) = 1;
				tmp(0,2) = -1;
				tmp(1,0) = 7;
				tmp(1,1) = 7;
				tmp(1,2) = -1;
				tmp(1,0) = 10;
				tmp(1,1) = 10;
				tmp(1,2) = -1;
				/*0,pi,pi,pi*/
				tmp(0,0) = 4;
				tmp(0,1) = 4;
				tmp(0,2) = -1;
				tmp(1,0) = 7;
				tmp(1,1) = 7;
				tmp(1,2) = -1;
				tmp(1,0) = 10;
				tmp(1,1) = 10;
				tmp(1,2) = -1;
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
				/*0,0,0*/
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
				/*pi,pi,pi,0*/
				tmp(1,2) = -1;
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
				/*pi,pi,pi,0*/
				tmp(2,2) = -1;
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
				/*}*/
			}break;
		default:{  std::cerr<<"void LadderFree:: bla"<<std::endl; }
	}
}
/*}*/

/*{method needed for checking*/
void LadderFree::check(){
	//compute_H();
	//std::cout<<"liens :"<<std::endl;
	//std::cout<<links_<<std::endl;
//
	//for(unsigned int i(0);i<n_;i++){
		//std::cout << "get_neighbourg. right - top/bot - left" << std::endl;
		//std::cout<<"i="<<i<<std::endl;
		//std::cout<<get_neighbourg(i)<<std::endl;// shows the links
	//}
	//std::cout<<"Hamiltonien"<<std::endl;
	//std::cout<<H_<<std::endl;

	lattice("./","lattice");
	//plot_band_structure();
}

void LadderFree::lattice(std::string const& path, std::string const& filename){
	compute_H();
	Matrix<int> nb;
	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(path,filename);
	ps.begin(-9,-2,16,2,filename);

	unsigned int n_plot((n_<32?n_:spuc_));
	for(unsigned int i(0);i<n_plot;i++) {
		xy0(0) = i/2;
		xy0(1) = i%2;
		ps.put(xy0(0),xy0(1)+(i%2?0.2:-0.2),"\\tiny{"+my::tostring(i)+"}");
		nb = get_neighbourg(i);

		if(H_(i,nb(0,0))){/*x-link*/
			xy1(0) = nb(0,0)/2;
			xy1(1) = nb(0,0)%2;

			if( H_(i,nb(0,0)) < 0){ color = "red"; }
			else { color = "black"; }

			if(xy1(0)<xy0(0)){
				xy1(0) = xy0(0)+1;
				linestyle="dashed";
			} else { linestyle="solid"; }

			linewidth = my::tostring(std::abs(H_(i,nb(0,0))))+"mm";

			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
			ps.put((xy0(0)+xy1(0))/2.0,xy0(1)-(i%2?0.2:-0.2),"\\tiny{"+my::tostring(H_(i,nb(0,0)))+"}");
		}
		if(i%2 && H_(i,nb(1,0))){/*y-link*/
			xy1(0) = nb(1,0)/2;
			xy1(1) = nb(1,0)%2;

			if( H_(i,nb(1,0)) < 0){ color = "red"; }
			else { color = "black"; }

			linestyle="solid";

			linewidth = my::tostring(std::abs(H_(i,nb(1,0))))+"mm";

			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
			ps.put(xy0(0)+0.2,(xy0(1)+xy1(1))/2.0,"\\tiny{"+my::tostring(H_(i,nb(1,0)))+"}");
		}
	}

	if(n_plot>spuc_){
		if(n_plot==n_){
			Matrix<double> polygon(4,2);
			polygon(0,0)=-0.3;
			polygon(0,1)=-0.3;
			polygon(1,0)=n_/2-0.1;
			polygon(1,1)=-0.3;
			polygon(2,0)=n_/2-0.1;
			polygon(2,1)=1.3;
			polygon(3,0)=-0.3;
			polygon(3,1)=1.3;
			ps.polygon(polygon,"linecolor=green");
		}

		Matrix<double> polygon(4,2);
		polygon(0,0)=-0.3;
		polygon(0,1)=-0.3;
		polygon(1,0)=spuc_/2-0.1;
		polygon(1,1)=-0.3;
		polygon(2,0)=spuc_/2-0.1;
		polygon(2,1)=1.3;
		polygon(3,0)=-0.3;
		polygon(3,1)=1.3;
		ps.polygon(polygon,"linecolor=blue");
	}

	ps.end(true,true,true);
}
/*}*/

/*{method needed for analysing*/
std::string LadderFree::extract_level_6(){
	(*data_write_)<<N_<<" "<<m_<<" "<<n_<<" "<<bc_<<" "<<asin(J_(1))<<" "<<E_<<IOFiles::endl;

	return filename_;
}
/*}*/
