#ifndef DEF_KAGOMEFERMI
#define DEF_KAGOMEFERMI

#include "Kagome.hpp"

template<typename Type>
class KagomeFermi: public Kagome<Type>{
	public:
		KagomeFermi(System const& s);
		~KagomeFermi() = default;

		void create();
		void check();

	protected:
		void compute_H();
		void display_results();
		void lattice();

		Matrix<double> set_ab() const;
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};

template<typename Type>
KagomeFermi<Type>::KagomeFermi(System const& s):
	System(s),
	Kagome<Type>(set_ab(),3,"kagome-fermi")
{
	if(this->status_==3){ this->init_lattice(); }
	if(this->status_==2){
		this->init_fermionic();

		this->system_info_.text("KagomeFermi :");
		this->system_info_.text(" Each color has the same Hamiltonian.");
	}
}

/*{method needed for running*/
template<typename Type>
void KagomeFermi<Type>::compute_H(){
	this->H_.set(this->n_,this->n_,0);

	double t(-1.0);
	for(unsigned int i(0);i<this->obs_[0].nlinks(); i++){
		this->H_(this->obs_[0](i,0),this->obs_[0](i,1)) = (this->obs_[0](i,4)?this->bc_*t:t);
	}
	this->H_ += this->H_.transpose();
}

template<typename Type>
Matrix<double> KagomeFermi<Type>::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 1.0;
	tmp(1,1) = sqrt(3.0);
	return tmp;
}

template<typename Type>
unsigned int KagomeFermi<Type>::match_pos_in_ab(Vector<double> const& x) const {
	Vector<double> match(2,0);
	if(my::are_equal(x,match)){ return 0; }
	match(0) = 0.5;
	if(my::are_equal(x,match)){ return 1; }
	match(1) = 0.5;
	if(my::are_equal(x,match)){ return 2; }
	return 3;
}
/*}*/

/*{method needed for checking*/
template<typename Type>
void KagomeFermi<Type>::lattice(){
	compute_H();

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(this->info_+this->path_+this->dir_,this->filename_);
	ps.begin(-20,-20,20,20,this->filename_);
	ps.polygon(this->cluster_vertex_,"linecolor=green");
	ps.polygon(this->draw_unit_cell(),"linecolor=black");

	Type t;
	unsigned int s0;
	unsigned int s1;
	for(unsigned int i(0);i<this->obs_[0].nlinks();i++){
		s0 = this->obs_[0](i,0);
		xy0 = this->x_[s0];

		s1 = this->obs_[0](i,1);
		xy1 = this->x_[s1];

		t = this->H_(s0,s1);
		if(std::abs(t)>1e-4){
			if((xy0-xy1).norm_squared()>1.0001){
				linestyle = "dashed";
				xy1 = (xy0+this->dir_nn_[this->obs_[0](i,3)]).chop();
				ps.put(xy1(0)-0.20,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
			} else{ linestyle = "solid"; }

			if(my::real(t)>0){ color = "blue"; }
			else             { color = "red"; }
			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}
		if(i%2){ ps.put(xy0(0)-0.20,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); }
	}
	ps.line("-",this->boundary_vertex_[0](0),this->boundary_vertex_[0](1),this->boundary_vertex_[1](0),this->boundary_vertex_[1](1),"linecolor=yellow");
	ps.line("-",this->boundary_vertex_[3](0),this->boundary_vertex_[3](1),this->boundary_vertex_[0](0),this->boundary_vertex_[0](1),"linecolor=yellow");
	ps.end(true,true,true);
}

template<typename Type>
void KagomeFermi<Type>::check(){
	this->info_ = "";
	this->path_ = "";
	this->dir_  = "./";
	this->filename_ ="kagome-fermi";
	display_results();
}

template<typename Type>
void KagomeFermi<Type>::display_results(){
	lattice();
}
/*}*/
#endif
