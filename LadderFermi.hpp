#ifndef DEF_LADDERFERMI
#define DEF_LADDERFERMI

#include "Ladder.hpp"

/*{Description*/
/*!Creates a ladder with uniform hopping parameter
 *
 *  => Fermi Ladder<=
 *
 * */
/*}*/
template<typename Type>
class LadderFermi: public Ladder<Type>{
	public:
		LadderFermi(System const& s);
		~LadderFermi() = default;

		void create();
		void check();

	private:
		void compute_H();
		void display_results();
};

template<typename Type>
LadderFermi<Type>::LadderFermi(System const& s):
	System(s),
	Ladder<Type>(1,"ladder-fermi")
{
	if(this->status_==2){
		this->init_fermionic();

		this->system_info_.text("LadderFermi :");
		this->system_info_.item("Each color has the same Hamiltonian.");
		this->system_info_.item("Uniform real hopping term.");
	}
}

/*{method needed for running*/
template<typename Type>
void LadderFermi<Type>::compute_H(){
	this->H_.set(this->n_,this->n_,0);
	for(unsigned int i(0);i<this->obs_[0].nlinks(); i++){
		this->H_(this->obs_[0](i,0),this->obs_[0](i,1)) = (this->obs_[0](i,4)?this->bc_:1);
	}
	this->H_ += this->H_.transpose();
}
/*}*/

/*{method needed for checking*/
template<typename Type>
void LadderFermi<Type>::check(){
	this->info_ = "";
	this->path_ = "";
	this->dir_  = "./";
	this->filename_ ="ladder-fermi";
	display_results();

	//this->plot_band_structure();
}

template<typename Type>
void LadderFermi<Type>::display_results(){
	compute_H();
	Matrix<int> nb;
	std::string color("black");
	std::string linestyle("solid");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(this->path_,this->filename_);
	ps.begin(-9,-10,16,10,this->filename_);
	for(unsigned int i(0);i<this->n_;i++) {
		xy0(0) = i/2;
		xy0(1) = i%2;
		ps.put(xy0(0)-0.20,xy0(1)+0.15,my::tostring(i));
		nb = this->get_neighbourg(i);

		if(nb(0,1)<0){ color = "red"; }
		else { color = "black"; }
		xy1(0) = nb(0,0)/2;
		xy1(1) = nb(0,0)%2;
		if(xy1(0)<xy0(0)){
			xy1(0) = xy0(0)+1;
			linestyle="dashed";
		} else{ linestyle="solid"; }
		/*x-link*/  ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linecolor="+color+",linestyle="+linestyle);
		if(i%2){
			color = "black";
			linestyle="solid";
			xy1(0) = nb(1,0)/2;
			xy1(1) = nb(1,0)%2;
			/*y-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linecolor="+color+",linestyle="+linestyle);
		}
	}

	Matrix<double> polygon(4,2);
	polygon(0,0)=-0.1;
	polygon(0,1)=-0.1;
	polygon(1,0)=this->n_/2-0.1;
	polygon(1,1)=-0.1;
	polygon(2,0)=this->n_/2-0.1;
	polygon(2,1)=1.1;
	polygon(3,0)=-0.1;
	polygon(3,1)=1.1;
	ps.polygon(polygon,"linecolor=green");

	polygon(0,0)=-0.1;
	polygon(0,1)=-0.1;
	polygon(1,0)=0.9;
	polygon(1,1)=-0.1;
	polygon(2,0)=0.9;
	polygon(2,1)=1.1;
	polygon(3,0)=-0.1;
	polygon(3,1)=1.1;
	ps.polygon(polygon,"linecolor=blue");

	ps.end(true,true,true);
}
/*}*/
#endif
