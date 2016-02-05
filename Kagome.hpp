#ifndef DEF_KAGOME
#define DEF_KAGOME

#include "System2DBis.hpp"

template<typename Type>
class Kagome: public System2DBis<Type>{
	public:
		/*!Constructor that organises the n=3L^2 sites (L integer)*/
		Kagome(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Kagome()=0;

	protected:
		void init_lattice();
		void set_obs(int nobs);

	private:
		double L_;

		Matrix<double> set_geometry(unsigned int const& n, unsigned int const& spuc);
		bool reset_pos_in_lattice(Vector<double>& x) const;
		Vector<double> get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d) const;
};

/*{constructor*/
template<typename Type>
Kagome<Type>::Kagome(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename):
	System2DBis<Type>(set_geometry(( (!this->obs_.size() || !this->obs_[0].nlinks()) ?this->n_:0),spuc),ab,spuc,4,filename)
{}

template<typename Type>
Kagome<Type>::~Kagome() = default;
/*}*/

/*{protected methods*/
template<typename Type>
void Kagome<Type>::init_lattice(){
	if(!this->obs_.size() || !this->obs_[0].nlinks()){
		/*{!the directions are given in the cartesian basis for the
		 * sublattice with even site number
		 * 
		 * f--g
		 *  \/  
		 *  d   
		 *  /\
		 * a--b
		 *
		 * a-b = (1,0)
		 * b-d = (-1,sqrt(3))/2
		 * a-d = (1,sqrt(3))/2
		 *}*/
		this->dir_nn_[0](0) = 1.0;
		this->dir_nn_[0](1) = 0.0;

		this->dir_nn_[1](0) =-0.5;
		this->dir_nn_[1](1) = sqrt(3.0)/2.0;

		this->dir_nn_[2](0) = 0.5;
		this->dir_nn_[2](1) = sqrt(3.0)/2.0;

		this->x_[0].set(2);
		this->x_[0](0) = -0.499;
		this->x_[0](1) = -0.499;

		PSTricks ps("./","test");
		ps.begin(-20,-20,20,20,"bla");
		ps.polygon(this->cluster_vertex_,"linecolor=green");
		ps.put(this->x_[0](0),this->x_[0](1),"\\tiny{"+my::tostring(0)+"}"); 

		Vector<double> x_loop(this->x_[0]);
		for(unsigned int i(1);i<this->n_;i++){
			switch(i%3){
				case 0:
					{ this->x_[i] = this->x_[i-1] + this->dir_nn_[0]; }break;
				case 1:
					{ this->x_[i] = this->x_[i-1] + this->dir_nn_[2]; }break;
				case 2:
					{ this->x_[i] = this->x_[i-1] - this->dir_nn_[1]; }break;
			}
			reset_pos_in_lattice(this->x_[i]);
			if(my::are_equal(this->x_[i],x_loop)){
				this->x_[i] += this->dir_nn_[2]*2.0;
				reset_pos_in_lattice(this->x_[i]);
				x_loop = this->x_[i];
			}
			this->x_[i] = this->x_[i].chop();
			ps.put(this->x_[i](0),this->x_[i](1),"\\tiny{"+my::tostring(i)+"}"); 
		}
		ps.end(true,true,true);

		this->boundary_vertex_[0] = (this->dir_nn_[2]-this->dir_nn_[0])*0.25 - (this->dir_nn_[0]+this->dir_nn_[2])*L_/1.5;
		this->boundary_vertex_[1] = this->boundary_vertex_[0] + this->dir_nn_[0]*L_*2.0;
		this->boundary_vertex_[2] = this->boundary_vertex_[1] + this->dir_nn_[2]*L_*2.0;
		this->boundary_vertex_[3] = this->boundary_vertex_[0] + this->dir_nn_[2]*L_*2.0;

		this->set_nn_links(Vector<unsigned int>(3,2));

		/*!sets the bond energy if it has not been set yet*/
		if(this->obs_[0].nlinks() != this->J_.size()){
			if(this->J_.size() == 1){ this->J_.set(this->obs_[0].nlinks(),this->J_(0)); }
			else { std::cerr<<__PRETTY_FUNCTION__<<" : setting J_ is problematic"<<std::endl; }
		}
	}
}

template<typename Type>
void Kagome<Type>::set_obs(int nobs){
	if(nobs>1){ /*the long range correlation*/
		/*missing bond energy*/
		this->obs_.push_back(Observable("Long range correlations",2,this->n_,this->n_));
		for(unsigned int i(0);i<this->n_;i++){
			this->obs_[2](i,0) = 0;
			this->obs_[2](i,1) = i;
			this->obs_[2](i,2) = i;
		}
	}
}
/*}*/

/*{private methods*/
template<typename Type>
Matrix<double> Kagome<Type>::set_geometry(unsigned int const& n, unsigned int const& spuc){
	if(n){
		L_ = sqrt(n/9.0);
		if(my::are_equal(L_,floor(L_))){
			double a(sqrt(3.0)/2);
			Matrix<double> tmp(7,2);
			tmp(0,0) =-0.5*L_;
			tmp(0,1) =-a*L_;
			tmp(1,0) =-tmp(0,0);
			tmp(1,1) = tmp(0,1);
			tmp(2,0) = L_;
			tmp(2,1) = 0.0;
			tmp(3,0) =-tmp(0,0);
			tmp(3,1) =-tmp(0,1);
			tmp(4,0) =-tmp(1,0);
			tmp(4,1) =-tmp(1,1);
			tmp(5,0) =-tmp(2,0);
			tmp(5,1) =-tmp(2,1);
			tmp(6,0) = tmp(0,0);
			tmp(6,1) = tmp(0,1);
			return tmp;
		}
		L_ = sqrt(n/3.0);
		if(my::are_equal(L_,floor(L_))){
			double a(sqrt(3.0)/2.0);
			Matrix<double> tmp(7,2);
			tmp(0,0) = 0.0;
			tmp(0,1) =-L_/a;
			tmp(1,0) = L_;
			tmp(1,1) =-L_*0.5;
			tmp(2,0) = tmp(1,0);
			tmp(2,1) =-tmp(1,1);
			tmp(3,0) =-tmp(0,0);
			tmp(3,1) =-tmp(0,1);
			tmp(4,0) =-tmp(1,0);
			tmp(4,1) =-tmp(1,1);
			tmp(5,0) =-tmp(2,0);
			tmp(5,1) =-tmp(2,1);
			tmp(6,0) = tmp(0,0);
			tmp(6,1) = tmp(0,1);
			return tmp;
		}
	}
	(void)(spuc);	
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown geometry (possible sizes)"<<std::endl; 
	return Matrix<double>();
}

template<typename Type>
bool Kagome<Type>::reset_pos_in_lattice(Vector<double>& x) const {
	if(this->pos_out_of_lattice(x)){
		if(!this->ref_(3)){
			double t(tan(M_PI/6.0)*x(0)/x(1));
			if(x(0)>0){
				if(std::abs(t)>1){ x-= this->dir_nn_[0]*L_*2.0; }
				else {
					if(t>0){       x-= this->dir_nn_[2]*L_*2.0; }
					else   {       x+= this->dir_nn_[1]*L_*2.0; }
				}
			} else {
				if(std::abs(t)>1){ x+= this->dir_nn_[0]*L_*2.0; }
				else {
					if(t>0){       x+= this->dir_nn_[2]*L_*2.0; }
					else   {       x-= this->dir_nn_[1]*L_*2.0; }
				}
			}
		} else {
			double t(tan(M_PI/3.0)*x(0)/x(1));
			if(x(1)>0){
				if(std::abs(t)<1){ x-=(this->dir_nn_[1]-this->dir_nn_[2])*L_*0.5; }
				else {
					if(t>0){       x-=(this->dir_nn_[0]-this->dir_nn_[2])*L_*0.5; }
					else   {       x-= this->dir_nn_[0]*L_*0.5; }
				}
			} else {
				if(std::abs(t)<1){ x+=(this->dir_nn_[1]-this->dir_nn_[2])*L_*0.5; }
				else {
					if(t>0){       x+=(this->dir_nn_[0]-this->dir_nn_[2])*L_*0.5; }
					else   {       x+=(this->dir_nn_[1]-this->dir_nn_[0])*L_*0.5; }
				}
			}
		}
		return true;
	} else { return false; }
}

template<typename Type>
Vector<double> Kagome<Type>::get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d) const {
	switch(i%3){
		case 0:
			{
				switch(d){
					case 0: { return  this->dir_nn_[0]; }break;
					case 1: { return  this->dir_nn_[2]; }break;
					case 2: { return -this->dir_nn_[0]; }break;
					case 3: { return -this->dir_nn_[2]; }break;
				}
			}break;
		case 1:
			{
				switch(d){
					case 0: { return  this->dir_nn_[2]; }break;
					case 1: { return  this->dir_nn_[1]; }break;
					case 2: { return -this->dir_nn_[2]; }break;
					case 3: { return -this->dir_nn_[1]; }break;
				}
			}break;
		case 2:
			{
				switch(d){
					case 0: { return  this->dir_nn_[0]; }break;
					case 1: { return  this->dir_nn_[1]; }break;
					case 2: { return -this->dir_nn_[0]; }break;
					case 3: { return -this->dir_nn_[1]; }break;
				}
			}break;
	}
	return Vector<double>();
}
/*}*/
#endif
