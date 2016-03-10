#ifndef DEF_KAGOME
#define DEF_KAGOME

#include "System2D.hpp"

template<typename Type>
class Kagome: public System2D<Type>{
	public:
		/*!Constructor that organises the n=9L^2 or 9L^2 sites (L integer)*/
		Kagome(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Kagome()=0;

	protected:
		void init_lattice();

	private:
		double L_;

		Matrix<double> set_geometry(unsigned int const& n, unsigned int const& spuc);
		bool reset_pos_in_lattice(Vector<double>& x) const;
		Vector<double> get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d, int& nn_dir) const;
};

/*{constructor*/
template<typename Type>
Kagome<Type>::Kagome(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename):
	System2D<Type>(set_geometry(( (!this->obs_.size() || !this->obs_[0].nlinks()) ?this->n_:0),spuc),ab,spuc,4,6,filename)
{}

template<typename Type>
Kagome<Type>::~Kagome() = default;
/*}*/

/*{protected methods*/
template<typename Type>
void Kagome<Type>::init_lattice(){
	if(!this->ref_(4)){ this->status_ = 2; }
	else {
		if(this->dir_nn_){
			/*{!the directions are given in the cartesian basis
			 *
			 * (-1,sqrt(3))/2  (1,sqrt(3))/2
			 *               \ /
			 *        (-1,0)--x--(1,0)
			 *               / \
			 * (-1,-sqrt(3))/2 (1,-sqrt(3))/2
			 *}*/
			this->dir_nn_[0](0) = 1.0;
			this->dir_nn_[0](1) = 0.0;

			this->dir_nn_[1](0) = 0.5;
			this->dir_nn_[1](1) = sqrt(3.0)/2.0;

			this->dir_nn_[2](0) =-0.5;
			this->dir_nn_[2](1) = sqrt(3.0)/2.0;

			this->dir_nn_[3](0) =-1.0;
			this->dir_nn_[3](1) = 0.0;

			this->dir_nn_[4](0) =-0.5;
			this->dir_nn_[4](1) =-sqrt(3.0)/2.0;

			this->dir_nn_[5](0) = 0.5;
			this->dir_nn_[5](1) =-sqrt(3.0)/2.0;

			this->x_[0].set(2);
			this->x_[0](0) =-0.5001; 
			if(this->ref_(3)){ this->x_[0](1) =-0.5001; }
			else             { this->x_[0](1) =-0.0001; }

			Vector<double> x_loop(this->x_[0]);
			for(unsigned int i(1);i<this->n_;i++){
				switch(i%3){
					case 0: { this->x_[i] = this->x_[i-1] + this->dir_nn_[0]; }break;
					case 1: { this->x_[i] = this->x_[i-1] + this->dir_nn_[1]; }break;
					case 2: { this->x_[i] = this->x_[i-1] + this->dir_nn_[5]; }break;
				}
				reset_pos_in_lattice(this->x_[i]);
				if(my::are_equal(this->x_[i],x_loop)){
					this->x_[i] += this->dir_nn_[1]*2.0;
					reset_pos_in_lattice(this->x_[i]);
					x_loop = this->x_[i];
				}
				this->x_[i] = this->x_[i].chop();
			}

			if(this->ref_(3)){
				this->equivalent_vertex_[0] = this->dir_nn_[4]*L_*2.0 + this->dir_nn_[3]*0.2;
				this->equivalent_vertex_[1] = this->dir_nn_[0]*L_*2.0 + this->dir_nn_[3]*0.2;
				this->equivalent_vertex_[2] = this->dir_nn_[2]*L_*2.0 + this->dir_nn_[3]*0.2;
			} else {
				this->equivalent_vertex_[0] = (this->dir_nn_[3]+this->dir_nn_[4])*L_/1.5 + (this->dir_nn_[4]+this->dir_nn_[3])*0.2;
				this->equivalent_vertex_[1] = (this->dir_nn_[5]+this->dir_nn_[0])*L_/1.5 + (this->dir_nn_[4]+this->dir_nn_[3])*0.2;
				this->equivalent_vertex_[2] = (this->dir_nn_[1]+this->dir_nn_[2])*L_/1.5 + (this->dir_nn_[4]+this->dir_nn_[3])*0.2;
			}

			if(this->unit_cell_allowed()){
				this->status_ = 2; 

				if(this->ref_(4)==2){ this->create_energy_obs(Vector<unsigned int>(3,2)); }
				else { this->ref_(4) = 0; }

				/*!sets the bond energy if it has not been set yet*/
				if(this->obs_[0].nlinks() != this->J_.size()){
					if(this->J_.size() == 1){ this->J_.set(this->obs_[0].nlinks(),this->J_(0)); }
					else { std::cerr<<__PRETTY_FUNCTION__<<" : setting J_ is problematic"<<std::endl; }
				}
			}
		} else { std::cerr<<__PRETTY_FUNCTION__<<" required memory has not been allocated"<<std::endl; }
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
			tmp(0,0) =-L_;
			tmp(0,1) =-L_*a*2.0;
			tmp(1,0) =-tmp(0,0);
			tmp(1,1) = tmp(0,1);
			tmp(2,0) = L_*2.0;
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
			tmp(1,1) = tmp(0,1)/2.0;
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

		std::cerr<<__PRETTY_FUNCTION__<<" : unknown geometry (possible sizes)"<<std::endl;
		std::vector<unsigned int> v;
		unsigned int m;
		for(unsigned int i(2);i<20;i++){
			m = 3*i*i;
			if(!(m%spuc)){ v.push_back(m); }
			m = (3*i)*(3*i);
			if(!(m%spuc)){ v.push_back(m); }
		}
		std::sort(v.begin(),v.end(),std::less<unsigned int>());
		v.erase(std::unique(v.begin(),v.end()),v.end());
		for(unsigned int i(0);i<v.size();i++){ std::cerr<<"n="<<v[i]<<std::endl; }
		std::cerr<<"n=9*l*l or 3*l*l"<<std::endl;
	}
	return Matrix<double>();
}

template<typename Type>
bool Kagome<Type>::reset_pos_in_lattice(Vector<double>& x) const {
	if(this->pos_out_of_lattice(x)){
		if(this->ref_(3)){
			double t(tan(M_PI/3.0)*x(0)/x(1));
			if(x(1)>0){
				if(std::abs(t)<1){ x+= (this->dir_nn_[4]+this->dir_nn_[5])*L_*2.0; }
				else {
					if(t>0){       x+= (this->dir_nn_[4]+this->dir_nn_[3])*L_*2.0; }
					else   {       x+= (this->dir_nn_[5]+this->dir_nn_[0])*L_*2.0; }
				}
			} else {
				if(std::abs(t)<1){ x+= (this->dir_nn_[1]+this->dir_nn_[2])*L_*2.0; }
				else {
					if(t>0){       x+= (this->dir_nn_[1]+this->dir_nn_[0])*L_*2.0; }
					else   {       x+= (this->dir_nn_[2]+this->dir_nn_[3])*L_*2.0; }
				}
			}
		} else {
			double t(tan(M_PI/6.0)*x(0)/x(1));
			if(x(0)>0){
				if(std::abs(t)>1){ x+= this->dir_nn_[3]*L_*2.0; }
				else {
					if(t>0){       x+= this->dir_nn_[4]*L_*2.0; }
					else   {       x+= this->dir_nn_[2]*L_*2.0; }
				}
			} else {
				if(std::abs(t)>1){ x+= this->dir_nn_[0]*L_*2.0; }
				else {
					if(t>0){       x+= this->dir_nn_[1]*L_*2.0; }
					else   {       x+= this->dir_nn_[5]*L_*2.0; }
				}
			}
		}
		return true;
	} else { return false; }
}

template<typename Type>
Vector<double> Kagome<Type>::get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d, int& nn_dir) const {
	switch(i%3){
		case 0:
			{
				switch(d){
					case 0: { nn_dir = 0; }break;
					case 1: { nn_dir = 1; }break;
					case 2: { nn_dir = 3; }break;
					case 3: { nn_dir = 4; }break;
				}
			}break;
		case 1:
			{
				switch(d){
					case 0: { nn_dir = 1; }break;
					case 1: { nn_dir = 2; }break;
					case 2: { nn_dir = 4; }break;
					case 3: { nn_dir = 5; }break;
				}
			}break;
		case 2:
			{
				switch(d){
					case 0: { nn_dir = 0; }break;
					case 1: { nn_dir = 2; }break;
					case 2: { nn_dir = 3; }break;
					case 3: { nn_dir = 5; }break;
				}
			}break;
	}
	return this->dir_nn_[nn_dir];
}
/*}*/
#endif
