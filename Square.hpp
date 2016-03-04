#ifndef DEF_SQUARE
#define DEF_SQUARE

#include "System2D.hpp"

template<typename Type>
class Square: public System2D<Type>{
	public:
		/*!Constructor that organises the n=p^2+q^2 sites (p,q integer)*/
		Square(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Square()=0;

	protected:
		void init_lattice();

	private:
		unsigned int p_;
		unsigned int q_;

		Matrix<double> set_geometry(unsigned int const& n, unsigned int const& spuc, unsigned int& ref3);
		bool reset_pos_in_lattice(Vector<double>& x) const;
		Vector<double> get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d, int& nn_dir) const;
};

/*{constructor*/
template<typename Type>
Square<Type>::Square(Matrix<double> const& ab, unsigned int const& spuc, std::string const& filename):
	System2D<Type>(set_geometry(( (!this->obs_.size() || !this->obs_[0].nlinks()) ?this->n_:0),spuc,this->ref_(3)),ab,spuc,4,4,filename)
{
	this->filename_ += +"-p"+my::tostring(p_)+"-q"+my::tostring(q_);
}

template<typename Type>
Square<Type>::~Square() = default;
/*}*/

/*{protected methods*/
template<typename Type>
void Square<Type>::init_lattice(){
	if( this->obs_.size() && this->obs_[0].nlinks() ){ this->status_ = 2; }
	else {
		if(this->dir_nn_){
			/*{!the directions are given in the cartesian basis
			 *
			 *       (0,1)
			 *         |
			 * (-1,0)--x--(1,0)
			 *         |
			 *      (0,-1)
			 *}*/
			this->dir_nn_[0](0) = 1.0;
			this->dir_nn_[0](1) = 0.0;

			this->dir_nn_[1](0) = 0.0;
			this->dir_nn_[1](1) = 1.0;

			this->dir_nn_[2](0) =-1.0;
			this->dir_nn_[2](1) = 0.0;

			this->dir_nn_[3](0) = 0.0;
			this->dir_nn_[3](1) =-1.0;

			this->x_[0] = this->dir_nn_[0]*0.5*(p_-q_)+this->dir_nn_[1]*0.5*(p_+q_);
			this->x_[0] = this->x_[0]/sqrt(this->x_[0].norm_squared())*0.001;

			Vector<double> x_loop(this->x_[0]);
			for(unsigned int i(1);i<this->n_;i++){
				this->x_[i] = this->x_[i-1] + this->dir_nn_[0];
				reset_pos_in_lattice(this->x_[i]);
				if(my::are_equal(this->x_[i],x_loop)){
					this->x_[i] += this->dir_nn_[1];
					reset_pos_in_lattice(this->x_[i]);
					x_loop = this->x_[i];
				}
				this->x_[i] = this->x_[i].chop();
			}

			this->equivalent_vertex_[0] = this->dir_nn_[0]*0.5*(p_-q_)+this->dir_nn_[1]*0.5*(p_+q_);
			this->equivalent_vertex_[1] = this->dir_nn_[1]*0.5*(p_-q_)+this->dir_nn_[2]*0.5*(p_+q_);
			this->equivalent_vertex_[2] = this->dir_nn_[3]*0.5*(p_-q_)+this->dir_nn_[0]*0.5*(p_+q_);

			if(this->unit_cell_allowed()){ 
				this->status_ = 2; 

				this->set_nn_links(Vector<unsigned int>(1,2));

				/*!sets the bond energy if it has not been set yet*/
				if(this->obs_[0].nlinks() != this->J_.size()){
					if(this->J_.size() == 1){ this->J_.set(this->obs_[0].nlinks(),this->J_(0)); }
					else { std::cerr<<__PRETTY_FUNCTION__<<" : setting J_ is problematic"<<std::endl; }
				}
			}
		}
	}
}
/*}*/

/*{private methods*/
template<typename Type>
Matrix<double> Square<Type>::set_geometry(unsigned int const& n, unsigned int const& spuc, unsigned int& ref3){
	if(n){
		p_ = sqrt(n);
		Matrix<double> tmp;
		if(!ref3 && my::are_equal(sqrt(n),p_)){
			tmp.set(5,2);
			q_ = 0;
		} else {
			for(unsigned int p(0);p<=sqrt(n);p++){
				for(unsigned int q(0);q<p+1;q++){
					if(p*p+q*q==n){
						tmp.set(5,2);
						p_ = p;
						q_ = q;
						p = n;
						q = n;
					}
				}
			}
			if(!q_){ ref3 = 0; }
		}
		if(tmp.ptr()){
			tmp(0,0) =-0.5*(p_-q_);
			tmp(0,1) =-0.5*(p_+q_);
			tmp(1,0) = 0.5*(p_+q_);
			tmp(1,1) =-0.5*(p_-q_);
			tmp(2,0) = 0.5*(p_-q_);
			tmp(2,1) = 0.5*(p_+q_);
			tmp(3,0) =-0.5*(p_+q_);
			tmp(3,1) = 0.5*(p_-q_);
			tmp(4,0) = tmp(0,0);
			tmp(4,1) = tmp(0,1);
			return tmp;
		}

		std::cerr<<__PRETTY_FUNCTION__<<" : unknown geometry (possible sizes)"<<std::endl;
		std::vector<unsigned int> v;
		unsigned int m;
		for(unsigned int p(2);p<2*sqrt(n);p++){
			for(unsigned int q(0);q<p+1;q++){
				m = p*p+q*q;
				if(!(m%spuc)){ v.push_back(m); }
			}
		}
		std::sort(v.begin(),v.end(),std::less<unsigned int>());
		v.erase(std::unique(v.begin(),v.end()),v.end());
		for(unsigned int i(0);i<v.size();i++){ std::cerr<<"n="<<v[i]<<std::endl; }
		std::cerr<<"n=p*p+q*q"<<std::endl;
	}
	return Matrix<double>();
}

template<typename Type>
bool Square<Type>::reset_pos_in_lattice(Vector<double>& x) const {
	if(this->pos_out_of_lattice(x)){
		double a(this->cluster_vertex_(0,0)*x(0)+this->cluster_vertex_(0,1)*x(1));
		double b(this->cluster_vertex_(1,0)*x(0)+this->cluster_vertex_(1,1)*x(1));
		if(a>0){
			if(b>0){ x += this->dir_nn_[1]*p_+this->dir_nn_[2]*q_; }
			else   { x += this->dir_nn_[0]*p_+this->dir_nn_[1]*q_; }
		} else {
			if(b>0){ x += this->dir_nn_[2]*p_+this->dir_nn_[3]*q_; }
			else   { x += this->dir_nn_[3]*p_+this->dir_nn_[0]*q_; }
		}
		reset_pos_in_lattice(x); //misplaced site (e.g. n=18, site 4)
		return true;
	} else { return false; }
}

template<typename Type>
Vector<double> Square<Type>::get_relative_neighbourg_position(unsigned int const& i, unsigned int const& d, int& nn_dir) const {
	nn_dir = d;
	(void)(i);
	return this->dir_nn_[d];
}
/*}*/
#endif
