#ifndef DEF_CHAINPOLYMERIZED
#define DEF_CHAINPOLYMERIZED

#include "Chain.hpp"

class ChainPolymerized: public Chain<double> {
	public:
		ChainPolymerized(unsigned int const& N, unsigned int const& n, unsigned int const& m, int const& bc, Vector<unsigned int> const& ref);
		~ChainPolymerized();

		void create(double const& x, unsigned int const& type);
		void check();
		void save(IOFiles& w) const;
		std::string get_filename() const { return filename_ + "-delta" + tostring(delta_);}
		
	private:
		void compute_P(Matrix<double>& P);
		void compute_T();

		double delta_;
};
#endif
