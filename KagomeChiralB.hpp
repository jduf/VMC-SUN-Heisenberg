#ifndef DEF_KAGOMECHiralBB
#define DEF_KAGOMECHiralBB

#include "Kagome.hpp"

class KagomeChiralB: public Kagome<std::complex<double> >{
	public:
		KagomeChiralB(System const& s, double const& phi);
		~KagomeChiralB() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	protected:
		double const phi_; //!< flux per triangular plaquette

		void compute_H();
		void display_results();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;
};
#endif
