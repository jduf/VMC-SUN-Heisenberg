#ifndef DEF_TRIANGLECHIRALSG
#define DEF_TRIANGLECHIRALSG

#include "Triangle.hpp"

class TriangleChiralSG: public Triangle<std::complex<double> >{
	public:
		TriangleChiralSG(System const& s, double const& phi);
		~TriangleChiralSG() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		double const phi_; //!< flux per triangular plaquette

		void compute_H();
		void display_results();
		void bond_energy_obs();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const { (void)(x); return 0; }

		std::string get_mc_run_command() const;
};
#endif
