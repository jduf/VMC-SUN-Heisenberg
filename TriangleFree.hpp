#ifndef DEF_TRIANGLEFREE
#define DEF_TRIANGLEFREE

#include "Triangle.hpp"

class TriangleFree: public Triangle<double>{
	public:
		TriangleFree(System const& s, Vector<double> const& t, Vector<double> const& mu);
		~TriangleFree() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	protected:
		Vector<double> const t_; //!< hopping terms
		Vector<double> const mu_;//!< chemical potentials

		void compute_H();
		void display_results();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
