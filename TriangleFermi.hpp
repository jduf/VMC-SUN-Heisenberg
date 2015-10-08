#ifndef DEF_TRIANGLEFERMI
#define DEF_TRIANGLEFERMI

#include "Triangle.hpp"

class TriangleFermi: public Triangle<double>{
	public:
		TriangleFermi(System const& s);
		~TriangleFermi() = default;

		void create(unsigned int const& which_observables);
		void check();

	protected:
		void compute_H();
		void lattice(std::string const& path, std::string const& filename);

		unsigned int match_pos_in_ab(Vector<double> const& x) const;
		Matrix<double> set_ab();
};
#endif
