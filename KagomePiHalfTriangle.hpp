#ifndef DEF_KAGOMEPIHALFTRIANGLE
#define DEF_KAGOMEPIHALFTRIANGLE

#include "Kagome.hpp"

class KagomePiHalfTriangle: public Kagome<std::complex<double> >{
	public:
		KagomePiHalfTriangle(System const& s, Vector<double> const& t, Vector<double> const& phi);
		~KagomePiHalfTriangle() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		Vector<double> const t_; //!< hopping terms
		Vector<double> const phi_;//!< phases

		void compute_H();
		void display_results();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;
};
#endif

