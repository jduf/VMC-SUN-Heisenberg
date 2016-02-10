#ifndef DEF_SQUAREDIMERIZEDBAR
#define DEF_SQUAREDIMERIZEDBAR

#include "Square.hpp"

class SquareDimerizedBar: public Square<double>{
	public:
		SquareDimerizedBar(System const& s, Vector<double> const& t);
		~SquareDimerizedBar() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	protected:
		Vector<double> const t_; //!< hopping terms

		void compute_H();
		void display_results();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
