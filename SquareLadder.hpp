#ifndef DEF_SQUARELADDER
#define DEF_SQUARELADDER

#include "Square.hpp"

class SquareLadder: public Square<double>{
	public:
		SquareLadder(System const& s, Vector<double> const& t);
		~SquareLadder() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		Vector<double> const t_; //!< hopping terms

		void compute_H();
		void display_results();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;

		std::string get_mc_run_command() const;
};
#endif
