#ifndef DEF_SQUAREVCS
#define DEF_SQUAREVCS

#include "Square.hpp"

class SquareVCS: public Square<double>{
	public:
		SquareVCS(System const& s, Vector<double> const& t);
		~SquareVCS() = default;

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
