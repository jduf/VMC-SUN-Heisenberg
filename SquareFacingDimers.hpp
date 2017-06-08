#ifndef DEF_SQUAREFACINGDIMERS
#define DEF_SQUAREFACINGDIMERS

#include "Square.hpp"

class SquareFacingDimers: public Square<double>{
	public:
		SquareFacingDimers(System const& s, Vector<double> const& t);
		~SquareFacingDimers() = default;

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

		std::string extract_level_9();
};
#endif
