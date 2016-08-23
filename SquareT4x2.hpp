#ifndef DEF_SQUARET4X2
#define DEF_SQUARET4X2

#include "Square.hpp"

class SquareT4x2: public Square<double>{
	public:
		SquareT4x2(System const& s, Vector<double> const& t);
		~SquareT4x2() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		Vector<double> const t_; //!< hopping terms

		void compute_H();
		void display_results();
		void param_fit_therm_limit(std::string& f, std::string& param, std::string& range);

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;

		std::string get_mc_run_command() const;
};
#endif
