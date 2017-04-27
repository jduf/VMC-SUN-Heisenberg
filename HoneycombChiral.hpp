#ifndef DEF_HONEYCOMBCHIRAL
#define DEF_HONEYCOMBCHIRAL

#include "Honeycomb.hpp"

class HoneycombChiral: public Honeycomb<std::complex<double> >{
	public:
		HoneycombChiral(System const& s, double const& phi);
		~HoneycombChiral() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		double const phi_; //!< flux per hexagonal plaquette

		void compute_H();
		void display_results();
		void param_fit_therm_limit(std::string& f, std::string& param, std::string& range);
		
		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;

		std::string extract_level_2();
};
#endif
