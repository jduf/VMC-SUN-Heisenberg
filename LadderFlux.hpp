#ifndef DEF_LADDERFLUX
#define DEF_LADDERFLUX

#include "Ladder.hpp"

/*{Description*/
/*!Creates a ladder with uniform hopping parameter
 *
 *  => Flux Ladder<=
 *
 * */
/*}*/
class LadderFlux: public Ladder<std::complex<double> >{
	public:
		LadderFlux(System const& s, Vector<double> const& t, Vector<double> const& flux);
		~LadderFlux() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		Vector<double> t_;
		Vector<double> flux_;

		void compute_H();
		unsigned int set_spuc(Vector<double> const& t, Vector<double> const& flux);

		void display_results();
		void plot(bool const& create_image);
		void lattice();

		std::string extract_level_6();
};
#endif
