#ifndef DEF_LADDERFREEFLUX
#define DEF_LADDERFREEFLUX

#include "Ladder.hpp"

/*{Description*/
/*!Creates a ladder with uniform hopping parameter
 *
 *  => FreeFlux Ladder<=
 *
 * */
/*}*/
class LadderFreeFlux: public Ladder<std::complex<double> >{
	public:
		LadderFreeFlux(System const& s, Vector<double> const& t, Vector<double> const& flux);
		~LadderFreeFlux() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();
		void get_wf_symmetries(std::vector<Matrix<int> >& sym) const;

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
