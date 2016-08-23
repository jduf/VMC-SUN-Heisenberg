#ifndef DEF_LADDERFREEREAL
#define DEF_LADDERFREEREAL

#include "Ladder.hpp"

class LadderFree: public Ladder<double>{
	public:
		LadderFree(System const&s, Vector<double> const& t, Vector<double> const& mu);
		~LadderFree() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		Vector<double> const t_; //!< hopping terms
		Vector<double> const mu_;//!< chemical potentials

		void compute_H();
		unsigned int set_spuc(Vector<double> const& t, Vector<double> const& mu, unsigned int const& spuc);

		void display_results();

		std::string extract_level_6();
};
#endif
