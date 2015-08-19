#ifndef DEF_KAGOMEVBC
#define DEF_KAGOMEVBC

#include "Kagome.hpp"

class KagomeVBC: public Kagome<std::complex<double> >{
	public:
		KagomeVBC(System const& s);
		~KagomeVBC() = default;

		void create();
		void check();

	protected:
		void compute_H();
		void lattice(std::string const& path);
		std::string extract_level_7();
		std::string extract_level_6();
		std::string extract_level_4();
};
#endif
