#ifndef DEF_CHAIN
#define DEF_CHAIN

#include "CreateSystem.hpp"

class Chain: public CreateSystem<double>{
	public:
		Chain(Parseur& P);
		~Chain();

	private:
		void compute_H();
		void compute_P();
		void compute_T();
		void compute_EVec();
		void compute_band_structure();
		void save(std::string filename);

		Matrix<double> Px_;
};
#endif
