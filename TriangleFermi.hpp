#ifndef DEF_TRIANGLEFERMI
#define DEF_TRIANGLEFERMI

#include "Triangle.hpp"

class TriangleFermi: public Triangle<double>{
	public:
		TriangleFermi(Parseur& P);
		~TriangleFermi();

	protected:
		void compute_T();

		void save(std::string filename);
};
#endif


