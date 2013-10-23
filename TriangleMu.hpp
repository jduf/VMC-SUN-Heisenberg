#ifndef DEF_TRIANGLEMU
#define DEF_TRIANGLEMU

#include "Triangle.hpp"

class TriangleMu: public Triangle<double>{
	public:
		TriangleMu(Parseur& P);
		~TriangleMu();

	protected:
		double mu_;
		
		void compute_T(unsigned int spin);
		void save(std::string filename);
};
#endif

