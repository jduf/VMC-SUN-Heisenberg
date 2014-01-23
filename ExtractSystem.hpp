#ifndef DEF_EXTRACTSYSTEM
#define DEF_EXTRACTSYSTEM

#include "Container.hpp"

class ExtractSystem{
	public:
		ExtractSystem(std::string filename);

		void extract();
		void extract(Container& input, Container& param);
		void print();

		bool use_complex();

	private:
		FileParser file; 
		Container data;
		Vector<unsigned int> ref;
};
#endif
