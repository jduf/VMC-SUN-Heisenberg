#ifndef DEF_EXTRACTSYSTEM
#define DEF_EXTRACTSYSTEM

#include "Container.hpp"

/*!Extract the System from a file created with GenericSystem.hpp. The extracted
 * System is stored in two containers, one needed by MonteCarlo.hpp (input) and
 * one needed to save the simulation (param)
*/
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
