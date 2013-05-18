#ifndef DEF_CHRONO
#define DEF_CHRONO

#include<iostream>
#include<ctime>

class Chrono{
	public:
		/*!sets  tick and t to 0 */
		Chrono();

		/*!returns the time elapsed bewteen the call of tic() and tac() in second*/
		double const& get_t();
		/*!starts the clock*/
		void tic(); 
		/*!stops the clock*/
		void tac();
		/*!returns true if t>sec*/
		bool time_limit_reached(double sec) const;

	private:
		clock_t tics; //!< number of ticks bewteen the calls of tic() and tac()
		double t; //!< time elapsed in second bewteen the calls of tic() and tac()
};

/*!operators that takes an instance of the Chrono class and returns its t argument*/
std::ostream& operator<<(std::ostream& flux, Chrono& t);

#endif
