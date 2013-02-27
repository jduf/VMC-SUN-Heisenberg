#ifndef DEF_CHRONO
#define DEF_CHRONO

#include<iostream>
#include<ctime>

class Chrono{
	public:
		/*!sets  tick and t to 0 */
		Chrono();

		/*!returns the time elapsed bewteen the call of tic() and tac() in second*/
		inline double const& get_time() const { return t;};
		/*!starts the clock*/
		void tic(); 
		/*!stops the clock*/
		void tac();

	private:
		clock_t tics; //!< number of ticks bewteen the call of tic() and tac()
		double t; //!< time elapsed in second bewteen the call of tic() and tac()
};

/*constructor*/
/*{*/
Chrono::Chrono()
	:tics(0),t(0.0)
{ }
/*}*/

/*methods*/
/*{*/
void Chrono::tic(){
	tics = std::clock();
}

void Chrono::tac(){
	tics = (std::clock()-t);
	t = tics/double(CLOCKS_PER_SEC);
}
/*}*/

/* operators that takes an instance of the Chrono class and returns its t argument*/
std::ostream& operator<<(std::ostream& flux, Chrono const& t){
	flux << t.get_time();
	return flux;
}
#endif
