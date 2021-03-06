#ifndef DEF_OBSERVABLES
#define DEF_OBSERVABLES

#include "Sampling.hpp"

class Observable{
	public:
		/*!Constructor*/
		Observable(std::string const& name, unsigned int const& type, unsigned int const& nval, Matrix<int> const& links, unsigned int const& modulo=0, unsigned int const& B=50, unsigned int const& b=5, bool const& conv=false);
		/*!Constructor*/
		Observable(std::string const& name, unsigned int const& type, unsigned int const& nval, unsigned int const& nlinks, unsigned int const& modulo=0, unsigned int const& B=50, unsigned int const& b=5, bool const& conv=false);
		/*!Constructor that reads from file*/
		Observable(IOFiles& r);
		/*!Copy constructor*/
		Observable(Observable const&);
		/*!Move constructor*/
		Observable(Observable&& obs);
		/*!Default destructor*/
		virtual ~Observable();
		/*!Assignment operator that uses the copy-and-swap idiom*/
		Observable& operator=(Observable obs);
		/*{Forbidden*/
		Observable() = delete;
		/*}*/

		Matrix<int> const& get_links() const { return links_; }
		std::string const& get_name() const { return name_; }
		unsigned int const& nval() const { return nval_; }
		unsigned int const& nlinks() const { return links_.row(); }
		unsigned int const& get_type() const { return type_; }

		void merge(Observable& obs);
		void delete_binning();
		void complete_analysis();

		void set_x(double const& val);
		void add(unsigned int const& i, double const& val);
		void add_sample();

		Data<double> const& operator[](unsigned int const& i) const { assert(i<nval_); return val_[i]; }
		Data<double>& operator[](unsigned int const& i){ assert(i<nval_); return val_[i]; }

		int const& operator()(unsigned int const& i, unsigned int const& j) const { return links_(i,j); }
		int& operator()(unsigned int const& i, unsigned int const& j){ return links_(i,j); }

		void print(bool const& all) const;
		void write(IOFiles& w) const;

		/*!Provides an easy solution to measure the energy and the bond energy at the same time*/
		void combine_measurement(bool const& combine);
		void reset();
		void remove_links(){ links_.set(); }

		static unsigned int const number_of_observables_defined = 5; //!< Number of different type of observables defined

	protected:
		std::string name_;
		unsigned int type_;
		unsigned int modulo_;
		unsigned int nval_;
		Matrix<int> links_;
		Data<double>* val_;

	private:
		void set(unsigned int const& B, unsigned int const& b, bool const& conv);
		void swap_to_assign(Observable& ds1, Observable& ds2);
};

IOFiles& operator<<(IOFiles& w, Observable const& obs);
IOFiles& operator>>(IOFiles& r, Observable& obs);
std::ostream& operator<<(std::ostream& flux, Observable const& obs);
#endif
