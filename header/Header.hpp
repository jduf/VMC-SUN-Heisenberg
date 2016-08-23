#ifndef DEF_HEADER
#define DEF_HEADER

#include "Time.hpp"
#include "RST.hpp"

class Header:public RST{
	public:
		/*!Constructor*/
		Header() = default;
		/*!Destructor*/
		virtual ~Header() = default;
		/*{Forbidden*/
		Header(Header const&) = delete;
		Header(Header&&) = delete;
		Header& operator=(Header) = delete;
		/*}*/

		/*!Initializes the header with a title and the time of creation*/
		void init(std::string const& s);

		/*!Adds some text s to the Header*/
		void add(std::string const& s);
		void add(std::string const& name, double const& val);
		void add(std::string const& name, bool const& val);
		void add(std::string const& name, unsigned int const& val);
		void add(std::string const& name, int const& val);
		void add(std::string const& name, std::string const& val);
		void add(std::string const& name, std::complex<double> const& val);
		template<typename Type>
			void add(std::string const& name, Type const& val);

	private:
		/*!Returns the time when a Header is created*/
		std::string when();
};

std::ostream& operator<<(std::ostream& flux, Header const& h);

template<typename Type>
void Header::add(std::string const& name, Type const& t){
	def(name,t.header_def());
}
#endif
