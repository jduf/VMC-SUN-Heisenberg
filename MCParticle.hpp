#ifndef DEF_MCPARTICLE
#define DEF_MCPARTICLE

#include "MCSim.hpp"
#include "PSO.hpp"
#include "List.hpp"

class MCParticle: public Particle{
	public:
		/*!Default constructor*/
		MCParticle() = default;
		/*!Default destructor*/
		~MCParticle() = default;
		/*{Forbid copy*/
		MCParticle(MCParticle const&) = delete;
		MCParticle(MCParticle&&) = delete;
		MCParticle& operator=(MCParticle) = delete;
		/*}*/

		void move(Vector<double> const& bx_all);
		void print() const;

		bool update(std::shared_ptr<MCSim> const& new_elem);
		void add_to_history(std::shared_ptr<MCSim> const& new_elem);
		bool select_new_best();
		void clear_history(){ history_.set(); }
		void set_ps(Vector<double>* ps){ ps_ = ps; }
		Vector<double> get_param() const;

	private:
		List<MCSim> history_;
		unsigned int const update_now_ = 10;
		unsigned int Nupdate_ = 0;
		Vector<double>* ps_ = NULL;

		/*!sets Particle::bx_ to the correct index of ps_*/
		void set_bx_via(Vector<double> const& param);
};
#endif
