#ifndef DEF_MCPARTICLE
#define DEF_MCPARTICLE

#include "MCSim.hpp"
#include "PSO.hpp"

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

		void set_ps(Vector<double> const* const ps){ ps_ = ps; }
		void set_symmetry(Matrix<int> const& sym){ sym_ = sym; }

		/*!Could remove the test within this function if I don't see any bug*/
		Vector<double> get_param() const;

		//might also crash if x_[i] is set to 0
		void move(Vector<double> const& bx_all);
		bool update(std::shared_ptr<MCSim> const& new_elem);
		void clear_history(){ history_.set(); }
		bool select_new_best();

		void print() const;

	private:
		List<MCSim> history_;
		Matrix<int> sym_;
		Vector<double> const* ps_ 	   = NULL;
		unsigned int Nupdate_		   = 0;
		unsigned int const update_now_ = 10;

		/*!sets Particle::bx_ to the correct index of ps_*/
		void set_bx_via(Vector<double> param);
		/*!use update() to add a MCSim because also update fbx_ and bx_ */
		void add_to_history(std::shared_ptr<MCSim> const& new_elem){ history_.add_end(new_elem); }
};
#endif
