#ifndef DEF_MCPARTICLE
#define DEF_MCPARTICLE

#include "PSO.hpp"
#include "MCSim.hpp"

class MCParticle: public Particle{
	public:
		/*Default constructor*/
		MCParticle() = default;
		/*Default destructor*/
		~MCParticle() = default;
		/*Forbid copy*/
		MCParticle(MCParticle const&) = delete;
		MCParticle(MCParticle&&) = delete;

		void move(Vector<double> const& bx_all);
		bool update(std::shared_ptr<MCSim> const& new_elem);

		void print() const;

	private:
		List<MCSim> history_;

		static void fuse(MCSim& list_elem, MCSim& new_elem);
};
#endif
