#ifndef DEF_PARTICLEONGRID
#define DEF_PARTICLEONGRID

#include "List.hpp"
#include "PSO.hpp"

class Measure{
	public:
		Measure(Vector<double> const& x, double const& fx):fx_(fx),x_(x),N_(0){};

		void print(std::ostream& flux) const
		{ flux<<x_<<" : "<<fx_<<std::endl; }

		static bool func(Measure const& a, Measure const& b);

		double fx_;
		Vector<double> x_;
		unsigned int N_;
};

std::ostream& operator<<(std::ostream& flux, Measure const& m);

class ParticleOnGrid: public Particle {
	public:
		ParticleOnGrid() = default;
		~ParticleOnGrid();

		void move(Vector<double> const& bx_all);

		void add_history(std::shared_ptr<Measure> const& m);
		void update(double fx);

	protected:
		List<Measure> particle_history_;
};
#endif
