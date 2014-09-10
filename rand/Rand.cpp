#include "Rand.hpp"

RandDouble::RandDouble(double const& min_inclusive, double const& max_exclusive):
	mt_(std::random_device()()),
	dist_(min_inclusive,max_exclusive)
{}


RandUnsignedInt::RandUnsignedInt(unsigned int const& min_inclusive, unsigned int const& max_inclusive):
	mt_(std::random_device()()),
	dist_(min_inclusive,max_inclusive)
{}

