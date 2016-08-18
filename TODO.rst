TODO
====

BUGS
----

URGENT
------

+ understand why RandArray is 10% slower than Rand

IDEA
----

+ replace "if" by assert() for checks that avoid a cashing mistake
+ try to use std::copy to copy the content of Vector and Matrix
+ change the unsigned int into size_t for size in matrix,vector...
+ Directory::sort() : implement a better algorithm
+ Lapack.hpp : find a way to determine optimal blocksize -> lwork
+ RST::text : automatically cut lines at the correct position
+ Swarm : try to replace std::vector<std::shared_ptr<Particle> > by Particle*
  or std::vector<std::unique_ptr<Particle> >
+ Variable has a clone method that doesn't use unique_ptr (e.g. MCSystem).
  therefore memory is allocated in GenericVariable and released in Container.
  this is not really good
