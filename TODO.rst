TODO
====

BUGS
----

+ can't add Vector to Parseur...

URGENT
------

+ PSO check the number of *new* value that is computed and stop if the ratio is
  too low
+ Need to check how MCParticles updates its best known value. check also
  MCParticle::set_bx_via(...). If everything works fine, could remove the found
  variable in MCParticle::set_bx_via

IDEA
----

+ in MCSystem, the bond energy is computed with the value of J : JSiSj. but
  when the results are displayed, the factor J is removed so that it actually
  displays the correlations and not the energy. so it might be smarter to
  compute SiSj instead of JSiSj
+ << >> should exactly the same for cout or IOFiles but print and write do
  different things
+ rename SystemFermionic<->FermionicSystem
+ could define once and for all create() because it's always doing the same thing
+ could define set_ab for a given lattice in the lattice main class, then when
  for instance Square is constructed, the constructor could just take an index
  to that unit cell. likewise, unit_cell_index could be just defined
  accordingly and use a pointer to select the correct method to call. might bug
  for SquareChiral
+ display_results() should take an argument so it plots the whole lattice when
  needed
+ get rid of my::norm_squared
