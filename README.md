<hr>
<h5>DSQMKryST<\h5>

Massively parallel implementation and approaches to simulate quantum dynamics using Krylov subspace techniques

Current version    : v0.1 (March 22, 2017)

<br><hr>
<h3>Synopsis</h3>

DSQMKryST is an application based on implemented distributed memory parallel algorithms in order to provide a computational framework suitable 
for massively parallel supercomputers to study the dynamics of one-dimensional quantum systems.

The main idea is to approximate the solutions to the time-dependent Schroedinger equation by means of Krylov subspace techniques. 

Currently only calculates estimation on the Loschmidt echo, but adding estimation to other quantum observables 
can be integrated with ease (such as magnetization and imbalance, or even entanglement entropy). The modifications required to
accomplish this can be seen [here](https://github.com/mbrenesn/LGT/tree/master)

The application depends on [PETSc](https://www.mcs.anl.gov/petsc/), [SLEPc](http://slepc.upv.es) and [Boost](http://www.boost.org)

<br><hr>
<h3>Get Started</h3>
