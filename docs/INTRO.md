<hr>
<h2>DSQMKryST</h2>

Massively parallel implementation and approaches to simulate quantum dynamics using Krylov subspace techniques

Current version    : v0.2 (January, 2018)

<br><hr>
<h3>Synopsis</h3>

DSQMKryST is an application based on implemented distributed memory parallel algorithms in order to provide a computational framework suitable for massively parallel supercomputers to study the dynamics of one-dimensional quantum systems.

The main idea is to approximate the solutions to the time-dependent Schroedinger equation by means of Krylov subspace techniques using a carefully constructed and distributed basis for the Hilbert subspace and a sparse matrix representation of the Hamiltonian of the system. 

Currently calculates estimations on the Loschmidt echo to show the basic scope of the application, but adding estimation to other quantum observables can be integrated with ease (such as magnetization and imbalance, or even entanglement entropy).

The application depends on [PETSc](https://www.mcs.anl.gov/petsc/), [SLEPc](http://slepc.upv.es) and [Boost](http://www.boost.org)

<br><hr>
<h3>Approaches</h3>

This in-depth documentation refers to the DSQMKryST project and it's composed of two different modules: NodeComm and RingComm. Further information about the difference between each of the approaches documented here can be found on the manuscript located in /docs.
