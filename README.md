<hr>
<h2>DSQMKryST</h2>

Massively parallel implementation and approaches to simulate quantum dynamics using Krylov subspace techniques

Current version    : v0.1 (March 22, 2017)

<br><hr>
<h3>Synopsis</h3>

DSQMKryST is an application based on implemented distributed memory parallel algorithms in order to provide a computational framework suitable for massively parallel supercomputers to study the dynamics of one-dimensional quantum systems.

The main idea is to approximate the solutions to the time-dependent Schroedinger equation by means of Krylov subspace techniques. 

Currently calculates estimations on the Loschmidt echo to show the basic scope of the application, but adding estimation to other quantum observables can be integrated with ease (such as magnetization and imbalance, or even entanglement entropy). The modifications required to accomplish this can be seen [here](https://github.com/mbrenesn/LGT/tree/master)

The application depends on [PETSc](https://www.mcs.anl.gov/petsc/), [SLEPc](http://slepc.upv.es) and [Boost](http://www.boost.org)

<br><hr>
<h3>Get Started</h3>

<h5>Instructions for a Linux cluster</h5>

You'll need PETSc, SLEPc and Boost to start with. These may already be installed in your cluster or personal computer, however I encourage you to download and install these libraries on your own as the PETSc build required to compile and execute the application properly requires custom settings - configuration and compilation of these libraries it's easy and gives you the freedom to customize as you want your installation.

You probably won't have root access to the cluster, so go ahead and download the source files in home directory or another directory of your choosing:

```bash
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.7.5.tar.gz
wget http://slepc.upv.es/download/download.php?filename=slepc-3.7.3.tar.gz
wget https://sourceforge.net/projects/boost/files/boost/1.63.0/boost_1_63_0.tar.gz
```

Might want to check for later versions or different URL's if the ones posted are broken or no longer valid.
Extract the files:

```bash
tar -xzf [petsc-compressed-file]
tar -xzf [slepc-compressed-file]
tar -xzf [boost-compressed-file]
```

PETSc and SLEPc require 

<h5>PETSc build</h5>

<h5>SLEPc build</h5>

<h5>Boost</h5>
