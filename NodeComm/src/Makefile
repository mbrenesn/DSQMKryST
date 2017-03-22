default : NodeComm.x

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

CLINKER=mpicxx
CXX=mpicxx
CXXFLAGS=-g -O3 -mavx -DNDEBUG
LD=mpicxx
LDFLAGS=-g -O3 -mavx -DNDEBUG

NodeComm.x : Environment.o Basis.o SparseOp.o TimeEvo.o main.o
	-${CLINKER} $(LDFLAGS) $^ -o aubry_NC.x ${SLEPC_SYS_LIB}

TimeEvo.o : TimeEvo.cc
	$(CXX) $(CXXFLAGS) -c $< -fPIC  -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fvisibility=hidden -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(BOOST_DIR)

SparseOp.o : SparseOp.cc
	$(CXX) $(CXXFLAGS) -c $< -fPIC  -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fvisibility=hidden -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(BOOST_DIR)

Environment.o : Environment.cc
	$(CXX) $(CXXFLAGS) -c $< -fPIC  -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fvisibility=hidden -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(BOOST_DIR)

Basis.o : Basis.cc
	$(CXX) $(CXXFLAGS) -c $< -fPIC  -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fvisibility=hidden -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(BOOST_DIR)

main.o : main.cc
	$(CXX) $(CXXFLAGS) -c $< -fPIC  -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fvisibility=hidden -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(BOOST_DIR)

%.o : %.cc
	$(CXX) $(CXXFLAGS) -c $< -I$(BOOST_DIR)

wipe : 
	rm -r *.o *.x
