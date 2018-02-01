default : aubry_NC.x

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

CC_FILES := $(wildcard src/*/*.cc)
OBJ_FILES := $(addprefix obj/,$(notdir $(CC_FILES:.cc=.o)))

CLINKER=mpiicpc
CXX=mpiicpc
CXXFLAGS=-g -O3 -mavx -DNDEBUG
LD=mpiicpc
LDFLAGS=-g -O3 -mavx -DNDEBUG


aubry_NC.x : $(OBJ_FILES)
	-${CLINKER} $(LDFLAGS) $^ -o $@ ${SLEPC_SYS_LIB}

obj/%.o : src/*/%.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $< -fPIC -wd1572 -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fvisibility=hidden -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(BOOST_DIR)

wipe : 
	rm -r obj/*.o *.x
