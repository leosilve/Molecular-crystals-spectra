TOPDIR=..


EXTLIBS_INCLUDE_DIR=boost_1_57_0 clapack/INCLUDE
INCLUDE_DIR=$(TOPDIR)/include

CC=gcc
CXXFLAGS=-I $(INCLUDE_DIR) -I src $(addprefix -I$(TOPDIR)/ExtLibs/, $(EXTLIBS_INCLUDE_DIR)) -Wall
LDFLAGS = -framework Accelerate -stdlib=libc++ -lc++ 

all: models modeld bands bandd spectras spectrad

models: src/PhononCloud.o src/lattice.o $(INCLUDE_DIR)/$(INCLUDE_DIR)/N_Functions.o src/sym_op.o src/pos.o $(INCLUDE_DIR)/physics.o $(INCLUDE_DIR)/OPutils.o \
	src/ME.o src/vector_per.o src/QuantumState.o src/mpstate.o src/SymmetricDoubleWell.o src/OP_model.o \
	$(INCLUDE_DIR)/Transfer_matrix.o src/models.o 
	$(CC) $(LDFLAGS) $^ -lm -o models.exe

modeld: src/PhononCloud.o src/lattice.o $(INCLUDE_DIR)/N_Functions.o src/sym_op.o src/pos.o $(INCLUDE_DIR)/physics.o $(INCLUDE_DIR)/OPutils.o \
	src/ME.o src/vector_per.o src/QuantumState.o src/mpstate.o src/SymmetricDoubleWell.o src/OP_model.o \
	$(INCLUDE_DIR)/Transfer_matrix.o src/modeld.o 
	$(CC) $(LDFLAGS) $^ -lm -o modeld.exe

bands: src/lattice.o $(INCLUDE_DIR)/N_Functions.o src/sym_op.o src/pos.o $(INCLUDE_DIR)/physics.o $(INCLUDE_DIR)/OPutils.o \
	src/bands.o
	$(CC) $(LDFLAGS) $^ -lm -o bands.exe

bandd: src/lattice.o $(INCLUDE_DIR)/N_Functions.o src/sym_op.o src/pos.o $(INCLUDE_DIR)/physics.o $(INCLUDE_DIR)/OPutils.o \
	src/bandd.o
	$(CC) $(LDFLAGS) $^ -lm -o bandd.exe

spectras: $(INCLUDE_DIR)/N_Functions.o $(INCLUDE_DIR)/physics.o $(INCLUDE_DIR)/OPutils.o $(INCLUDE_DIR)/Transfer_matrix.o src/spectras.o
	$(CC) $(LDFLAGS) $^ -lm -o spectras.exe

spectrad: $(INCLUDE_DIR)/N_Functions.o $(INCLUDE_DIR)/physics.o $(INCLUDE_DIR)/OPutils.o $(INCLUDE_DIR)/Transfer_matrix.o src/spectrad.o
	$(CC) $(LDFLAGS) $^ -lm -o spectrad.exe

PC: src/PhononCloud.o src/lattice.o $(INCLUDE_DIR)/N_Functions.o src/sym_op.o src/pos.o $(INCLUDE_DIR)/physics.o $(INCLUDE_DIR)/OPutils.o \
	src/ME.o src/vector_per.o src/QuantumState.o src/mpstate.o src/SymmetricDoubleWell.o src/PC_main.o
	$(CC) $(LDFLAGS) $^ -lm -o PC.exe 

test: src/PhononCloud.o src/lattice.o $(INCLUDE_DIR)/N_Functions.o src/sym_op.o src/pos.o $(INCLUDE_DIR)/physics.o $(INCLUDE_DIR)/OPutils.o \
	src/ME.o src/vector_per.o src/QuantumState.o src/mpstate.o src/SymmetricDoubleWell.o src/OP_model.o \
	$(INCLUDE_DIR)/Transfer_matrix.o src/test_main.o 
	$(CC) $(LDFLAGS) $^ -lm -o test.exe


	
%.o: %.cpp
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f src/*.o *.exe

.PHONY: utils QS

utils: $(INCLUDE_DIR)/OPutils.cpp
	$(CC) $(CXXFLAGS) -c $< -o $(INCLUDE_DIR)/OPutils.o
	
QS: src/QuantumState.cpp
	$(CC) $(CXXFLAGS) -c $< -o src/QuantumState.o

trans: $(INCLUDE_DIR)/Transfer_matrix.cpp
	$(CC) $(CXXFLAGS) -c $< -o $(INCLUDE_DIR)/Transfer_matrix.o
