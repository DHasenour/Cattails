
INITIAL  = readtable
HYDRO    = euler-helmeos
REACTIONS = pyna13
OUTPUT   = ascii
REPORT	 = report
#paul.h needs REACTIONS directory info

UNAME = $(shell uname)
ifeq ($(UNAME),Linux)
H55 = /home/install/app/hdf5
endif
ifeq ($(UNAME),Darwin)
H55 = /opt/homebrew
endif

CC = mpicc
FLAGS = -O3 -g

GFORT = gfortran
CXX = mpicxx

INC = -I$(H55)/include

OBJ = main.o mpisetup.o profiler.o readpar.o domain.o gridsetup.o geometry.o exchange.o misc.o timestep.o onestep.o riemann.o boundary.o gravity.o nozzle.o plm.o helmeos.o $(INITIAL).o $(OUTPUT).o $(REPORT).o $(HYDRO).o reacstep.o nse_solver.o get_rhs_jacobn.o actual_network_data.o

default: cattails.exe

%.o: %.c paul.h
	$(CC) $(FLAGS) $(INC) -c $<

helmeos.o: Hydro/Helmeos/helmeos.f90 
	$(GFORT) -c Hydro/Helmeos/helmeos.f90

$(TIMESTEP).o: Timestep/$(TIMESTEP).c paul.h
	$(CC) $(FLAGS) $(INC) -c Timestep/$(TIMESTEP).c

$(INITIAL).o : Initial/$(INITIAL).c paul.h
	$(CC) $(FLAGS) $(INC) -c Initial/$(INITIAL).c

$(HYDRO).o : Hydro/$(HYDRO).c paul.h
	$(CC) $(FLAGS) $(INC) -c Hydro/$(HYDRO).c

$(OUTPUT).o : Output/$(OUTPUT).c paul.h
	$(CC) $(FLAGS) $(INC) -c Output/$(OUTPUT).c

$(REPORT).o : Report/$(REPORT).c paul.h
	$(CC) $(FLAGS) $(INC) -c Report/$(REPORT).c

reacstep.o : Reactions/$(REACTIONS)/reacstep.c paul.h
	$(CC) $(FLAGS) -c Reactions/$(REACTIONS)/reacstep.c

nse_solver.o : Reactions/$(REACTIONS)/nse_solver.c paul.h
	$(CC) $(FLAGS) -c Reactions/$(REACTIONS)/nse_solver.c

get_rhs_jacobn.o : Reactions/$(REACTIONS)/get_rhs_jacobn.cpp
	$(CXX) -std=c++20 $(FLAGS) -I Reactions/$(REACTIONS)/headers -D SCREENING -c Reactions/$(REACTIONS)/get_rhs_jacobn.cpp

actual_network_data.o : Reactions/$(REACTIONS)/actual_network_data.cpp
	$(CXX) -std=c++20 $(FLAGS) -I Reactions/$(REACTIONS)/headers -D SCREENING -c Reactions/$(REACTIONS)/actual_network_data.cpp
#-D SCREENING
cattails.exe: $(OBJ) paul.h
	$(CXX) -std=c++20 $(FLAGS) -D SCREENING -o cattails.exe $(OBJ) -lgfortran

clean:
	rm -f *.o cattails.exe
