
HYDRO    = euler-helmeos
REACTIONS = pyna13
#paul.h needs REACTIONS directory info

CC = gcc
GFORT = gfortran

OBJ = WD_maker.o helmeos.o $(HYDRO).o

default: WD_maker.exe

%.o: %.c
	$(CC) -c $<

helmeos.o: ../../Hydro/Helmeos/helmeos.f90 
	$(GFORT) -c ../../Hydro/Helmeos/helmeos.f90

$(HYDRO).o : ../../Hydro/$(HYDRO).c ../../paul.h
	$(CC) -c ../../Hydro/$(HYDRO).c

WD_maker.exe: $(OBJ)
	$(GFORT) -o WD_maker.exe $(OBJ)

clean:
	rm -f *.o WD_maker.exe
