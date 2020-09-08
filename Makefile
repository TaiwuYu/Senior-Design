CC=g++ 

HEADERFILE=Gchemical.h Gfftw3.h Ggradient.h GInitConfig.h Gmt.h Gelastic.h

OBJ=Gchemical.o Gfftw3.o Ggradient.o GInitConfig.o Gmt.o Pf-MT.o Gelastic.o

LDFLAGS  = -lfftw3  
CFLAGS = ""

PROJECT_NAME=Hphase_model

%.o: %.c  $(HEADERFILE)
	$(CC) -c -o $@ $< $(LDFLAGS) 

all: ${PROJECT_NAME}

${PROJECT_NAME}: $(OBJ) $(HEADERFILE)
	$(CC)  -o $@.exe  $(OBJ) $(LDFLAGS)

clean:
	$(RM) -rf *.o Hphase_model.exe a.out

cleanall:
	$(RM) -rf *.bin *.dat *.o Hphase_model.exe a.out *~ *.bin *.vtk *.dat

cleandata:
	$(RM) -rf output_files/*.bin output_files/*.vtk output_files/*.dat
