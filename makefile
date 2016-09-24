FORTRAN = gfortran
SOURCES = driver.f UMAT.f
DEBUG = -g
OBJECTS = $(SOURCES:.f=.o)
EXECUTABLE = mTest.x

all : $(EXECUTABLE) $(OBJECTS)

$(EXECUTABLE) : $(OBJECTS)
	$(FORTRAN) $(OBJECTS) -o $@ 

$(OBJECTS):
	$(FORTRAN) $(SOURCES) -c   

clean :
	rm -f *.out $(EXECUTABLE) *.o

