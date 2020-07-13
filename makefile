CCD.exe:ccd.o main.o
	mpicc $(CFLAGS) -o CCD.exe ccd.o main.o
ccd.o:ccd.c ccd.h
	mpicc $(CFLAGS) -c ccd.c
main.o:main.c ccd.h
	mpicc $(CFLAGS) -c main.c

clean:
	rm *.o -f
