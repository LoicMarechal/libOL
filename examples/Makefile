
CC     = gcc
FC     = gfortran
CFLAGS = -Ofast
FFLAGS = -Ofast


all: test_libol1_c # test_libol1_f

test_libol1_c: test_libol1_c.o
	$(CC) $(CFLAGS) $< $(HOME)/lib/$(ARCHI)/libol1.a $(HOME)/lib/$(ARCHI)/libmesh6.a -lm -o $@

test_libol1_c.o: test_libol1.c $(HOME)/include/libmesh6.h $(HOME)/include/libol1.h
	$(CC) $(CFLAGS) -c -I $(HOME)/include test_libol1.c -o $@

#test_libol1_f: test_libol1_f.o
#	$(FC) $(FFLAGS) $< $(HOME)/lib/$(ARCHI)/libol1.a $(HOME)/lib/$(ARCHI)/libmesh6_fortran_api.a $(HOME)/lib/$(ARCHI)/libmesh6.a -lm -o $@

#test_libol1_f.o: test_libol1.f $(HOME)/include/libmesh6.h $(HOME)/include/libol1.h
#	$(FC) $(FFLAGS) -c -I $(HOME)/include test_libol1.f -o $@

clean:
	rm *.o test_libol1_c #test_libol1_f
