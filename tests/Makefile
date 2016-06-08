
CC     = /usr/bin/gcc
FC     = gfortran
CFLAGS = -Ofast -flto
FFLAGS = -Ofast


all: test_libol1 test_fortran test_lemma bug_lemma

test_libol1: test_libol1.o
	$(CC) $(CFLAGS) $< $(HOME)/lib/$(ARCHI)/libol1.a $(HOME)/lib/$(ARCHI)/libmesh6.a -lm -o $@

test_libol1.o: test_libol1.c $(HOME)/include/libmesh6.h $(HOME)/include/libol1.h
	$(CC) $(CFLAGS) -c -I $(HOME)/include test_libol1.c

test_fortran: test_fortran.o
	$(FC) $(FFLAGS) $< $(HOME)/lib/$(ARCHI)/libol1.a $(HOME)/lib/$(ARCHI)/libmesh6_fortran_api.a $(HOME)/lib/$(ARCHI)/libmesh6.a -lm -o $@

test_fortran.o: test_fortran.f $(HOME)/include/libmesh6.h $(HOME)/include/libol1.h
	$(FC) $(FFLAGS) -c -I $(HOME)/include test_fortran.f 

test_lemma: test_lemma.o
	$(CC) $(CFLAGS) $< $(HOME)/lib/$(ARCHI)/libol1.a $(HOME)/lib/$(ARCHI)/libmesh6.a -lm -o $@

test_lemma.o: test_lemma.c $(HOME)/include/libmesh6.h $(HOME)/include/libol1.h
	$(CC) $(CFLAGS) -c -I $(HOME)/include test_lemma.c

bug_lemma: bug_lemma.o
	$(CC) $(CFLAGS) $< $(HOME)/lib/$(ARCHI)/libol1.a $(HOME)/lib/$(ARCHI)/libmesh6.a -lm -o $@

bug_lemma.o: bug_lemma.c $(HOME)/include/libmesh6.h $(HOME)/include/libol1.h
	$(CC) $(CFLAGS) -c -I $(HOME)/include bug_lemma.c

clean:
	rm *.o test_libol1 test_fortran test_lemma bug_lemma
