CC     = gcc
CPP    = g++

CFLAGS = -O2 -Wall -fpermissive
SPEED  = -ffast-math
# SPEED  = -O3 -finline-functions  -fexpensive-optimizations -mcpu=i686 -march=i686 -mpentiumpro -fomit-frame-pointer 
#CFLAGS  = -O3 -finline-functions  -fexpensive-optimizations -mcpu=i686 -march=i686 -mpentiumpro -fomit-frame-pointer 
LDFLAGS= -lm




VERS = 0.58a
SRC  = *.c *.cpp *.h Makefile README COPYING

CPP_OBJ_FILES = bridge.o matrix.o estimate.o maximize.o simulate.o errorchk.o allfreq.o model.o linkage.o
OBJ_FILES =  cestring.o distr.o fileio.o nrutil.o pedipet.o ranlib.o lists.o qnorm.c

#
# This is the makefile and it should work
# The program consists of a mix of C and C++ code, 
# that should be linked together to a single program: pedipet
#
# Claus Ekstrøm 2000

.SUFFIXES: .o .cpp .c

# Creates the final program
pedipet: $(OBJ_FILES) $(CPP_OBJ_FILES) 
	$(CPP) $(LDFLAGS) -o pedipet $(OBJ_FILES) $(CPP_OBJ_FILES)

# Orders to compile C++ code
.cpp.o:	$< $*.h
	$(CPP) $(CFLAGS) $(SPEED) -c -o $@ $<

# Orders to compile C code
# Use special flag for some header files so we can link C++ and C code
.c.o:	$< $*.h
	$(CC) $(CFLAGS) $(SPEED) -c -D GCC -o $@ $<

clean:
	rm -f core *.o pedipet *~
	rm -f marker.info
	rm -f mibdred.*
	rm -f pedindex.*
	rm -f freq.info
	rm -f ibd_dist.out
	rm -f haplo.dump


###
###  Here comes Maintainer stuff -- ignore!
###

# Orders to make a backup to floppy. Remember to mount first
diskcopy: 
	cp *.c /mnt/floppy
	cp *.cpp /mnt/floppy
	cp *.h /mnt/floppy
	cp Makefile /mnt/floppy
	cp example.dat /mnt/floppy
	cp example.par /mnt/floppy
	cp example.run /mnt/floppy
	gzip pedipet
	cp pedipet.gz /mnt/floppy
	gunzip pedipet.gz

distribution: pedipet
	@ls $(SRC) | sed s:^:pedipet-$(VERS)/: >MANIFEST
	@(cd ..; ln -s pedipet pedipet-$(VERS))
	(cd ..; tar -czvf pedipet/pedipet-$(VERS).tgz `cat pedipet/MANIFEST`)
	@(cd ..; rm pedipet-$(VERS))
