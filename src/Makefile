#
# National Center for Supercomputing Applications
# University of Illinois at Urbana-Champaign
#
# Large-Scale Agent-Based Social Simulation
# Les Gasser, NCSA Fellow
#
# Author: Santiago Nunez-Corrales
#

IDIR=../include
LDIR=../lib
CC = gcc
CFLAGS=-O2 -I$(IDIR)
LDFLAGS=-shared
ODIR=../obj
RANSRC=../rand
TARGET_LIB=libpse.so

_RNGDEPS = rnglib.h ranlib.h
RNGDEPS = $(patsubst %,$(IDIR)/%,$(_RNGDEPS))

_RNGOBJ = rnglib.o ranlib.o
RNGOBJ = $(patsubst %,$(ODIR)/%,$(_RNGOBJ))

_PSEDEPS = pse.h
PSEDEPS = $(patsubst %,$(IDIR)/%,$(_PSEDEPS))

_PSEOBJ = psedict.o
PSEOBJ = $(patsubst %,$(ODIR)/%,$(_PSEOBJ))

_PSEDICTDEPS = psedict.h
PSEDICTDEPS = $(patsubst %,$(IDIR)/%,$(_PSEDICTDEPS))

_PSEDICTOBJ = psedict.o
PSEDICTOBJ = $(patsubst %,$(ODIR)/%,$(_PSEDICTOBJ))

$(ODIR)/%.o: $(RANSRC)/%.c $(RNGDEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: %.c $(PSEDEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
	
$(ODIR)/%.o: %.c $(PSEDICTDEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
	
$(LDIR)/$(TARGET_LIB): $(RNGOBJ) $(PSEOBJ)
	$(CC) ${LDFLAGS} -o $@ $^

.PHONY: all
all: $(LDIR)/$(TARGET_LIB)	

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ $(LDIR)/*.so
