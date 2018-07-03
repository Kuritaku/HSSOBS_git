###############################################################################
#
# Makefile grand control
#
###############################################################################
TOPDIR = .

TARGET = hssdriver

FC     = ifort 
FCFLAGS= -assume byterecl -traceback -convert big_endian
#FCFLAGS= -assume byterecl -traceback -check all -CB -fpe0

ICOGDIR=$(TOPDIR)/icohgrid

OBJECTS = mod_hssconst.o  \
					mod_common_prepbufr.o \
					mod_avefid.o  \
					mod_int2fname.o  \
					mod_ico2latlon.o \
					mod_setup_zlayer.o \
					mod_conv_vcoords_fnl.o \
					mod_qs.o  \
					mod_llselect.o mod_getrgn.o \
					mod_dlatlon.o mod_distdeg.o  mod_gridsearch.o \
					mod_hssll.o \
					mod_hssrgn.o \
					mod_compave.o \
					mod_icoid2ll.o \
					prg_hssobs_driver.o

all: $(TARGET)
			./$(TARGET)

%.o: %.f90
			$(FC) -c $(FCFLAGS) $<

%.mod: %.f90 %.o
			@:

$(TARGET): $(OBJECTS)
			$(FC) -o $@ $(OBJECTS) $(FCFLAGS)


clean:
		rm -rf  $(ICOGDIR)
		rm -f   *.o  *.mod *.cnf  hssdriver
