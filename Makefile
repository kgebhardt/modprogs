F77=g77
FFLAGS = -c -O3
LFLAGS = -O3
HOSTLIBS= -lm
PGPLOT =  /u/gebhardt/lib/pglib/libpgplot.a -L/usr/X11R6/lib -lX11
LIBDIR  = /a/gebhardt/mod/modprogs/
QUEST   = /u/gebhardt/lib/libquest/libquest.o
NUMREC  = /u/gebhardt/lib/numrec/numrec.a

# RULES:
.SUFFIXES: .o .f
.f.o:
	$(F77) $(FFLAGS) $<

all: gden gdenconv library sos2vor vor2vphase model vlook plot.libsos seecor seecor2 clean

clean : 
	rm *.o

#------------------------------------------------------------------------------#
# makefile for library

library : $(LIBDIR)library.o $(LIBDIR)area.o $(LIBDIR)binners.o $(LIBDIR)binset.o $(LIBDIR)dataread.o $(LIBDIR)density.o $(LIBDIR)densitywrite.o $(LIBDIR)derivs.o $(LIBDIR)ekin.o $(LIBDIR)energy.o $(LIBDIR)force.o $(LIBDIR)forcefix.o $(LIBDIR)galaxyread.o $(LIBDIR)librarian.o $(LIBDIR)moment.o $(LIBDIR)orbit1.o $(LIBDIR)orbit2.o $(LIBDIR)phasewrite.o $(LIBDIR)phasvol.o $(LIBDIR)pl.o $(LIBDIR)potential.o $(LIBDIR)ran1.o $(LIBDIR)rk4.o $(LIBDIR)rtbis.o $(LIBDIR)sort2.o $(LIBDIR)sosmom.o $(LIBDIR)spacel.o $(LIBDIR)step.o $(LIBDIR)tables.o $(LIBDIR)tablesread.o $(LIBDIR)tableswrite.o $(LIBDIR)totmlwrite.o $(LIBDIR)haloread.o $(LIBDIR)halotables.o $(LIBDIR)addtables.o $(LIBDIR)potwrite.o $(LIBDIR)halodens.o
	$(F77) $(LFLAGS) -o library.x library.o area.o binners.o binset.o dataread.o density.o densitywrite.o derivs.o ekin.o energy.o force.o forcefix.o galaxyread.o librarian.o moment.o orbit1.o orbit2.o phasewrite.o phasvol.o pl.o potential.o ran1.o rk4.o rtbis.o sort2.o sosmom.o spacel.o step.o tables.o tablesread.o tableswrite.o totmlwrite.o haloread.o halotables.o addtables.o potwrite.o halodens.o

library.o area.o binners.o binset.o dataread.o density.o densitywrite.o derivs.o ekin.o energy.o force.o forcefix.o galaxyread.o librarian.o moment.o orbit1.o orbit2.o phasewrite.o phasvol.o potential.o spacel.o step.o tables.o tablesread.o tableswrite.o totmlwrite.o haloread.o halotables.o adtables.o potwrite.o halodens.o: libdefs.h


#------------------------------------------------------------------------------#
# makefile for oneorbit

oneorbit : $(LIBDIR)oneorbit.o $(LIBDIR)binners.o $(LIBDIR)binset.o $(LIBDIR)dataread.o $(LIBDIR)derivs.o $(LIBDIR)ekin.o $(LIBDIR)energy.o $(LIBDIR)force.o $(LIBDIR)forcefix.o $(LIBDIR)forcewrite.o $(LIBDIR)haloread.o $(LIBDIR)pl.o $(LIBDIR)galaxyread.o $(LIBDIR)potential.o $(LIBDIR)density.o $(LIBDIR)rk4.o $(LIBDIR)rtbis.o $(LIBDIR)sort2.o $(LIBDIR)sosmom.o $(LIBDIR)spacel.o $(LIBDIR)step.o $(LIBDIR)tables.o $(LIBDIR)tablesread.o $(LIBDIR)tableswrite.o
	$(F77) $(LFLAGS) -o oneorbit.x oneorbit.o binners.o binset.o dataread.o derivs.o ekin.o energy.o force.o forcefix.o forcewrite.o haloread.o pl.o galaxyread.o potential.o density.o rk4.o rtbis.o sort2.o sosmom.o spacel.o step.o tables.o tablesread.o tableswrite.o $(HOSTLIBS) $(PGPLOT) ; rm *.o

oneorbit.o binset.o dataread.o derivs.o ekin.o energy.o galaxyread.o potential.o density.o rk4.o rtbis.o sort2.o spacel.o step.o tables.o tablesread.o tableswrite.o : libdefs.h

#-----------------------------------------------------------------------------#
# makefile for model

model : $(LIBDIR)model.o $(LIBDIR)binners.o $(LIBDIR)binset.o $(LIBDIR)cmatrix.o $(LIBDIR)comparewrite.o $(LIBDIR)comparer.o $(LIBDIR)d3coarse.o $(LIBDIR)d3read.o $(LIBDIR)entropy.o $(LIBDIR)filter.o $(LIBDIR)galreadm.o $(LIBDIR)getfwhm.o $(LIBDIR)getlowright.o $(LIBDIR)irivibin.o $(LIBDIR)iterread.o $(LIBDIR)iterwrite.o $(LIBDIR)linpack.o $(LIBDIR)minmax.o $(LIBDIR)moment.o $(LIBDIR)phaseread.o $(LIBDIR)qual.o $(LIBDIR)smcoarse.o $(LIBDIR)spear.o $(LIBDIR)spline.o $(LIBDIR)summ.o $(LIBDIR)totmlread.o $(LIBDIR)vdataread.o $(LIBDIR)velbin.o $(LIBDIR)velmtod.o $(LIBDIR)vlibread.o $(LIBDIR)weightread.o $(LIBDIR)weightwrite.o $(LIBDIR)xlibread.o
	$(F77) $(LFLAGS) -o model.x model.o binners.o binset.o cmatrix.o comparewrite.o comparer.o d3coarse.o d3read.o entropy.o filter.o galreadm.o getfwhm.o getlowright.o irivibin.o iterread.o iterwrite.o linpack.o minmax.o moment.o phaseread.o qual.o smcoarse.o spear.o spline.o summ.o totmlread.o vdataread.o velbin.o velmtod.o vlibread.o weightread.o weightwrite.o xlibread.o

binset.o cmatrix.o comparewrite.o comparer.o d3coarse.o d3read.o entropy.o filter.o galreadm.o getfwhm.o getlowright.o irivibin.o iterread.o iterwrite.o minmax.o model.o phaseread.o qual.o smcoarse.o spear.o spline.o summ.o totmlread.o vdataread.o velbin.o velmtod.o vlibread.o weightread.o weightwrite.o xlibread.o : moddefs.h

# makefile for vlook

vlook : $(LIBDIR)vlook.o $(LIBDIR)binners.o $(LIBDIR)binset.o $(LIBDIR)chiprob.o $(LIBDIR)cmatrix.o $(LIBDIR)comparewrite.o $(LIBDIR)comparer.o $(LIBDIR)entropy.o $(LIBDIR)filter.o $(LIBDIR)fitherm.o $(LIBDIR)galreadm.o $(LIBDIR)getfwhm.o $(LIBDIR)getlowright.o $(LIBDIR)irivibin.o $(LIBDIR)iterread.o $(LIBDIR)iterwrite.o $(LIBDIR)linpack.o $(LIBDIR)minmax.o $(LIBDIR)moment.o $(LIBDIR)phaseread.o $(LIBDIR)qual.o $(LIBDIR)smcoarse.o $(LIBDIR)spacem.o $(LIBDIR)spear.o $(LIBDIR)spline.o $(LIBDIR)summ.o $(LIBDIR)totmlread.o $(LIBDIR)vaniread.o $(LIBDIR)vdataread.o $(LIBDIR)velbin.o $(LIBDIR)velmtod.o $(LIBDIR)vlibread.o $(LIBDIR)weightread.o $(LIBDIR)weightwrite.o $(LIBDIR)xlibread.o
	$(F77) $(LFLAGS) -o vlook.x vlook.o binners.o binset.o chiprob.o cmatrix.o comparewrite.o comparer.o entropy.o filter.o fitherm.o galreadm.o getfwhm.o getlowright.o irivibin.o iterread.o iterwrite.o linpack.o minmax.o moment.o phaseread.o qual.o spear.o smcoarse.o spacem.o spline.o summ.o totmlread.o vaniread.o vdataread.o velbin.o velmtod.o vlibread.o weightread.o weightwrite.o xlibread.o $(PGPLOT)

binset.o cmatrix.o comparewrite.o comparer.o entropy.o filter.o galreadm.o getfwhm.o getlowright.o irivibin.o iterread.o iterwrite.o minmax.o vlook.o phaseread.o qual.o spear.o smcoarse.o spacem.o spline.o summ.o totmlread.o vaniread.o vdataread.o velbin.o velmtod.o vlibread.o weightread.o weightwrite.o xlibread.o : moddefs.h

#------------------------------------------------------------------------------#
# makefile for plot.libsos

plot.libsos : $(LIBDIR)plot.libsos.o $(LIBDIR)binners.o $(LIBDIR)binset.o $(LIBDIR)galaxyread.o $(LIBDIR)pgs.o $(LIBDIR)totmlread.o
	 $(F77) $(LFLAGS) -o plot.libsos.x plot.libsos.o binners.o binset.o galaxyread.o pgs.o totmlread.o $(PGPLOT) $(NUMREC)

plot.libsos.o binners.o binset.o dataread.o galaxyread.o totmlread.o : libdefs.h

# makefile for sos2vor

sos2vor : $(LIBDIR)sos2vor.o $(LIBDIR)binners.o $(LIBDIR)binset.o $(LIBDIR)dataread.o $(LIBDIR)density.o $(LIBDIR)force.o $(LIBDIR)forcefix.o $(LIBDIR)galaxyread.o $(LIBDIR)pl.o $(LIBDIR)potential.o $(LIBDIR)ran1.o $(LIBDIR)sort2.o $(LIBDIR)spacel.o $(LIBDIR)tablesread.o $(LIBDIR)haloread.o $(LIBDIR)halodens.o
	$(F77) $(LFLAGS) -o sos2vor.x sos2vor.o binners.o binset.o dataread.o density.o force.o forcefix.o galaxyread.o pl.o potential.o ran1.o sort2.o spacel.o tablesread.o haloread.o halodens.o $(NUMREC)

sos2vor.o binners.o binset.o dataread.o density.o force.o forcefix.o galaxyread.o potential.o spacel.o tablesread.o haloread.o halodens.o : libdefs.h

# makefile for vor2vphase

vor2vphase : $(LIBDIR)vor2vphase.o $(LIBDIR)binners.o $(LIBDIR)binset.o $(LIBDIR)density.o $(LIBDIR)ekin.o $(LIBDIR)energy.o $(LIBDIR)force.o $(LIBDIR)forcefix.o $(LIBDIR)galaxyread.o $(LIBDIR)potential.o $(LIBDIR)sort2.o $(LIBDIR)spacel.o $(LIBDIR)tablesread.o $(LIBDIR)totmlread.o $(LIBDIR)haloread.o $(LIBDIR)pl.o $(LIBDIR)dataread.o
	$(F77) $(LFLAGS) -o vor2vphase.x vor2vphase.o binners.o binset.o density.o ekin.o energy.o force.o forcefix.o galaxyread.o potential.o sort2.o spacel.o tablesread.o totmlread.o haloread.o pl.o dataread.o $(NUMREC)

vor2vphase.o binners.o binset.o density.o ekin.o energy.o force.o forcefix.o galaxyread.o potential.o spacel.o tablesread.o totmlread.o halo read.o pl.o dataread.o : libdefs.h

# make for seecor

seecor : $(LIBDIR)seecor.o $(LIBDIR)binners.o $(LIBDIR)binset.o $(LIBDIR)convolve.o $(LIBDIR)convolve2.o $(LIBDIR)galreadm.o $(LIBDIR)getsum.o $(LIBDIR)irivibin.o $(LIBDIR)velbin.o
	$(F77) $(LFLAGS) -o seecor.x seecor.o binners.o binset.o convolve.o convolve2.o galreadm.o getsum.o irivibin.o velbin.o $(HOSTLIBS)

seecor.o binners.o binset.o convolve.o convolve2.o galreadm.o getsum.o irivibin.o velbin.o : moddefs.h

# make for seecor2

seecor2 : $(LIBDIR)seecor2.o $(LIBDIR)binners.o $(LIBDIR)binset.o $(LIBDIR)convolve.o $(LIBDIR)convolve2.o $(LIBDIR)galreadm.o $(LIBDIR)getsum.o $(LIBDIR)irivibin.o $(LIBDIR)velbin.o
	$(F77) $(LFLAGS) -o seecor2.x seecor2.o binners.o binset.o convolve.o convolve2.o galreadm.o getsum.o irivibin.o velbin.o $(HOSTLIBS)

seecor2.o binners.o binset.o convolve.o convolve2.o galreadm.o getsum.o irivibin.o velbin.o : moddefs.h

plotall : $(LIBDIR)plotall.o
	$(F77) $(LFLAGS) -o plotall plotall.o $(HOSTLIBS) $(PGPLOT) ; rm *.o
plotall2 : $(LIBDIR)plotall2.o
	$(F77) $(LFLAGS) -o plotall2 plotall2.o $(HOSTLIBS) $(PGPLOT) ; rm *.o
porbit : $(LIBDIR)porbit.o
	$(F77) $(LFLAGS) -o porbit.x porbit.o $(HOSTLIBS) $(PGPLOT) ; rm *.o
porbit2 : $(LIBDIR)porbit2.o
	$(F77) $(LFLAGS) -o porbit2.x porbit2.o $(HOSTLIBS) $(PGPLOT) ; rm *.o

gden : $(LIBDIR)gden.o $(LIBDIR)halodens.o $(LIBDIR)halowrite.o
	$(F77) $(LFLAGS) -o gden gden.o halodens.o halowrite.o $(HOSTLIBS) $(QUEST)
gden.o : bothdefs.h

gdenconv : $(LIBDIR)gdenconv.o
	$(F77) $(LFLAGS) -o gdenconv gdenconv.o $(HOSTLIBS)

plotbin : $(LIBDIR)plotbin.o
	$(F77) $(LFLAGS) -o plotbin plotbin.o $(HOSTLIBS) $(PGPLOT) ; rm *.o
plotbin3 : $(LIBDIR)plotbin3.o
	$(F77) $(LFLAGS) -o plotbin3 plotbin3.o $(HOSTLIBS) $(PGPLOT) ; rm *.o
