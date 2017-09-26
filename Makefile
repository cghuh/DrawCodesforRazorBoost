#

PROGRAM 	= StackedPlot
PROGRAM1 	= StackedPlot_2D
CC		= g++
LD		= g++
CXXFLAGS        = -g --exceptions -I$(ROOTSYS)/include -std=c++11
LDFLAGS         = $(LIBS) $(GLIBS) 
#-B$(LINKTIME_DIR)/usr/X11R6/lib -B$(LINKTIME_DIR)/lib -B$(LINKTIME_DIR)/usr/lib -Wl,-rpath-link -Wl,$(LINKTIME_DIR)/usr/X11R6/lib -Wl,-rpath-link -Wl,$(LINKTIME_DIR)/lib -Wl,-rpath-link -Wl,$(LINKTIME_DIR)/usr/lib

ROOTLIBS	= $(shell root-config --libs)
ROOTGLIBS	= $(shell root-config --glibs)

LIBS            = $(ROOTLIBS) -lMinuit
GLIBS           = $(ROOTGLIBS) -L/usr/X11R6/lib -lXext #-lX11

OBJS		= StackedPlot.o
OBJS1		= StackedPlots_2D.o
INPUTS		= $(PWD)/StackedPlot.C
INPUTS1		= $(PWD)/StackedPlots_2D.C


all: $(OBJS)
	$(CC) $(OPT) $(LDFLAGS) -o $(PROGRAM) -g $(OBJS) $(LIBS)
#all: $(OBJS1)
	$(CC) $(OPT) $(LDFLAGS) -o $(PROGRAM1) -g $(OBJS1) $(LIBS)


ControlPlots_Macro.o: $(INPUTS)
	$(CC) $(CXXFLAGS) -c StackedPlot.C
StackedPlots_2D.o: $(INPUTS1)
	$(CC) $(CXXFLAGS) -c StackedPlots_2D.C


clean:
	-rm -f $(PROGRAM)
	-rm -f $(PROGRAM1)




run: $(PROGRAM) $(OBJS)
	gmake all
	./$(PROGRAM) 1
