# Uncomment on ieng6
# export PUB= /home/linux/ieng6/cs160f/public
# Uncomment on Bang
 export PUB=/share/class/public/cse160-fa13

# If you want to compile with MPI enabled,
# uncomment this line

# mpi = 1
include $(PUB)/Arch/arch.gnu-4.7_c++11.generic

#
# Add symbol table information for gdb
#
ifeq ($(dbg), 1)
        CFLAGS += -g
        LDFLAGS += -g
        C++FLAGS += -g
endif   

#
# Add debugging output controlled via
# conditional compilation with the DEBUG macro
# Add symbol table information for gdb/cachegrind
ifeq ($(debug), 1)
        CFLAGS += -g -DDEBUG
        LDFLAGS += -g
        C++FLAGS += -g -DDEBUG
endif   

app:		apf

OBJECTS = apf.o solve.o Plotting.o cmdLine.o Report.o
ifneq ($(mpi),1)
OBJECTS += Timer.o
endif

apf:	        $(OBJECTS) 
		$(C++LINK) $(LDFLAGS) -o $@ $(OBJECTS)  $(LDLIBS)

.PHONY: clean
clean:	
	$(RM) *.o apf;
	$(RM) core;
