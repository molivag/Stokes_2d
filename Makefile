
#========== Definicion de variables ==========
#	compiler
FC = ifort # /opt/intel/oneapi/compiler/2021.2.0/mac/bin/intel64/ifort

#OBJ_DIR = ./obj	#FOR OBJECTS FILES
#BIN_DIR = ./bin  #FOR EXECUTABLE
#SRC_DIR = .

#	compiler flags
#	standard
#CFLAGS = -std=f2008ts
#	warning flags
CFLAGS += -warn all
#	debugging options
CFLAGS += -check all -CB -traceback 
#	GDB flag
CFLAGS += -g
#	source files
SRCS = mod_library stokes

OBJS = $(SRCS:=.o)

#	executable 
MAIN = stokes
#========== Fin variables ===========

#	compile project
all : $(MAIN)
#	@echo Compiling files . . . . .
#	@echo Making objects  . . . . . 
#	@echo Building an executable . . . . . 
	@echo '======================'
	@echo Compilation completed . . . . .
	@echo '======================'

$(MAIN) : $(OBJS)
	@$(FC) $(CFLAGS) -mkl $(OBJS) -o $(MAIN)

.SUFFIXES : .f90 .o
#.o.f90 :Dos opciones, cual sera la correcta?

.f90.o :
	@$(FC) $(CFLAGS) -c $<

#	Regla ficticia, es decir que no tiene dependencias (phony rules)
clean :
	@$(RM) *.o *.mod $(MAIN)
#	clean no tiene dependencias pero si objetivos
	@echo Everything is clean
