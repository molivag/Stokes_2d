#========== Definicion de variables ==========
# Actualizado el 15/07/2021
#	compiler
FC = ifort

#OBJ_DIR = ./obj	#FOR OBJECTS FILES
#BIN_DIR = ./bin  #FOR EXECUTABLE
#SRC_DIR = .

#	compiler flags
#	standard
#CFLAGS = -std=f2008ts
#	debugger option
CFLAGS += -g
#	warning flags
CFLAGS += -warn all
#	optimization flags
CFLAGS += -O0 -heap-arrays
#	error finding options
CFLAGS += -traceback -check all -CB -fp-stack-check

#	source files
SRCS = mod_parameters mod_elemKinds mod_library stokes

OBJS = $(SRCS:=.o)

#	executable 
MAIN = stokes.exe
#========== Fin variables ===========

#	compile project
all : $(MAIN)
#	@echo Compiling files . . . . .
#	@echo Making objects  . . . . . 
#	@echo Building an executable . . . . . 
	@echo ' '
	@echo '======================'
	@echo Compilation completed . . . . .
	@echo '======================'

$(MAIN) : $(OBJS)
	@$(FC) $(CFLAGS) -g -mkl -O0 $(OBJS) -o $(MAIN) 
# @$(FC) $(CFLAGS) $(OBJS) -o $(MAIN)
#These three flags are needed to correct execution of the LAPACK library
#If the libray is comented (not used) in the code, must them desactivated the flegs

.SUFFIXES : .f90 .o
#.o.f90 :Dos opciones, cual sera la correcta?

.f90.o :
	@$(FC) $(CFLAGS) -c $<

#	Regla ficticia, es decir que no tiene dependencias (phony rules)
clean :
	@$(RM) *.o *.mod *.exe $(MAIN)
#	clean no tiene dependencias pero si objetivos
	@echo '||'
	@echo ' - Everything is clean'
	@echo '||'
