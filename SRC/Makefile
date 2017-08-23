
# The library name
LIBRARY_TAG = PCRTM
LIBRARY_NAME = lib$(LIBRARY_TAG).a

HOME = ..

SOURCE1	= \
${HOME}/SRC/PCRTM_TYPE_KIND.f90  \
${HOME}/SRC/PCRTM_CONSTANTS.f90 \
${HOME}/SRC/PCRTM_MATH_UTILITY.f90 \
${HOME}/SRC/PCRTM_FILE_UTILITIES.f90  \
${HOME}/SRC/PCRTM_ATMOSPHERE_DEFINE.f90\
${HOME}/SRC/PCRTM_ATM_ABSORPTION_DEFINE.f90\
${HOME}/SRC/PCRTM_ATM_ABSORPTION_IO.f90 \
${HOME}/SRC/PCRTM_ATMOSPHERE_LAYER.f90\
${HOME}/SRC/PCRTM_SOLAR_DEFINE.f90\
${HOME}/SRC/PCRTM_SOLAR_PARAMETER.f90\
${HOME}/SRC/PCRTM_JACOBIAN.f90\
${HOME}/SRC/PCRTM_CLOUD_LUT_IO.f90\
${HOME}/SRC/PCRTM_CLOUD_DEFINE.f90\
${HOME}/SRC/PCRTM_CLOUD_LUT_INTERP.f90\
${HOME}/SRC/PCRTM_CALC_CLOUD.f90\
${HOME}/SRC/PCRTM_RT_SOLUTION_DEFINE.f90 \
${HOME}/SRC/PCRTM_ATM_ABSORPTION.f90\
${HOME}/SRC/PCRTM_PC_SOLUTION.f90 \
${HOME}/SRC/PCRTM_STRUCT_MOD.f90 \
${HOME}/SRC/PCRTM_TR_SOLUTION.f90 \
${HOME}/SRC/CLEAR_PCRTM.f90\
${HOME}/SRC/INIT_PCRTM.f90 \
${HOME}/SRC/PCRTM_CALC_RAD.f90 \
${HOME}/SRC/PCRTM_SOLAR_RT.f90\
${HOME}/SRC/PCRTM_FORWARD_MODEL.f90\



.SUFFIXES:	.o .f90

INTEL_OPTIM = " f90 = ifort "\
              " FF_FLAGS = -O3 \
                           -warn errors \
                           -free \
      	                   -assume byterecl" \

INTEL_DEBUG = " f90 = ifort " \
              " FF_FLAGS = -g \
                           -check bounds \
                           -traceback \
                           -warn errors \
                           -free \
                           -assume byterecl \
                           -fpstkchk \
                           -mp " \

GFORTRAN_OPTIM = " f90 = gfortran "\
                 " FF_FLAGS =-c -O3 \
                             -ffast-math \
                             -ffpe-trap=overflow,zero \
                             -ffree-form \
                             -fno-second-underscore \
                             -frecord-marker=4 \
                             -funroll-loops \
                             -static \
                             -Wall " \

GFORTRAN_DEBUG = " f90 = gfortran" \
                 " FF_FLAGS = -g \
                              -fbounds-check \
                              -ffpe-trap=overflow,zero \
                              -ffree-form \
                              -fno-second-underscore \
                              -frecord-marker=4 \
                              -ggdb \
	                      -static \
                              -Wall " \

OBJECTS1 = ${SOURCE1:.f90=.o}

.f90.o:	  	${SOURCE1}        
	${f90} ${FF_FLAGS} -c $<

intel: $(OBJECTS1)
	make -f Makefile library  $(INTEL_OPTIM)

gfortran : $(OBJECTS1)
	make -f Makefile library  $(GFORTRAN_OPTIM)

intel_debug: $(OBJECTS1)
	make -f Makefile library  $(INTEL_DEBUG)

gfortran_debug : $(OBJECTS1)
	make -f Makefile library  $(GFORTRAN_DEBUG)

clean:
	rm -f  ${OBJECTS1} *.mod *.o *.a

# Define build dependencies
include make.dependencies

# Build the library 
library: $(OBJECTS1)
	ar crvs $(LIBRARY_NAME) $(OBJECTS1)