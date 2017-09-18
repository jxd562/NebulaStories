# make file for the box model of the solar nebula for Oxygen MIF

########################
# Setting definitions
########################
#FC = /opt/intel/fce/9.0/bin/ifort
#FC = /usr/local/intel-2013.3/bin/ifort
FC = /usr/local/bin/gfortran


FCFLAG = -r8 -i8 -align all -i-static
FCFLAG =
OUTPUT = oxiclean
PDOTO =  OXYGEN.o

# Main file
PPHOT = oxygen

OBPATH = $(PDOTO)/$(PPHOT).o

################################
# Loading compiled subroutines
################################
$(OUTPUT) : $(OBPATH)
	$(FC) $(FCFLAG) -o $(OUTPUT) $(OBPATH)
	
################################
# Object compile lines
################################
$(PDOTO)/$(PPHOT).o : $(PPHOT).f
	$(FC) $(FCFLAG) -c $(PPHOT).f
	\mv $(PPHOT).o $(PDOTO)
