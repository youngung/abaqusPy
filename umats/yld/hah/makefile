## make file to compile hah


## global variables
OBJDIR=objs
OBJS=bauschinger.o bauschinger_lib.o crossh.o hah.o \
	hah_lib.o hah_update.o latent.o


## Compiler options
CMP=gfortran
FLAGS=-Wall -g -std=gnu -fcheck=all -fstack-arrays


# targets below


######################################################################
# executables

hah_test: bauschinger_lib.o bauschinger.o hah_test.o
	$(CMP) $(FLAGS) bauschinger_lib.o bauschinger.o hah_test.o  -o hah_test


######################################################################
# object files


hah_test.o: hah_test.f
	$(CMP) $(FLAGS) -c hah_test.f


bauschinger_lib.o: bauschinger_lib.f
	$(CMP) $(FLAGS) -c bauschinger_lib.f



bauschinger.o: bauschinger.f
	$(CMP) $(FLAGS) -c bauschinger.f




# hah.o: hah.f
# 	$(CMP) $(FLAGS) -c hah.f

# hah_lib.o: hah_lib.f
# 	$(CMP) $(FLAGS) -c hah_lib.f





.PHONY : clean
clean :
	-rm $(OBJS)