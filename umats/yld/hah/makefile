## make file to compile hah

## global variables
OBJDIR=objs

## Compiler options
CMP=gfortran
FLAGS=-g -std=gnu -fcheck=all -fstack-arrays


# targets below


######################################################################
# executables

hah_test: hah_test.o bauschinger_lib.o hah_lib.o bauschinger.o hah.o hah_yieldsurface.o
	$(CMP) $(FLAGS) hah_test.o bauschinger_lib.o hah_lib.o bauschinger.o hah.o hah_yieldsurface.o -o hah_test

######################################################################
# object files


hah_test.o: hah_test.f
	$(CMP) $(FLAGS) -c hah_test.f

bauschinger_lib.o: bauschinger_lib.f
	$(CMP) $(FLAGS) -c bauschinger_lib.f

hah_lib.o: hah_lib.f
	$(CMP) $(FLAGS) -c hah_lib.f

bauschinger.o: bauschinger.f
	$(CMP) $(FLAGS) -c bauschinger.f

hah.o: hah.f
	$(CMP) $(FLAGS) -c hah.f

hah_yieldsurface.o: hah_yieldsurface.f
	$(CMP) $(FLAGS) -c hah_yieldsurface.f



.PHONY : clean
clean :
	-rm *.o