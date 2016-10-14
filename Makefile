## make file for yld, yld testing

CMP=gfortran
FLAGS=-g -fno-automatic -Wall -fcheck=all -Waliasing -Warray-bounds -fbacktrace -fstack-arrays

objects_hah_test=hah_test.o bauschinger_lib.o \
	hah_lib.o bauschinger.o hah.o \
	hah_yieldsurface.o yld.o vm.o cnv.o hill48.o yld2000_2d.o\
	hah_update.o algb.o dev.o lib_write.o is.o lib.o yld_lib.o\
	latent.o
objects_yld_test=yld.o yld_uten_test.o yld2000_2d.o bauschinger.o \
	hah.o hah_update.o dev.o algb.o cnv.o vm.o hill48.o \
	hah_yieldsurface.o hah_lib.o lib_write.o is.o lib.o latent.o

######################################################################
# Fortran executables
hah_test: $(objects_hah_test)
	$(CMP) $(FLAGS) $(objects_hah_test) -o hah_test

yld_test: $(objects_yld_test)
	$(CMP) $(FLAGS) $(objects_yld_test) -o yld_test
######################################################################

yld_uten_test.o: umats/yld/test.f
	$(CMP) $(FLAGS) -c umats/yld/test.f -o yld_uten_test.o
yld.o: umats/yld/yld.f
	$(CMP) $(FLAGS) -c umats/yld/yld.f
vm.o:  umats/yld/vm.f
	$(CMP) $(FLAGS) -c umats/yld/vm.f
hill48.o:  umats/yld/hill48.f
	$(CMP) $(FLAGS) -c umats/yld/hill48.f
yld2000_2d.o:  umats/yld/yld2000_2d.f
	$(CMP) $(FLAGS) -c umats/yld/yld2000_2d.f
hah_update.o: umats/yld/hah/hah_update.f
	$(CMP) $(FLAGS) -c umats/yld/hah/hah_update.f
dev.o: umats/lib/dev.f
	$(CMP) $(FLAGS) -c umats/lib/dev.f
algb.o: umats/lib/algb.f
	$(CMP) $(FLAGS) -c umats/lib/algb.f
cnv.o: umats/lib/cnv.f
	$(CMP) $(FLAGS) -c umats/lib/cnv.f
hah_test.o: umats/yld/hah/hah_test.f
	$(CMP) $(FLAGS) -c umats/yld/hah/hah_test.f
bauschinger_lib.o: umats/yld/hah/bauschinger_lib.f
	$(CMP) $(FLAGS) -c umats/yld/hah/bauschinger_lib.f
bauschinger.o: umats/yld/hah/bauschinger.f
	$(CMP) $(FLAGS) -c umats/yld/hah/bauschinger.f
hah.o: umats/yld/hah/hah.f
	$(CMP) $(FLAGS) -c umats/yld/hah/hah.f
hah_lib.o: umats/yld/hah/hah_lib.f
	$(CMP) $(FLAGS) -c umats/yld/hah/hah_lib.f
hah_yieldsurface.o: umats/yld/hah/hah_yieldsurface.f
	$(CMP) $(FLAGS) -c umats/yld/hah/hah_yieldsurface.f
lib_write.o: umats/lib/lib_write.f
	$(CMP) $(FLAGS) -c umats/lib/lib_write.f
is.o: umats/lib/is.f
	$(CMP) $(FLAGS) -c umats/lib/is.f
lib.o: umats/lib/lib.f
	$(CMP) $(FLAGS) -c umats/lib/lib.f
yld_lib.o: umats/yld/yld_lib.f
	$(CMP) $(FLAGS) -c umats/yld/yld_lib.f
latent.o: umats/yld/hah/latent.f
	$(CMP) $(FLAGS) -c umats/yld/hah/latent.f


.PHONY: all clean
all: hah_test yld_test
clean :
	-rm hah_test yld_test *.o
