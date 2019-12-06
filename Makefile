## make file for yld, yld testing

FC=gfortran
OBJDIR:=objdir

# -fno-automatic

FLAGS=-g -fcheck=all -Waliasing -Warray-bounds \
	-fbacktrace -fstack-arrays -finit-local-zero -Winteger-division\
	-Wfatal-errors # -Werror -Wall
objects_hah_test=$(addprefix $(OBJDIR)/,hah_test.o bauschinger_lib.o \
	hah_lib.o bauschinger.o hah.o \
	hah_yieldsurface.o yld.o vm.o cnv.o hill48.o yld2000_2d.o\
	hah_update.o algb.o dev.o lib_write.o is.o lib.o yld_lib.o\
	latent.o microd.o crossh.o deriv_lib.o elast.o)
objects_hahd_test=$(addprefix $(OBJDIR)/,hah_test.o bauschinger_lib.o \
	hah_lib.o bauschinger.o hah.o \
	hah_yieldsurface.o yld.o vm.o cnv.o hill48.o yld2000_2d.o\
	hah_update.o algb.o dev.o lib_write.o is.o lib.o yld_lib.o\
	latent.o microd.o deriv_lib.o elast.o crossh.o)
objects_yld_test=$(addprefix $(OBJDIR)/,yld.o yld_uten_test.o yld2000_2d.o bauschinger.o \
	hah.o hah_update.o dev.o algb.o cnv.o vm.o hill48.o \
	hah_yieldsurface.o hah_lib.o lib_write.o is.o lib.o latent.o deriv_lib.o \
	crossh.o elast.o microd.o)

objs+=$(objects_hah_test)
objs+=$(objects_hahd_test)
objs+=$(objects_yld_test)

######################################################################
# Fortran executables
hah_test: $(objects_hah_test)
	$(FC) $(FLAGS) -O3 $(objects_hah_test) -o hah_test
hahd_test: $(objects_hahd_test)
	$(FC) $(FLAGS) -O3 $(objects_hahd_test) -o hahd_test
yld_test: $(objects_yld_test)
	$(FC) $(FLAGS) $(objects_yld_test) -o yld_test
######################################################################

$(objects_hah_test): | $(OBJDIR)
$(objects_hahd_test): | $(OBJDIR)
$(objects_yld_test): | $(OBJDIR)

$(OBJDIR):
	mkdir $(OBJDIR)


$(OBJDIR)/microd.o: umats/yld/hah/microd.f
	$(FC) $(FLAGS) -c umats/yld/hah/microd.f  -o $(OBJDIR)/microd.o
$(OBJDIR)/yld_uten_test.o: umats/yld/test.f
	$(FC) $(FLAGS) -c umats/yld/test.f -o $(OBJDIR)/yld_uten_test.o
$(OBJDIR)/yld.o: umats/yld/yld.f
	$(FC) $(FLAGS) -c umats/yld/yld.f -o $(OBJDIR)/yld.o
$(OBJDIR)/vm.o:  umats/yld/vm.f
	$(FC) $(FLAGS) -c umats/yld/vm.f -o $(OBJDIR)/vm.o
$(OBJDIR)/hill48.o:  umats/yld/hill48.f
	$(FC) $(FLAGS) -c umats/yld/hill48.f -o $(OBJDIR)/hill48.o
$(OBJDIR)/yld2000_2d.o:  umats/yld/yld2000_2d.f
	$(FC) $(FLAGS) -c umats/yld/yld2000_2d.f -o $(OBJDIR)/yld2000_2d.o
$(OBJDIR)/hah_update.o: umats/yld/hah/hah_update.f
	$(FC) $(FLAGS) -c umats/yld/hah/hah_update.f -o $(OBJDIR)/hah_update.o
$(OBJDIR)/dev.o: umats/lib/dev.f
	$(FC) $(FLAGS) -c umats/lib/dev.f -o $(OBJDIR)/dev.o
$(OBJDIR)/algb.o: umats/lib/algb.f
	$(FC) $(FLAGS) -c umats/lib/algb.f -o $(OBJDIR)/algb.o
$(OBJDIR)/cnv.o: umats/lib/cnv.f
	$(FC) $(FLAGS) -c umats/lib/cnv.f -o $(OBJDIR)/cnv.o
$(OBJDIR)/hah_test.o: umats/yld/hah/hah_test.f
	$(FC) $(FLAGS) -c umats/yld/hah/hah_test.f -o $(OBJDIR)/hah_test.o
$(OBJDIR)/bauschinger_lib.o: umats/yld/hah/bauschinger_lib.f
	$(FC) $(FLAGS) -c umats/yld/hah/bauschinger_lib.f -o $(OBJDIR)/bauschinger_lib.o
$(OBJDIR)/bauschinger.o: umats/yld/hah/bauschinger.f
	$(FC) $(FLAGS) -c umats/yld/hah/bauschinger.f -o $(OBJDIR)/bauschinger.o
$(OBJDIR)/hah.o: umats/yld/hah/hah.f
	$(FC) $(FLAGS) -c umats/yld/hah/hah.f -o $(OBJDIR)/hah.o
$(OBJDIR)/hah_lib.o: umats/yld/hah/hah_lib.f
	$(FC) $(FLAGS) -c umats/yld/hah/hah_lib.f -o $(OBJDIR)/hah_lib.o
$(OBJDIR)/hah_yieldsurface.o: umats/yld/hah/hah_yieldsurface.f
	$(FC) $(FLAGS) -c umats/yld/hah/hah_yieldsurface.f -o $(OBJDIR)/hah_yieldsurface.o
$(OBJDIR)/lib_write.o: umats/lib/lib_write.f
	$(FC) $(FLAGS) -c umats/lib/lib_write.f -o $(OBJDIR)/lib_write.o
$(OBJDIR)/is.o: umats/lib/is.f
	$(FC) $(FLAGS) -c umats/lib/is.f -o $(OBJDIR)/is.o
$(OBJDIR)/lib.o: umats/lib/lib.f
	$(FC) $(FLAGS) -c umats/lib/lib.f -o $(OBJDIR)/lib.o
$(OBJDIR)/yld_lib.o: umats/yld/yld_lib.f
	$(FC) $(FLAGS) -c umats/yld/yld_lib.f -o $(OBJDIR)/yld_lib.o
$(OBJDIR)/latent.o: umats/yld/hah/latent.f
	$(FC) $(FLAGS) -c umats/yld/hah/latent.f -o $(OBJDIR)/latent.o
$(OBJDIR)/crossh.o: umats/yld/hah/crossh.f
	$(FC) $(FLAGS) -c umats/yld/hah/crossh.f -o $(OBJDIR)/crossh.o
$(OBJDIR)/deriv_lib.o: umats/yld/hah/deriv_lib.f
	$(FC) $(FLAGS) -c umats/yld/hah/deriv_lib.f -o $(OBJDIR)/deriv_lib.o
$(OBJDIR)/elast.o: umats/lib/elast.f
	$(FC) $(FLAGS) -c umats/lib/elast.f -o $(OBJDIR)/elast.o


.PHONY: all clean
all: hah_test hahd_test yld_test
clean:
	-rm hah_test yld_test hahd_test $(OBJDIR)/*.o
