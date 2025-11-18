#!/usr/bin/make

#main building variables
DSRC    = .
DOBJ    = obj/
DMOD    = mod/
DEXE    = ./
LIBS    =
FC      = ifort
OPTSC   = -c -qopenmp -module mod
OPTSL   =  -oqpenmp -module mod
VPATH   = $(DSRC) $(DOBJ) $(DMOD)
MKDIRS  = $(DOBJ) $(DMOD) $(DEXE)
LCEXES  = $(shell echo $(EXES) | tr '[:upper:]' '[:lower:]')
EXESPO  = $(addsuffix .o,$(LCEXES))
EXESOBJ = $(addprefix $(DOBJ),$(EXESPO))

#auxiliary variables
COTEXT  = "Compiling $(<F)"
LITEXT  = "Assembling $@"

#building rules
$(DEXE)MAIN: $(MKDIRS) $(DOBJ)main.o \
	$(DOBJ)labor1.o \
	$(DOBJ)labor2.o \
	$(DOBJ)labor3.o \
	$(DOBJ)labors.o \
	$(DOBJ)lsupply.o \
	$(DOBJ)simulation.o \
	$(DOBJ)solveactivelife.o \
	$(DOBJ)solvefirstactive.o \
	$(DOBJ)solveinretirement.o \
	$(DOBJ)statistics.o \
	$(DOBJ)dogleg.o \
	$(DOBJ)dpmpar.o \
	$(DOBJ)enorm.o \
	$(DOBJ)fdjac1.o \
	$(DOBJ)hybrd.o \
	$(DOBJ)qform.o \
	$(DOBJ)qrfac.o \
	$(DOBJ)r1mpyq.o \
	$(DOBJ)r1updt.o
	@rm -f $(filter-out $(DOBJ)main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) MAIN

#compiling rules
$(DOBJ)bspline_kinds_module.o: ./bspline_kinds_module.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)bspline_sub_module.o: ./bspline_sub_module.f90 \
	$(DOBJ)bspline_kinds_module.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)glob0.o: ./glob0.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)goldensearch.o: ./GoldenSearch.f90 \
	$(DOBJ)policyfunctions.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)labor1.o: ./labor1.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)labor2.o: ./labor2.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)labor3.o: ./labor3.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)labors.o: ./labors.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)lsupply.o: ./lsupply.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o \
	$(DOBJ)hybrd_wrapper.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)main.o: ./main.f90 \
	$(DOBJ)utilities.o \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)tauchen.o \
	$(DOBJ)hybrd_wrapper.o \
	$(DOBJ)partest.o \
	$(DOBJ)solvelastactive.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)model_parameters.o: ./Model_Parameters.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)partest.o: ./partest.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o \
	$(DOBJ)bspline_sub_module.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)policyfunctions.o: ./PolicyFunctions.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)bspline_sub_module.o \
	$(DOBJ)utilities.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)simulation.o: ./Simulation.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)utilities.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)solveactivelife.o: ./SolveActiveLife.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o \
	$(DOBJ)bspline_sub_module.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)solvefirstactive.o: ./Solvefirstactive.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o \
	$(DOBJ)bspline_sub_module.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)solveinretirement.o: ./SolveInRetirement.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)utilities.o \
	$(DOBJ)bspline_sub_module.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)solvelastactive.o: ./Solvelastactive.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)glob0.o \
	$(DOBJ)utilities.o \
	$(DOBJ)bspline_sub_module.o \
	$(DOBJ)goldensearch.o \
	$(DOBJ)valuefunctions.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)statistics.o: ./Statistics.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)utilities.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)tauchen.o: ./tauchen.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)utilities.o: ./Utilities.f90 \
	$(DOBJ)model_parameters.o \
	$(DOBJ)glob0.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)valuefunctions.o: ./ValueFunctions.f90 \
	$(DOBJ)policyfunctions.o \
	$(DOBJ)bspline_sub_module.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)dogleg.o: ./minpack-hybrd/dogleg.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)dpmpar.o: ./minpack-hybrd/dpmpar.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)enorm.o: ./minpack-hybrd/enorm.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)fdjac1.o: ./minpack-hybrd/fdjac1.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)hybrd.o: ./minpack-hybrd/hybrd.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)hybrd_wrapper.o: ./minpack-hybrd/hybrd_wrapper.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)qform.o: ./minpack-hybrd/qform.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)qrfac.o: ./minpack-hybrd/qrfac.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)r1mpyq.o: ./minpack-hybrd/r1mpyq.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)r1updt.o: ./minpack-hybrd/r1updt.f
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

#phony auxiliary rules
.PHONY : $(MKDIRS)
$(MKDIRS):
	@mkdir -p $@
.PHONY : cleanobj
cleanobj:
	@echo deleting objects
	@rm -fr $(DOBJ)
.PHONY : cleanmod
cleanmod:
	@echo deleting mods
	@rm -fr $(DMOD)
.PHONY : cleanexe
cleanexe:
	@echo deleting exes
	@rm -f $(addprefix $(DEXE),$(EXES))
.PHONY : clean
clean: cleanobj cleanmod
.PHONY : cleanall
cleanall: clean cleanexe
