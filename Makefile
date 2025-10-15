# This file is a part of DUDI-heliocentric, the Fortran-95 implementation 
# of the two-body model for the dynamics of dust ejected from an atmosphereless
# body moving around the Sun
# Version 1.0.1
# This is free software. You can use and redistribute it 
# under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
# If you do, please cite the following paper
# Anastasiia Ershova and Jürgen Schmidt, 
# Two-body model for the spatial distribution of dust ejected from
# an atmosphereless body, 2021, A&A, 650, A186 

# Author: Anastasiia Ershova
# E-mail: vveyzaa@gmail.com

# -------- compiler --------
FC      ?= gfortran
PYTHON  ?= python3
LDFLAGS ?=

# phony
.PHONY: all clean distclean clean-warnings help list \
        dudihc phaethon_dudihc select_method_dudihc \
        example_image phaethon select \
        run-example run-phaethon_dudihc run-select_method \
        debug debug-example
        
DEBUG_FLAGS = -O0 -g -fimplicit-none -fbacktrace -fcheck=all \
              -ffpe-trap=invalid,zero,overflow -finit-real=snan -finit-integer=-2147483648

debug:
	$(MAKE) clean
	$(MAKE) -B all FFLAGS='$(DEBUG_FLAGS)'

debug-example:
	$(MAKE) clean
	$(MAKE) -B $(BINDIR)/dudihc FFLAGS='$(DEBUG_FLAGS)'
	$(BINDIR)/dudihc


# your normal flags
FFLAGS ?= -O3 -fimplicit-none -Wall -Wno-tabs -Wno-unused-variable

# strict flags
STRICT_WARNINGS = -O0 -g -fimplicit-none -Wall -Wextra -Wconversion -Wsurprising \
                  -Warray-temporaries -Wcharacter-truncation -Wreal-q-constant \
                  -Wtarget-lifetime -Wimplicit-interface

# ---- guard against FC=f77 (disable on clean/list/help) ----
NEEDS_FC_GUARD := $(filter-out clean distclean help list,$(MAKECMDGOALS))
ifneq ($(NEEDS_FC_GUARD),)
  ifneq (,$(findstring f77,$(FC)))
    $(info FC='$(FC)' looks like Fortran 77; switching to gfortran)
    override FC := gfortran
  endif
endif

# -------- layout --------
SRCDIR := src
EXDIR  := examples
BINDIR := bin
MODDIR := build
RESDIR := results

# -------- core sources (ORDERED: providers before users) --------
CORE_SOURCES := \
  $(SRCDIR)/const.f90 \
  $(SRCDIR)/nan_utils.f90 \
  $(SRCDIR)/define_types.f90 \
  $(SRCDIR)/help.f90 \
  $(SRCDIR)/distributions_fun.f90 \
  $(SRCDIR)/twobody_fun.f90 \
  $(SRCDIR)/data_in.f90 \
  $(SRCDIR)/data_out.f90 \
  $(SRCDIR)/DUDIhc.f90 \
  $(SRCDIR)/phaethon_input.f90

# -------- example mains --------
EXAMPLE_SRC   ?= $(EXDIR)/example.f90
PHAETHON_SRC  ?= $(EXDIR)/phaethon.f90
SELECT_SRC    ?= $(EXDIR)/select_method.f90

# -------- plot scripts --------
EXAMPLE_PYSCRIPT  ?= scripts/show_image.py
PHAETHON_PYSCRIPT ?= scripts/plot_Fig10.py

# -------- objects --------
CORE_OBJS   := $(patsubst %.f90,$(MODDIR)/%.o,$(CORE_SOURCES))
EXAMPLE_OBJ := $(patsubst %.f90,$(MODDIR)/%.o,$(EXAMPLE_SRC))
PHAETHON_OBJ:= $(patsubst %.f90,$(MODDIR)/%.o,$(PHAETHON_SRC))
SELECT_OBJ  := $(patsubst %.f90,$(MODDIR)/%.o,$(SELECT_SRC))

# Default goal builds all binaries
all: $(BINDIR)/dudihc $(BINDIR)/phaethon_dudihc $(BINDIR)/select_method_dudihc

# ------------ link steps ------------
$(BINDIR)/dudihc: $(CORE_OBJS) $(EXAMPLE_OBJ) | $(BINDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) $^ -o $@ $(LDFLAGS)
dudihc: $(BINDIR)/dudihc

$(BINDIR)/phaethon_dudihc: $(CORE_OBJS) $(PHAETHON_OBJ) | $(BINDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) $^ -o $@ $(LDFLAGS)
phaethon_dudihc: $(BINDIR)/phaethon_dudihc

$(BINDIR)/select_method_dudihc: $(CORE_OBJS) $(SELECT_OBJ) | $(BINDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) $^ -o $@ $(LDFLAGS)
select_method_dudihc: $(BINDIR)/select_method_dudihc

# ------------ compile step (pattern rule) ------------
$(MODDIR)/%.o: %.f90 | $(MODDIR)
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) -c $< -o $@

# Suppress compare-reals ONLY for nan_utils (intentional x /= x)
$(MODDIR)/$(SRCDIR)/nan_utils.o: FFLAGS += -Wno-compare-reals

# ------------ pipelines: build -> run -> plot ------------
example_image: $(BINDIR)/dudihc | $(RESDIR)
	@echo ">>> Running dudihc"
	$(BINDIR)/dudihc
	@echo ">>> Plot: $(PYTHON) $(EXAMPLE_PYSCRIPT)"
	$(PYTHON) $(EXAMPLE_PYSCRIPT) || true

phaethon: $(BINDIR)/phaethon_dudihc | $(RESDIR)
	@echo ">>> Running phaethon_dudihc"
	$(BINDIR)/phaethon_dudihc
	@echo ">>> Plot: $(PYTHON) $(PHAETHON_PYSCRIPT)"
	$(PYTHON) $(PHAETHON_PYSCRIPT) || true

select: $(BINDIR)/select_method_dudihc | $(RESDIR)
	@echo ">>> Running select_method_dudihc"
	$(BINDIR)/select_method_dudihc

# run-only helpers
run-example:           $(BINDIR)/dudihc            ; $(BINDIR)/dudihc
run-phaethon_dudihc:   $(BINDIR)/phaethon_dudihc   ; $(BINDIR)/phaethon_dudihc
run-select_method:     $(BINDIR)/select_method_dudihc ; $(BINDIR)/select_method_dudihc

# ------------ utilities ------------
help:
	@echo "Build-only:  dudihc | phaethon_dudihc | select_method_dudihc"
	@echo "Pipelines:   example_image | phaethon | select   (build -> run -> plot)"
	@echo "List files:  make list"
	@echo "Clean:       make clean   /  make distclean"

list:
	@echo "Core:   $(CORE_SOURCES)"
	@echo "Mains:  EXAMPLE=$(EXAMPLE_SRC)  PHAETHON=$(PHAETHON_SRC)  SELECT=$(SELECT_SRC)"
	@echo "Plots:  EXAMPLE=$(EXAMPLE_PYSCRIPT)  PHAETHON=$(PHAETHON_PYSCRIPT)"

# ensure dirs
$(BINDIR) $(MODDIR) $(RESDIR):
	mkdir -p $@

# cleaning
clean:
	rm -rf $(MODDIR)/*
distclean:
	rm -rf $(BINDIR) $(RESDIR)/* $(MODDIR)/*

# Strict warnings sweep: clean + rebuild (compile only)
clean-warnings:
	$(MAKE) clean
	$(MAKE) -B all FFLAGS='$(STRICT_WARNINGS)'
