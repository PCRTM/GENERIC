#================================================================================
#
# Makefile for PCRTM library build
#
#================================================================================

# The library name
LIBRARY_TAG = PCRTM
LIBRARY_NAME = lib$(LIBRARY_TAG).a

# The library directories
SOURCE_DIR = SRC
LIBRARY_DIR = lib
INCLUDE_DIR = include


# -------------------
# Compilation targets
# -------------------
# Default compiler build
all:
	cd $(SOURCE_DIR); make; cd ..

intel:
	cd $(SOURCE_DIR); make intel; cd ..

intel_debug:
	cd $(SOURCE_DIR); make intel_debug; cd ..

gfortran:
	cd $(SOURCE_DIR); make gfortran; cd ..

gfortran_debug:
	cd $(SOURCE_DIR); make gfortran_debug; cd ..


# ---------------------------------------------------------------------
# Install the library and include files in their respective directories
# ---------------------------------------------------------------------
install: install_lib install_include

install_lib:
	@if [ ! -d $(LIBRARY_DIR) ]; then \
	  mkdir $(LIBRARY_DIR); \
	fi; \
	cd $(SOURCE_DIR); \
	if [ -f $(LIBRARY_NAME) ]; then \
	  mv $(LIBRARY_NAME) ../$(LIBRARY_DIR); \
	fi ; \
	cd ..


install_include:
	@if [ ! -d $(INCLUDE_DIR) ]; then \
	  mkdir $(INCLUDE_DIR); \
	fi; \
	cd $(SOURCE_DIR); \
	for FILE in `ls -1 | egrep "\.(mod|MOD|m|M)$$"`; do \
	  if [ -f $$FILE ]; then \
	    mv -f $$FILE ../$(INCLUDE_DIR); \
	  fi ; \
	done ; \
	for FILE in `ls -1 | egrep "\.(o|O)$$"`; do \
	  if [ -f $$FILE ]; then \
	    rm -f $$FILE ; \
	  fi ; \
	done ; \
	cd ..

# --------
# Clean up
# --------
# Just the source directory
clean:
	@cd $(SOURCE_DIR); make clean; cd ..


# Everything
distclean: clean clean_lib clean_include

clean_lib:
	rm -r $(LIBRARY_DIR)

clean_include:
	rm -r $(INCLUDE_DIR)

