# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, MFEM Elasticity PA contributors.
.DEFAULT_GOAL := all

# MFEM configuration
MFEM_INSTALL_DIR ?= /usr/local
ifneq ($(filter-out clean distclean,$(or $(MAKECMDGOALS),all)),)
include $(MFEM_INSTALL_DIR)/share/mfem/config.mk
endif
MFEM_INC_DIR = $(MFEM_INSTALL_DIR)/include
MFEM_LIB_DIR = $(MFEM_INSTALL_DIR)/lib

# Compiler options
PLATFORM ?= generic
MY_FLAGS ?= -O3

ifeq ($(PLATFORM),generic)
    PLATFORM_FLAGS =
else ifeq ($(PLATFORM),amd)
    PLATFORM_FLAGS = -march=znver3 -mavx2 -mfma -mprefetchwt1
else ifeq ($(PLATFORM),arm)
    PLATFORM_FLAGS = -march=armv8.2-a+sve
else
    $(error PLATFORM must be generic, amd, or arm)
endif

# Source files
BUILD_DIR ?= build
OUTPUT_DIR = $(BUILD_DIR)/$(PLATFORM)
LIBRARY = $(OUTPUT_DIR)/libmfem_elasticity_pa.a
SOURCES = src/integ/MyElasticityIntegrator.cpp \
          src/kernels/MyElasticityKernels.cpp
OBJECTS = $(addprefix $(OUTPUT_DIR)/,$(SOURCES:.cpp=.o))
EXECUTABLE = $(OUTPUT_DIR)/myex2
DRIVER_OBJECT = $(OUTPUT_DIR)/myex2.o
DEPS = $(OBJECTS:.o=.d) $(DRIVER_OBJECT:.o=.d)

# Build
.PHONY: all myex2 clean distclean
all: $(LIBRARY)

ifeq ($(MFEM_USE_MPI),YES)
all: myex2
myex2: $(EXECUTABLE)

$(EXECUTABLE): $(DRIVER_OBJECT) $(LIBRARY)
	$(MFEM_CXX) $(MFEM_CXXFLAGS) $(MY_FLAGS) $(PLATFORM_FLAGS) $^ $(MFEM_LIBS) -o $@
else
myex2:
	$(error myex2 requires an MPI-enabled MFEM installation)
endif

$(LIBRARY): $(OBJECTS)
	$(AR) rcs $@ $^

$(OUTPUT_DIR)/%.o: %.cpp makefile
	mkdir -p $(dir $@)
	$(MFEM_CXX) $(MFEM_FLAGS) $(MY_FLAGS) $(PLATFORM_FLAGS) -MMD -MP -c $< -o $@

# Clean
clean:
	rm -rf -- "$(OUTPUT_DIR)"
	rm -f -- mesh.[0-9]* sol.[0-9]*

distclean:
	rm -rf -- "$(BUILD_DIR)/generic" "$(BUILD_DIR)/amd" "$(BUILD_DIR)/arm"
	rm -f -- mesh.[0-9]* sol.[0-9]*
	@rmdir -- "$(BUILD_DIR)" 2>/dev/null || true

-include $(DEPS)
