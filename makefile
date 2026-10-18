# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, MFEM Elasticity PA contributors.
.DEFAULT_GOAL := all

# MFEM configuration
MFEM_INSTALL_DIR ?= /usr/local
include $(MFEM_INSTALL_DIR)/share/mfem/config.mk
MFEM_INC_DIR = $(MFEM_INSTALL_DIR)/include

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
DEPS = $(OBJECTS:.o=.d)

# Build
.PHONY: all clean
all: $(LIBRARY)

$(LIBRARY): $(OBJECTS)
	$(AR) rcs $@ $^

$(OUTPUT_DIR)/%.o: %.cpp makefile
	mkdir -p $(dir $@)
	$(MFEM_CXX) $(MFEM_FLAGS) $(MY_FLAGS) $(PLATFORM_FLAGS) -MMD -MP -c $< -o $@

# Clean
clean:
	rm -f $(LIBRARY) $(OBJECTS) $(DEPS)

-include $(DEPS)
