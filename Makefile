CXX ?= g++
CPPFLAGS ?= -Isrc
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic -fopenmp
LDFLAGS ?= -fopenmp

TARGET := bin/2D-DWE
SOURCES := src/scheme_inertial.cpp src/main.cpp

ifeq ($(WITH_GDAL),1)
CPPFLAGS += -DDWE_WITH_GDAL $(shell gdal-config --cflags)
SOURCES += src/io_gdal.cpp
LDLIBS += $(shell gdal-config --libs)
endif

OBJECTS := $(patsubst src/%.cpp,build/obj/%.o,$(SOURCES))
DEPFILES := $(OBJECTS:.o=.d)

.PHONY: all clean test test-ea

all: $(TARGET)

$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

build/obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

test: $(TARGET)
	bash tests/verify.sh $(TARGET)

test-ea: $(TARGET)
	bash tests/verify_ea_test4.sh $(TARGET)

clean:
	rm -rf build bin

-include $(DEPFILES)
