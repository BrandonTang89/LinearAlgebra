# Makefile to compile main.cpp

# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -g --std=c++23

# Target executable
TARGET = main

# Source files
SRCS = main.cpp Vec.cpp Matrix.cpp ConjugateGradient.cpp GivensRotations.cpp
# Object files
OBJS = $(SRCS:.cpp=.o)

# Default target
all: $(TARGET)

# Rule to link object files and create the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Rule to compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean