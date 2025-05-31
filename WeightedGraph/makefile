# Compiler and flags
CC = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

# Source files
SRCS = Edge.cpp Knot.cpp WeightGraph.cpp scf2.cpp
OBJS = $(SRCS:.cpp=.o)

# Executable name
TARGET = scf2

# Dependencies
DEPS = Edge.h Knot.h WeightGraph.h

# Build the executable
$(TARGET): $(OBJS)
    $(CC) $(CXXFLAGS) $(OBJS) -o $(TARGET)

# Compile .cpp files to .o files
%.o: %.cpp $(DEPS)
    $(CC) $(CXXFLAGS) -c $< -o $@

# Clean command
clean:
    rm -f $(OBJS) $(TARGET)

# Run the program
run: $(TARGET)
    ./$(TARGET)

# Rebuild everything
rebuild: clean $(TARGET)

.PHONY: clean run rebuild
