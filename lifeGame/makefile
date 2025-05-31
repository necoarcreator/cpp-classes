# Compiler and flags
CC = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

# Source files
SRCS = Cell.cpp Game.cpp main.cpp ekzamen.cpp
OBJS = $(SRCS:.cpp=.o)

# Executable name
TARGET = ekzamen

# Dependencies
DEPS = Cell.h Game.h

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
