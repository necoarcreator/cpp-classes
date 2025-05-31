# Compiler and flags
CC = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

# Source files
SRCS = ActiveArty.cpp Calculator.cpp Container.cpp Interface.cpp main.cpp
OBJS = $(SRCS:.cpp=.o)

# Executable name
TARGET = main

# Dependencies
DEPS = ActiveArty.h Calculator.h Container.h Interface.h

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
