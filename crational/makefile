# Compiler and flags
CC = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

# Output executable name
TARGET = crational

# Source files
SRCS = main.cpp Crational.cpp CrationalExeption.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Header files
HDRS = Crational.h CrationalExeption.h

# All dependencies
DEPS = $(HDRS)

# Default target
all: $(TARGET)

# Rule to build the executable
$(TARGET): $(OBJS)
    $(CC) $(CXXFLAGS) $^ -o $@

# Rule to compile .cpp files into .o files
%.o: %.cpp $(DEPS)
    $(CC) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
    rm -f $(OBJS) $(TARGET)

# Rebuild rule
rebuild: clean all

# Help rule
help:
    @echo "Available targets:"
    @echo "  all     - Build the executable"
    @echo "  clean   - Remove object files and executable"
    @echo "  rebuild - Clean and rebuild everything"
    @echo "  help    - Show this help message"
