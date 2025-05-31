# Compiler and flags
CC = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

# Output executable name
TARGET = lifofifo_interf2

# Source files
SRCS = lifofifo_abstract.cpp Datastruct.cpp Queue.cpp Stack.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Header files
HDRS = Datastruct.h Queue.h Stack.h Struct.h

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
