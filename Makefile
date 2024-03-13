# Compiler
CC = gcc

# Compiler flags
CPPFLAGS += -Ihtslib/ -Iinclude/
CFLAGS = -g -Wall -O2  -std=c99
LDFLAGS  += $(LIBS) -lz -lm -lpthread

# Source files
SRCS = src/error.c src/file_handler.c src/logger.c src/main.c src/misc_p.c src/misc.c src/mod.c src/thread.c src/utils.c

# Object files directory
BUILD_DIR = build

# Object files
OBJS = $(addprefix $(BUILD_DIR)/, $(notdir $(SRCS:.c=.o)))

# Executable name
BINARY = realfreq

# Phony targets
.PHONY: all clean

# Default target
all: $(BINARY)

# Address sanitizer
ifdef asan
	CFLAGS += -fsanitize=address -fno-omit-frame-pointer
	LDFLAGS += -fsanitize=address -fno-omit-frame-pointer
endif

# Compile source files into object files
$(BUILD_DIR)/%.o: src/%.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

# Link object files into the executable
$(BINARY): $(OBJS) htslib/libhts.a
	$(CC) $(CFLAGS) $(OBJS) htslib/libhts.a $(LDFLAGS) -o $@

# Check if htslib is installed
htslib/libhts.a:
	@if test -e $(BUILD_DIR)/lib/libhts.a; then \
		echo "htslib found at htslib/libhts.a"; \
	else \
		echo "htslib not found at htslib/libhts.a"; \
		echo "Please run 'scripts/install-hts.sh' first"; \
		exit 1; \
	fi

# Clean up intermediate and executable files
clean:
	rm -f $(OBJS) $(BINARY)

# Run tests
test: $(BINARY)
	./test/test.sh

memtest: $(BINARY)
	./test/test.sh mem

# Create build directory if it doesn't exist
$(shell mkdir -p $(BUILD_DIR))
