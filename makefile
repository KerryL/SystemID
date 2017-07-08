# makefile (SystemID)
#
# Include the common definitions
include makefile.inc

# Name of the executable to compile and link
TARGET = systemID
TARGET_DEBUG = systemIDd

# Directories in which to search for source files
DIRS = \
	src/

# Source files
SRC = $(foreach dir, $(DIRS), $(wildcard $(dir)*.cpp))

# Object files
OBJS_DEBUG = $(addprefix $(OBJDIR_DEBUG),$(SRC:.cpp=.o))
OBJS_RELEASE = $(addprefix $(OBJDIR_RELEASE),$(SRC:.cpp=.o))

.PHONY: all debug clean

all: $(TARGET)
debug: $(TARGET_DEBUG)

$(TARGET): $(OBJS_RELEASE)
	$(MKDIR) $(BINDIR)
	$(CC) $(OBJS_RELEASE) $(LDFLAGS_RELEASE) -L$(LIBOUTDIR) $(addprefix -l,$(PSLIB)) -o $(BINDIR)$@

$(TARGET_DEBUG): $(OBJS_DEBUG)
	$(MKDIR) $(BINDIR)
	$(CC) $(OBJS_DEBUG) $(LDFLAGS_DEBUG) -L$(LIBOUTDIR) $(addprefix -l,$(PSLIB)) -o $(BINDIR)$@

$(OBJDIR_RELEASE)%.o: %.cpp
	$(MKDIR) $(dir $@)
	$(CC) $(CFLAGS_RELEASE) -c $< -o $@

$(OBJDIR_DEBUG)%.o: %.cpp
	$(MKDIR) $(dir $@)
	$(CC) $(CFLAGS_DEBUG) -c $< -o $@

clean:
	$(RM) -r $(OBJDIR)
	$(RM) $(BINDIR)$(TARGET)
