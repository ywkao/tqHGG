################################################################################## 
# 26. Nov. 2018, Editor: Yu-Wei Kao                                              #
# References:                                                                    # 
# 1) https://hiltmon.com/blog/2013/07/03/a-simple-c-plus-plus-project-structure/ #
# 2) http://mropengate.blogspot.com/2018/01/makefile.html                        # 
# 3) https://stackoverflow.com/questions/5950395/makefile-to-compile-multiple-c-programs
# 4) http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/
# 5) https://www.wooster.edu/_media/files/academics/areas/computer-science/resources/makefile-tut.pdf
##################################################################################

# CC := clang --analyze # and comment out the linker last line for sanity
CC       := g++ # This is the main compiler
SRCDIR   := src
BUILDDIR := build
TARGET1  := bin/preselection
TARGET2  := bin/selection
TARGET3  := bin/pustudy

SRCEXT   := cpp
SOURCES  := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS  := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS   := $(shell root-config --cflags) -g -O3 #-Wno-write-strings -D_FILE_OFFSET_BITS=64 -DDROP_CGAL #-Wall -Wextra
LIB      := $(shell root-config --libs) -lMinuit
INC      := -I include

all: ${TARGET1} ${TARGET2} ${TARGET3}

#$(TARGET): $(OBJECTS)
#	@echo " Linking..."
#	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)
$(TARGET1): build/preselection.o
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET1) $(LIB)"; $(CC) $^ -o $(TARGET1) $(LIB)

$(TARGET2): build/selection.o
	@echo " Linking for selection cpp..."
	@echo " $(CC) $^ -o $(TARGET2) $(LIB)"; $(CC) $^ -o $(TARGET2) $(LIB)

$(TARGET3): build/study_pileup_reweighting.o
	@echo " Linking for selection cpp..."
	@echo " $(CC) $^ -o $(TARGET3) $(LIB)"; $(CC) $^ -o $(TARGET3) $(LIB)

#$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
build/preselection.o: src/preselection.cpp
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

build/selection.o: src/selection.cpp
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

build/study_pileup_reweighting.o: src/study_pileup_reweighting.cpp
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET1) $(TARGET2) $(TARGET3)"; $(RM) -r $(BUILDDIR) $(TARGET1) $(TARGET2) $(TARGET3)



# Tests
tester:
	$(CC) $(CFLAGS) test/tester.cpp $(INC) $(LIB) -o bin/tester

# Spikes
ticket:
	$(CC) $(CFLAGS) spikes/ticket.cpp $(INC) $(LIB) -o bin/ticket

.PHONY: clean
