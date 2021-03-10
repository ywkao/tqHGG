################################################################################## 
# 26. Nov. 2018, Editor: Yu-Wei Kao                                              #
# References:                                                                    # 
# 1) https://hiltmon.com/blog/2013/07/03/a-simple-c-plus-plus-project-structure/ #
# 2) http://mropengate.blogspot.com/2018/01/makefile.html                        # 
# 3) https://stackoverflow.com/questions/5950395/makefile-to-compile-multiple-c-programs
# 4) http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/
# 5) https://www.wooster.edu/_media/files/academics/areas/computer-science/resources/makefile-tut.pdf
##################################################################################
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) -lMinuit
ROOTGLIBS     = $(shell root-config --glibs)

# CC := clang --analyze # and comment out the linker last line for sanity
CC       := g++ # This is the main compiler
SRCDIR   := src
BUILDDIR := build
TARGET   := bin/check
TARGET0  := bin/covarianceMatrixStudy
TARGET1  := bin/preselection
#TARGET1  := bin_2/preselection
#TARGET1  := bin_3/preselection
TARGET2  := bin/selection
TARGET3  := bin/generalChiSquareStudy_hadronic_exe
TARGET3b := bin/generalChiSquareStudy_leptonic_exe
TARGET4  := bin/preselection_npustudy
TARGET5  := bin/myKinFit

SRCEXT   := cpp
SOURCES  := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS  := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS   := $(shell root-config --cflags) -g -O3 #-Wno-write-strings -D_FILE_OFFSET_BITS=64 -DDROP_CGAL #-Wall -Wextra
LIB      := $(shell root-config --libs) -lMinuit
INC      := -I include


all: build/libHistFactory.so ${TARGET} ${TARGET0} ${TARGET1} ${TARGET2} $(TARGET3) $(TARGET3b) $(TARGET4)
#all: ${TARGET} ${TARGET0} ${TARGET1} ${TARGET2} $(TARGET3) $(TARGET4) $(TARGET5)


build/libHistFactory.so: include/hist_factory.h include/hist_factory.cpp
	      @echo "####### Building library build/libHistFactory.so"
		  @gcc -fPIC -shared $(ROOTCFLAGS) $(ROOTLIBS) -I. include/hist_factory.cpp -o build/libHistFactory.so


#$(TARGET): $(OBJECTS)
#	@echo " Linking..."
#	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)
$(TARGET): build/check.o
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(TARGET0): build/covarianceMatrixStudy.o
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET0) $(LIB)"; $(CC) $^ -o $(TARGET0) $(LIB)

$(TARGET1): build/preselection.o build/genMatching.o
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET1) $(LIB)"; $(CC) $^ -o $(TARGET1) $(LIB)

$(TARGET2): build/selection.o build/ctag_reshaping.o
	@echo " Linking for selection cpp..."
	@echo " $(CC) $^ -o $(TARGET2) -L/wk_cms2/ykao/CMSSW_9_4_10/src/2017/TopKinFit -lKinFit $(LIB)"; $(CC) $^ -o $(TARGET2) -L/wk_cms2/ykao/CMSSW_9_4_10/src/2017/TopKinFit -lKinFit $(LIB)

$(TARGET3): build/generalChiSquareStudy_hadronic_exe.o build/libHistFactory.so build/gen_reco_performance_helper.o
	@echo " Linking for generalChiSquareStudy_hadronic_exe cpp..."
	@echo " $(CC) $^ -o $(TARGET3) $(LIB)"; $(CC) $^ -o $(TARGET3) $(LIB)

$(TARGET3b): build/generalChiSquareStudy_leptonic_exe.o build/libHistFactory.so
	@echo " Linking for generalChiSquareStudy_leptonic_exe cpp..."
	@echo " $(CC) $^ -o $(TARGET3b) -L/wk_cms2/ykao/CMSSW_9_4_10/src/2017/TopKinFit -lKinFit $(LIB)"; $(CC) $^ -o $(TARGET3b) -L/wk_cms2/ykao/CMSSW_9_4_10/src/2017/TopKinFit -lKinFit $(LIB)

$(TARGET4): build/preselection_npustudy.o
	@echo " Linking for preselection_npustudy cpp..."
	@echo " $(CC) $^ -o $(TARGET4) $(LIB)"; $(CC) $^ -o $(TARGET4) $(LIB)

$(TARGET5): build/myKinFit.o
	@echo " Linking for myKinFit cpp..."
	@echo " $(CC) $^ -o $(TARGET5) -L/wk_cms2/ykao/CMSSW_9_4_10/src/2017/TopKinFit -lKinFit $(LIB)"; $(CC) $^ -o $(TARGET5) -L/wk_cms2/ykao/CMSSW_9_4_10/src/2017/TopKinFit -lKinFit $(LIB) 



#$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
build/check.o: src/check.cpp
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

build/covarianceMatrixStudy.o: src/covarianceMatrixStudy.cpp
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

build/preselection.o: src/preselection_exe.cpp
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

build/genMatching.o: include/genMatching.cc
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

build/ctag_reshaping.o: include/ctag_reshaping.cpp
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

build/gen_reco_performance_helper.o: include/gen_reco_performance_helper.cc include/gen_reco_performance_helper.h
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

build/selection.o: src/selection.cpp
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) -I /wk_cms2/ykao/CMSSW_9_4_10/src/2017/TopKinFit $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) -I /wk_cms2/ykao/CMSSW_9_4_10/src/2017/TopKinFit $(CFLAGS) $(INC) -c -o $@ $<

build/generalChiSquareStudy_hadronic_exe.o: src/generalChiSquareStudy_hadronic_exe.cpp
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

build/generalChiSquareStudy_leptonic_exe.o: src/generalChiSquareStudy_leptonic_exe.cpp
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) -I /wk_cms2/ykao/CMSSW_9_4_10/src/2017/TopKinFit $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) -I /wk_cms2/ykao/CMSSW_9_4_10/src/2017/TopKinFit $(CFLAGS) $(INC) -c -o $@ $<

build/preselection_npustudy.o: src/preselection_npustudy.cpp
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

build/myKinFit.o: src/myKinFit.cpp
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) -I /wk_cms2/ykao/CMSSW_9_4_10/src/2017/TopKinFit $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) -I /wk_cms2/ykao/CMSSW_9_4_10/src/2017/TopKinFit $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET) $(TARGET0) $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET3b) $(TARGET4) $(TARGET5)"
	$(RM) -r $(BUILDDIR) $(TARGET) $(TARGET0) $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET3b) $(TARGET4) $(TARGET5)



# Tests
tester:
	$(CC) $(CFLAGS) test/tester.cpp $(INC) $(LIB) -o bin/tester

# Spikes
ticket:
	$(CC) $(CFLAGS) spikes/ticket.cpp $(INC) $(LIB) -o bin/ticket

.PHONY: clean
