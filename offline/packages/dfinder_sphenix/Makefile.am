AUTOMAKE_OPTIONS = foreign subdir-objects

lib_LTLIBRARIES = libdfinder_sphenix.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs`

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -DHomogeneousField

pkginclude_HEADERS = \
  DFinderCombinatorics.h \
  DFinderConfig.h \
  DFinderFitter.h \
  DFinderFun4All.h \
  DFinderNtuple.h \
  DFinderSphenix.h

libdfinder_sphenix_la_SOURCES = \
  DFinderCombinatorics.cc \
  DFinderFitter.cc \
  DFinderNtuple.cc \
  DFinderSphenix.cc

libdfinder_sphenix_la_LIBADD = \
  -lcentrality \
  -lffaobjects \
  -lfun4all \
  -lglobalvertex \
  -lKFParticle \
  -lkfparticle_sphenix \
  -lphool \
  -lSubsysReco \
  -ltrackbase_historic \
  -ltrack_reco

check_PROGRAMS = test_combinatorics
TESTS = test_combinatorics

test_combinatorics_SOURCES = tests/test_combinatorics.cc
test_combinatorics_LDADD = libdfinder_sphenix.la
