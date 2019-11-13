#SUBDIRS = test
dist_doc_DATA = README.txt

AM_CFLAGS    = -I./src -Iinclude -I$(srcdir)/include
lib_LTLIBRARIES = libmusic.la
libmusic_la_SOURCES = src/lm.c
include_HEADERS = include/lm.h
libmusic_la_LDFLAGS = -version-info 1:0:0

bin_PROGRAMS = lmtool
lmtool_SOURCES = src/lmtool.c
lmtool_LDADD = libmusic.la

check_PROGRAMS = lmtest1 lmtest2 lmtest3 lmtest4 lmtest5 lmtest6 lmtest7 lmtest8 lmtest9 lmtest10 lmtest11
TESTS = $(check_PROGRAMS)

# MUSIC algorithm tests
LIBS += -llapack -lm

# Test internals

lmtest1_SOURCES = test/lmtest1.c src/lm.c lm.h
lmtest1_CFLAGS = -Iinclude

# Test SVD

lmtest2_SOURCES = test/lmtest2.c src/lm.c lm.h
lmtest2_CFLAGS = -Iinclude

lmtest3_SOURCES = test/lmtest3.c src/lm.c lm.h
lmtest3_CFLAGS = -Iinclude

lmtest4_SOURCES = test/lmtest4.c src/lm.c lm.h
lmtest4_CFLAGS = -Iinclude

# Test correlation

lmtest5_SOURCES = test/lmtest5.c src/lm.c lm.h
lmtest5_CFLAGS = -Iinclude

# Test MUSIC algorithm, lm_detect

lmtest6_SOURCES = test/lmtest6.c src/lm.c lm.h
lmtest6_CFLAGS = -Iinclude

lmtest7_SOURCES = test/lmtest7.c src/lm.c lm.h
lmtest7_CFLAGS = -Iinclude

lmtest8_SOURCES = test/lmtest8.c src/lm.c lm.h
lmtest8_CFLAGS = -Iinclude

# Test DTMF detect with MUSIC algorithm, lm_standard_init and lm_dtmf_init

lmtest9_SOURCES = test/lmtest9.c src/lm.c lm.h
lmtest9_CFLAGS = -Iinclude

# Test DTMF detect with MUSIC algorithm, lm_dtmf_detect, on reference data

lmtest10_SOURCES = test/lmtest10.c src/lm.c lm.h
lmtest10_CFLAGS = -Iinclude

# Test DTMF detect with MUSIC algorithm, lm_dtmf_detect, external file

lmtest11_SOURCES = test/lmtest11.c src/lm.c lm.h
lmtest11_CFLAGS = -Iinclude
