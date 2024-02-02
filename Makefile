#
# Note to myself: When using the GNU build tools instead of this simple
# Makefile, I must do these steps:
# 
#	1) Add #include "config.h" to the C source file.
#	2) Create the files "configure.ac", "Makefile.am".
#	3) Run the program "autoreconf --install".
#	4) ./configure
#	5) make
#	6) make distcheck
#
#
all: starchy_create starchy_search

CC = gcc
INCLUDE = .
CFLAGS = -Wall
PREFIX = /usr/local
CLIBS = -lm

starchy_search: libstarchy.a starchy_search.c
	$(CC) starchy_search.c libstarchy.a -I$(INCLUDE) $(CFLAGS) $(CLIBS) -o starchy_search

starchy_create: libstarchy.a starchy_create.c
	$(CC) starchy_create.c libstarchy.a -I$(INCLUDE) $(CFLAGS) $(CLIBS) -o starchy_create

libstarchy.a: starchy.o
	ar crv libstarchy.a starchy.o
	ranlib libstarchy.a

starchy.o: starchy.c starchy.h
	$(CC) -I$(INCLUDE) $(CFLAGS) $(CLIBS) -c starchy.c

install: libstarchy.a
	@if [ -d $(PREFIX) ]; then \
	   cp libstarchy.a $(PREFIX)/lib; \
	   echo "Installed in $(PREFIX)/lib"; \
	   cp starchy.h $(PREFIX)/include; \
	   cp starchy_index.member $(PREFIX)/share; \
	   cp starchy_index.star $(PREFIX)/share; \
	   cp starchy_index.ttree $(PREFIX)/share; \
	   cp starchy_index.utree $(PREFIX)/share; \
	else \
	   echo "Sorry, $(PREFIX) does not exist"; \
	fi

uninstall:
	-rm $(PREFIX)/lib/libstarchy.a
	-rm $(PREFIX)/include/starchy.h
	-rm $(PREFIX)/share/starchy_index.member
	-rm $(PREFIX)/share/starchy_index.star
	-rm $(PREFIX)/share/starchy_index.ttree
	-rm $(PREFIX)/share/starchy_index.utree

clean:
	-rm starchy.o
