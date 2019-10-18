#
# $Id: Makefile,v 1.3 2001/12/06 17:27:45 botzler Exp $
#

CC=gcc
#CFLAGS= -O3 -g
CFLAGS= -O2
LINK=gcc
LIBS=-lm

TARGET=varea

all: $(TARGET)

$(TARGET): varea.o triangle.o
	$(LINK) varea.o triangle.o -o $(TARGET) $(LIBS)

varea.o: varea.c triangle.h
	$(CC) $(CFLAGS) -c -o varea.o varea.c

triangle.o: triangle.c triangle.h
	$(CC) $(CFLAGS) -DTRILIBRARY -c -o triangle.o triangle.c

clean: 
	rm -f *~ *.o $(TARGET) core

install:
	cp $(TARGET) $(HOME)/bin/
