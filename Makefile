IDIR=.
CC=g++
#LIB=/usr/lib
#CFLAGS=-I${IDIR} -I/usr/local/includes -std=c++11 -L${LIB} -Wall -O3 -DNDEBUG
CFLAGS=-I${IDIR} -std=c++11 -Wall -O3 -DNDEBUG

OBJS = $(patsubst %.cpp,%.o,$(wildcard *.cpp))

all: ${OBJS} generate_dict sparsify threshold  
	echo "All made."

sparsify: ${OBJS} 
	${CC} -o $@ ${CFLAGS} $@.o library.o

generate_dict: generate_dict.o library.o
	${CC} -o $@ ${CFLAGS} $@.o library.o

threshold: threshold.o
	${CC} -o $@ ${CFLAGS} $@.o

%.o: %.cpp global.h
	${CC} ${CFLAGS} -c -o $@ $<

clean:
	rm -f sparsify generate_dict threshold ${OBJS}
	@echo "All cleaned up!"
