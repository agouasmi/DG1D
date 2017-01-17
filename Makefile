EXEC = DG1D
OBJECTS = main.o  
SOURCES = main.cpp Element.h DG_1D.h Problem_functions.h  
CFLAGS = -c -g -std=c++11 -I/usr/include -I/usr/dislin
LFLAGS = -L/usr/local/dislin -ldiscpp -ldislin -lXt -lm 
CC = g++


${EXEC}: ${OBJECTS}
	${CC} ${OBJECTS} ${LFLAGS} -o DG1D

${OBJECTS}: ${SOURCES}
	${CC} ${CFLAGS} ${SOURCES}

clean:
	rm -r ${OBJECTS} ${EXEC} 
	
run:
	./burgers Run/input.txt

