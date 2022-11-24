CXX = g++
CXXFLAGS = -std=c++0x -Wall -pedantic-errors -g

SRCS =  main.cpp EqSystemSolver.cpp
OBJS = ${SRCS:.cpp=.o}
HEADERS = EqSystemSolver.hpp

MAIN = main

all: ${MAIN}
        @echo   Simple compilter named myprog has been compiled

${MAIN}: ${OBJS}
        ${CXX} ${CXXFLAGS} ${OBJS} -o ${MAIN}

.cpp.o:
        ${CXX} ${CXXFLAGS} -c $< -o $@

clean:
        ${RM} ${PROGS} ${OBJS} *.o *~.
