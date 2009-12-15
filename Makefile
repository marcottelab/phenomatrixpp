CC=g++
OBJS=test.o oracle.o mindist.o avgmindist.o adjacency_list.o utilities.o hypergeometric.o genephene.o knearest.o marshall.o partialbayes.o distance.o euclidean.o
LIB=/lib/
DRACO_BOOST_FILESYSTEM=/lib/libboost_filesystem-gcc43-mt.a
BOOST_FILESYSTEM=/usr/local/lib/libboost_filesystem.a
DRACO_BOOST_PROGRAMOPTIONS=/lib/libboost_program_options-gcc43-mt.a
BOOST_PROGRAMOPTIONS=/usr/local/lib/libboost_program_options.a
DRACO_BOOST_SYSTEM=/lib/libboost_system-gcc43-mt.a
BOOST_SYSTEM=/usr/local/lib/libboost_system.a
DEBUGOPTS=-g3 -DDEBUG -Wall
OPTS=-O3 -Wall
BIN=phenomatrix
RAILS_DIR=~/NetBeansProjects/crossval/bin/lovelace/

all: ${OBJS}
	${CC} -o ${BIN} ${OBJS} -L${LIB} -l:${BOOST_FILESYSTEM} -l:${BOOST_SYSTEM} -l:${BOOST_PROGRAMOPTIONS} ${OPTS}

test.o: test.cpp mindist.o avgmindist.o partialbayes.o knearest.o ; ${CC} -c test.cpp ${OPTS}

euclidean.o: euclidean.cpp euclidean.h ; ${CC} -c euclidean.cpp ${OPTS}

hypergeometric.o: hypergeometric.cpp hypergeometric.h ; ${CC} -c hypergeometric.cpp ${OPTS}

distance.o: distance.cpp distance.h ; ${CC} -c distance.cpp ${OPTS}

oracle.o: adjacency_list.o oracle.cpp oracle.h genephene.o type_shield.h ; ${CC} -c oracle.cpp ${OPTS}

marshall.o: marshall.cpp marshall.h utilities.o genephene.o distance.o ; ${CC} -c marshall.cpp ${OPTS}

genephene.o: genephene.cpp genephene.h utilities.o adjacency_list.o hypergeometric.o ; ${CC} -c genephene.cpp ${OPTS}

mindist.o: mindist.cpp mindist.h oracle.o ; ${CC} -c mindist.cpp ${OPTS}

avgmindist.o: avgmindist.cpp avgmindist.h oracle.o ; ${CC} -c avgmindist.cpp ${OPTS}

adjacency_list.o: adjacency_list.cpp adjacency_list.h ; ${CC} -c adjacency_list.cpp ${OPTS}

utilities.o: utilities.cpp utilities.h constants.h ; ${CC} -c utilities.cpp ${OPTS}

knearest.o: knearest.cpp knearest.h oracle.o type_shield.h ; ${CC} -c knearest.cpp ${OPTS}

partialbayes.o: partialbayes.cpp partialbayes.h knearest.o ; ${CC} -c partialbayes.cpp ${OPTS}

rails:
	cp ${BIN} ${RAILS_DIR}

clean:
	rm -f *.o ${BIN}
