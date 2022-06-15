SHELL := /bin/bash
CPP_FLAG = -lpthread
CPP_FOLDER=cpp/src
CPP_FILE=${CPP_FOLDER}/*.cpp
OBJECT_FILE=${CPP_FOLDER}/*.o


all: lin_reg hawkq nonHawkq 

hawkq:	${OBJECT_FILE} hawk_q.cpp
	g++ hawk_q.cpp $(OBJECT_FILE) -lm -lpthread -o hawk_q
	
kmerTest: ${OBJECT_FILE} kmerTest.cpp
	g++ kmerTest.cpp $(OBJECT_FILE) -lm -lpthread -o kmerTest
	
nonHawkq: 
	g++ bonf_fasta.cpp -o bonf_fasta
	g++ kmersearch.cpp -o kmersearch
	g++ kmersummary.cpp -o kmersummary
	g++ preProcess.cpp -o preProcess
	g++ convertToFasta.cpp -o convertToFasta


lin_reg: ${OBJECT_FILE} lin_reg.cpp
	g++ lin_reg.cpp $(OBJECT_FILE) -lm -lpthread -o $@

log_reg_control.out: log_reg_control.o ${OBJECT_FILE}
	g++ $^ ${CPP_FLAG} -o $@

lin_reg.o: lin_reg.cpp
	g++ $^ -c -o $@

log_reg_control.o: log_reg_control.cpp
	g++ $^  -c -o $@

${OBJECT_FILE}: ${CPP_FILE}
	for f in `ls ${CPP_FILE}`;do echo $$f;g++ $$f -c -o $$f.o; done

clean:
	rm ${CPP_FOLDER}/*.o
	rm *.o
	rm *.out

