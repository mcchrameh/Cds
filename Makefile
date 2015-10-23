CC=mpic++   
CFLAGS=-c -Wall  -std=c++11  -g 
LDFLAGS= -lm -pg
SOURCES=testCahnH.cpp CahnH.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=testCahn

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@  $(LDFLAGS)

.cpp.o:
	$(CC)  $(CFLAGS) $< -o $@  $(LDFLAGS)


clean:	 
	rm -f $(OBJECTS)
	rm -f $(EXECUTABLE)

run:
	mpirun -np 2 ./$(EXECUTABLE) 
