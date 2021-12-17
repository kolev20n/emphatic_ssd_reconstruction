LIBS=`root-config --cflags` 
CFLAGS=`root-config --glibs` -lHistPainter
CC=g++
COPT=-g

default:
	$(CC) $(COPT) example.cpp -o run_example $(LIBS) $(CFLAGS)

clean:
	rm run_eample *.o