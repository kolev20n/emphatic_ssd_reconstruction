LIBS=`root-config --evelibs` 
CFLAGS=`root-config --cflags`
CC=g++
COPT=-g

default:
	$(CC) $(COPT) run_vis.cpp -o run_vis $(LIBS) $(CFLAGS)

clean:
	rm run_vis *.o