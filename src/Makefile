Cells : Cells.o mersenne.o Queue.o cartesian.o Simulation.o
	gcc -g -o ../Cells Cells.o mersenne.o Queue.o Simulation.o cartesian.o -lm

Cells.o : Cells.c
	gcc -Wall -c -g Cells.c -lm

Simulation.o : Simulation.c
	gcc -Wall -c -g Simulation.c -lm

Queue.o : Queue.c
	gcc -Wall -c -g Queue.c -lm

cartesian.o : cartesian.c
	gcc -Wall -c -g cartesian.c -lm

mersenne.o : mersenne.c 
	gcc -Wall -c -g mersenne.c -lm

test : Cells.o mersenne.o Queue.o cartesian.o Cells-test.o
	gcc -g -o test Cells.o mersenne.o Queue.o cartesian.o Cells-test.o -lm

Cells-test.o : Cells-test.c 
	gcc -Wall -c -g Cells-test.c -lm

NetDetails : Cells.o mersenne.o Queue.o cartesian.o NetDetails.o
	gcc -g -o NetDetails Cells.o mersenne.o Queue.o cartesian.o NetDetails.o -lm

NetDetails.o : NetDetails.c 
	gcc -Wall -c -g NetDetails.c -lm


clean :
	rm -f *.o test ../Cells
