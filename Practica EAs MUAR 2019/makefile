



#
# How to make .c into .o   
#
.c.o: ga.h
	cc  -c $<


ga-test: ga-test.o  
	gcc ga-test.c -o ga-test  -L./libga  -lGA -lm

clean:
	rm -f *~
	rm -f *.o


