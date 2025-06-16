all: ./exec/main# ./exec/autodiffusion ./exec/corr_func

OBJS = ./exec/verlet_periodic.o ./exec/in_cond.o ./exec/write_to_vtk.o

#CC = cc -O3 -ggdb
CC = cc -O3 -ggdb -Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp
#CC = /opt/homebrew/bin/gcc-15 -O3 -ggdb -fopenmp

./exec/main: ./exec/main.o $(OBJS)
	$(CC) -o $@ ./exec/main.o $(OBJS) -lm

./exec/autodiffusion: ./exec/autodiffusion.o $(OBJS)
	$(CC) -o $@ ./exec/autodiffusion.o $(OBJS) -lm

./exec/corr_func: ./exec/corr_func.o $(OBJS)
	$(CC) -o $@ ./exec/corr_func.o $(OBJS) -lm

./exec/main.o: main.c
	$(CC) -o $@ -c $<

./exec/autodiffusion.o: autodiffusion.c
	$(CC) -o $@ -c $<

./exec/corr_func.o: corr_func.c
	$(CC) -o $@ -c $<

./exec/verlet_periodic.o: verlet_periodic.c
	$(CC) -o $@ -c $<

./exec/in_cond.o: in_cond.c
	$(CC) -o $@ -c $<

./exec/write_to_vtk.o: write_to_vtk.c
	$(CC) -o $@ -c $<

clean:
	rm -f ./exec/main ./exec/autodiffusion ./*.o ./exec/*.o ./exec/corr_func
