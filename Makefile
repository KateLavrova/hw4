# По умолчанию линейный режим
all: make8

# Линейный режим
make2: 
	mpic++ main.cpp -o main
	mpiexec -n 2 ./main

# Параллельный режим 
make4:
	mpic++ main.cpp -o main
	mpiexec -n 4 ./main

# Параллельный режим с 2 потоками
make8:
	mpic++ main.cpp -o main
	mpiexec -n 8 ./main

# Очистка
clean:
	rm -f main
