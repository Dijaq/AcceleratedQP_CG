all:
	nvcc -I /usr/local/include main.cu -o mcuda -std=c++11 `pkg-config --cflags --libs opencv`
exec:
	g++ -I /usr/local/include main.cpp -o main -std=c++11 -lcudart `pkg-config --cflags --libs opencv`
	nvcc -I /usr/local/include practice.cu -o practice -std=c++11 -lglut -lGL -lGLEW `pkg-config --cflags --libs opencv` 	 
	g++ -I /usr/local/include practice.cpp -o practice -std=c++11 -lglut -lglfw -lGL -lGLEW `pkg-config --cflags --libs opencv` 	
	g++ -I /usr/local/include cow_IsoDist.cpp -o main2 -std=c++11