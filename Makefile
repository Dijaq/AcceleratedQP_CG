all:
	g++ -I /usr/local/include main.cpp -o main -std=c++11 `pkg-config --cflags --libs opencv`
exec:
	g++ -I /usr/local/include main2.cpp -o main2 -std=c++11