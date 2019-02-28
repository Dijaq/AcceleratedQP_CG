all:
	g++ -I /usr/local/include main.cpp -o main -std=c++11 `pkg-config --cflags --libs opencv`
exec:
	g++ -I /usr/local/include practice_v2.cpp -o practice -std=c++11