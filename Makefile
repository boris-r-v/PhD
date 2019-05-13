all:
	g++ --std=c++11 sum_temp.cpp main.cpp -o calc
	g++ --std=c++11 sum_temp.cpp main2.cpp -o calc2
	g++ --std=c++11 sum_temp.cpp main4.cpp -o calc4

main5:
	g++ --std=c++11 dt_heating_model.cpp main5.cpp -o calc5
	

