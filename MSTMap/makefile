all: MSTMap.exe

MSTMap.exe: main.o linkage_group_DH.o genetic_map_DH.o  linkage_group_RIL.o genetic_map_RIL.o MSTOpt.o
	g++ -o MSTMap.exe  main.o linkage_group_DH.o genetic_map_DH.o linkage_group_RIL.o genetic_map_RIL.o MSTOpt.o

main.o: main.cpp constants.h
	g++ -o main.o -c -g main.cpp

MSTOpt.o: MSTOpt.cpp MSTOpt.h
	g++ -o MSTOpt.o -c -g MSTOpt.cpp

genetic_map_DH.o: genetic_map_DH.cpp genetic_map_DH.h
	g++ -o genetic_map_DH.o -c -g genetic_map_DH.cpp 

linkage_group_DH.o : linkage_group_DH.h linkage_group_DH.cpp
	g++ -o linkage_group_DH.o  -c -g linkage_group_DH.cpp

genetic_map_RIL.o: genetic_map_RIL.cpp genetic_map_RIL.h
	g++ -o genetic_map_RIL.o -c -g genetic_map_RIL.cpp 

linkage_group_RIL.o : linkage_group_RIL.h linkage_group_RIL.cpp
	g++ -o linkage_group_RIL.o  -c -g linkage_group_RIL.cpp

clean: 
	rm -rf *.o; rm -rf *.exe; rm -rf *~; rm -rf *.out
	
