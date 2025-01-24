programme: qubit.cpp matrice.cpp rk4_qubit.cpp
	g++ -Wall -std=c++17 -o programme qubit.cpp matrice.cpp rk4_qubit.cpp \
	-I/opt/anaconda3/envs/Projet_C++/include/python3.13 \
	-L/opt/anaconda3/envs/Projet_C++/lib -lpython3.13
	install_name_tool -add_rpath /opt/anaconda3/envs/Projet_C++/lib ./programme

clean:
	rm -f *~

mrproper: clean
	rm -f programme