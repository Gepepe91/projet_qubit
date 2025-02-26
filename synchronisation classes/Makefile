programme: qubit.cpp tens_2q.cpp matrice.cpp main.cpp
	g++ -Wall -std=c++17 -o programme qubit.cpp tens_2q.cpp matrice.cpp main.cpp \
	-I/opt/anaconda3/envs/Projet_C++/include/python3.13 \
	-L/opt/anaconda3/envs/Projet_C++/lib -lpython3.13
	install_name_tool -add_rpath /opt/anaconda3/envs/Projet_C++/lib ./programme

clean:
	rm -f *~

mrproper: clean
	rm -f programme