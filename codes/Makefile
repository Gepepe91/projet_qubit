# Compilateur et options
CXX = g++
CXXFLAGS = -Wall -std=c++17 -Iinclude

# Répertoires
SRC_DIR = src
INC_DIR = include

# Fichiers sources et objets
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES = $(SRC_FILES:.cpp=.o)

# Nom du programme final
EXEC = programme

# Règle principale
$(EXEC): $(OBJ_FILES) main.o
	$(CXX) $(CXXFLAGS) -o $(EXEC) $^

# Compilation des fichiers sources
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compilation de main séparément
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

# Nettoyage des fichiers objets et fichiers temporaires
clean:
	rm -f $(SRC_DIR)/*.o main.o *~

# Nettoyage complet
mrproper: clean
	rm -f $(EXEC)
