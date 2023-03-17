# Configuration pour MSYS 2
# Peut être modifié à votre convenance...
MPICXX = mpic++
#CXXFLAGS = -std=c++17 
CXXFLAGS = -std=c++17 -fPIC
ifdef DEBUG
CXXFLAGS += -g -O0 -Wall -fbounds-check -pedantic -fsanitize=address -fopenmp
CXXFLAGS2 = CXXFLAGS
else
CXXFLAGS2 = ${CXXFLAGS} -O2 -march=native -Wall -fopenmp
CXXFLAGS += -O3 -march=native -Wall -fopenmp
endif
LIB=-lmpi -lpthread -lsfml-graphics -lsfml-window -lsfml-system



