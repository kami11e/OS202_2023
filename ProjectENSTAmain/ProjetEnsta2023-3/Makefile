# include à modifier selon sa configuration
# include make_msys.inc 
include make_linux.inc


ALL= vortexSimulation.exe

ALL2= vortexSimulationMPI.exe

ALL3= vortexSimulationMPI2.exe

ALL4 = vortexSimulationMPI3.exe


default:	help
all: $(ALL)
mpi: $(ALL2)
mpi2: $(ALL3)
mpi3: $(ALL4)

clean:
	@rm -fr objs/*.o *.exe src/*~ *.png

OBJS1= objs/vortex.o objs/screen.o objs/runge_kutta.o objs/cloud_of_points.o objs/cartesian_grid_of_speed.o \
      objs/vortexSimulation.o

OBJS2= objs/vortex.o objs/screen.o objs/runge_kutta.o objs/cloud_of_points.o objs/cartesian_grid_of_speed.o \
      objs/vortexSimulationMPI.o

OBJS3= objs/vortex.o objs/screen.o objs/runge_kutta.o objs/cloud_of_points.o objs/cartesian_grid_of_speed.o \
      objs/vortexSimulationMPI2.o

OBJS4= objs/vortex.o objs/screen.o objs/runge_kutta.o objs/cloud_of_points.o objs/cartesian_grid_of_speed.o \
      objs/vortexSimulationMPI3.o

objs/vortex.o:	src/point.hpp src/vector.hpp src/vortex.hpp src/vortex.cpp
	$(CXX) $(CXXFLAGS) -Isrc -I"/mnt/c/Libraries/SFML-2.5.1/include" -c -o $@ src/vortex.cpp

objs/cartesian_grid_of_speed.o: src/point.hpp src/vector.hpp src/vortex.hpp src/cartesian_grid_of_speed.hpp src/cartesian_grid_of_speed.cpp
	$(CXX) $(CXXFLAGS) -Isrc -I"/mnt/c/Libraries/SFML-2.5.1/include" -c -o $@ src/cartesian_grid_of_speed.cpp

objs/cloud_of_points.o: src/point.hpp src/rectangle.hpp src/cloud_of_points.hpp src/cloud_of_points.cpp 
	$(CXX) $(CXXFLAGS) -Isrc -I"/mnt/c/Libraries/SFML-2.5.1/include" -c -o $@ src/cloud_of_points.cpp 

objs/runge_kutta.o:	src/vortex.hpp src/cloud_of_points.hpp src/cartesian_grid_of_speed.hpp src/runge_kutta.hpp src/runge_kutta.cpp 
	$(MPICXX) $(CXXFLAGS) -Isrc -I"/mnt/c/Libraries/SFML-2.5.1/include" -c -o $@ src/runge_kutta.cpp

objs/screen.o:	src/vortex.hpp src/cloud_of_points.hpp src/cartesian_grid_of_speed.hpp src/screen.hpp src/screen.cpp
	$(CXX) $(CXXFLAGS) -Isrc -I"/mnt/c/Libraries/SFML-2.5.1/include" -c -o $@ src/screen.cpp

objs/vortexSimulation.o: src/cartesian_grid_of_speed.hpp src/vortex.hpp src/cloud_of_points.hpp src/runge_kutta.hpp src/screen.hpp src/vortexSimulation.cpp
	$(CXX) $(CXXFLAGS) -Isrc -I"/mnt/c/Libraries/SFML-2.5.1/include" -c -o $@ src/vortexSimulation.cpp

vortexSimulation.exe: $(OBJS1)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS1) -L /mnt/c/Libraries/SFML-2.5.1/lib $(LIB)

objs/vortexSimulationMPI.o: src/cartesian_grid_of_speed.hpp src/vortex.hpp src/cloud_of_points.hpp src/runge_kutta.hpp src/screen.hpp src/vortexSimulationMPI.cpp
	$(MPICXX) $(CXXFLAGS) -Isrc -I"/mnt/c/Libraries/SFML-2.5.1/include" -c -o $@ src/vortexSimulationMPI.cpp

vortexSimulationMPI.exe: $(OBJS2)
	$(MPICXX) $(CXXFLAGS) -o $@ $(OBJS2) -L /mnt/c/Libraries/SFML-2.5.1/lib $(LIB)

objs/vortexSimulationMPI2.o: src/cartesian_grid_of_speed.hpp src/vortex.hpp src/cloud_of_points.hpp src/runge_kutta.hpp src/screen.hpp src/vortexSimulationMPI2.cpp
	$(MPICXX) $(CXXFLAGS) -Isrc -I"/mnt/c/Libraries/SFML-2.5.1/include" -c -o $@ src/vortexSimulationMPI2.cpp

vortexSimulationMPI2.exe: $(OBJS3)
	$(MPICXX) $(CXXFLAGS) -o $@ $(OBJS3) -L /mnt/c/Libraries/SFML-2.5.1/lib $(LIB)

objs/vortexSimulationMPI3.o: src/cartesian_grid_of_speed.hpp src/vortex.hpp src/cloud_of_points.hpp src/runge_kutta.hpp src/screen.hpp src/vortexSimulationMPI3.cpp
	$(MPICXX) $(CXXFLAGS) -Isrc -I"/mnt/c/Libraries/SFML-2.5.1/include" -c -o $@ src/vortexSimulationMPI3.cpp

vortexSimulationMPI3.exe: $(OBJS4)
	$(MPICXX) $(CXXFLAGS) -o $@ $(OBJS4) -L /mnt/c/Libraries/SFML-2.5.1/lib $(LIB)

help:
	@echo "Available targets : "
	@echo "    all                           : compile all executables"
	@echo "    vortexSimulation.exe          : compile simple this executable"
	@echo "Add DEBUG=yes to compile in debug"
	@echo "Configuration :"
	@echo "    CXX      :    $(CXX)"
	@echo "    CXXFLAGS :    $(CXXFLAGS)"
