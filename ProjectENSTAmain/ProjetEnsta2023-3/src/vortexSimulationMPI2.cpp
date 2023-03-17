#include <SFML/Window/Keyboard.hpp>
#include <ios>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <chrono>
#include "cartesian_grid_of_speed.hpp"
#include "vortex.hpp"
#include "cloud_of_points.hpp"
#include "runge_kutta.hpp"
#include "screen.hpp"
extern "C"{
    #include <pthread.h>
}
#include <mpi.h>

auto readConfigFile( std::ifstream& input )
{
    using point=Simulation::Vortices::point;

    int isMobile;
    std::size_t nbVortices;
    Numeric::CartesianGridOfSpeed cartesianGrid;
    Geometry::CloudOfPoints cloudOfPoints;
    constexpr std::size_t maxBuffer = 8192;
    char buffer[maxBuffer];
    std::string sbuffer;
    std::stringstream ibuffer;
    // Lit la première ligne de commentaire :
    input.getline(buffer, maxBuffer); // Relit un commentaire
    input.getline(buffer, maxBuffer);// Lecture de la grille cartésienne
    sbuffer = std::string(buffer,maxBuffer);
    ibuffer = std::stringstream(sbuffer);
    double xleft, ybot, h;
    std::size_t nx, ny;
    ibuffer >> xleft >> ybot >> nx >> ny >> h;
    cartesianGrid = Numeric::CartesianGridOfSpeed({nx,ny}, point{xleft,ybot}, h);
    input.getline(buffer, maxBuffer); // Relit un commentaire
    input.getline(buffer, maxBuffer); // Lit mode de génération des particules
    sbuffer = std::string(buffer,maxBuffer);
    ibuffer = std::stringstream(sbuffer);
    int modeGeneration;
    ibuffer >> modeGeneration;
    if (modeGeneration == 0) // Génération sur toute la grille 
    {
        std::size_t nbPoints;
        ibuffer >> nbPoints;
        cloudOfPoints = Geometry::generatePointsIn(nbPoints, {cartesianGrid.getLeftBottomVertex(), cartesianGrid.getRightTopVertex()});
    }
    else 
    {
        std::size_t nbPoints;
        double xl, xr, yb, yt;
        ibuffer >> xl >> yb >> xr >> yt >> nbPoints;
        cloudOfPoints = Geometry::generatePointsIn(nbPoints, {point{xl,yb}, point{xr,yt}});
    }
    // Lit le nombre de vortex :
    input.getline(buffer, maxBuffer); // Relit un commentaire
    input.getline(buffer, maxBuffer); // Lit le nombre de vortex
    sbuffer = std::string(buffer, maxBuffer);
    ibuffer = std::stringstream(sbuffer);
    try {
        ibuffer >> nbVortices;        
    } catch(std::ios_base::failure& err)
    {
        std::cout << "Error " << err.what() << " found" << std::endl;
        std::cout << "Read line : " << sbuffer << std::endl;
        throw err;
    }
    Simulation::Vortices vortices(nbVortices, {cartesianGrid.getLeftBottomVertex(),
                                               cartesianGrid.getRightTopVertex()});
    input.getline(buffer, maxBuffer);// Relit un commentaire
    for (std::size_t iVortex=0; iVortex<nbVortices; ++iVortex)
    {
        input.getline(buffer, maxBuffer);
        double x,y,force;
        std::string sbuffer(buffer, maxBuffer);
        std::stringstream ibuffer(sbuffer);
        ibuffer >> x >> y >> force;
        vortices.setVortex(iVortex, point{x,y}, force);
    }
    input.getline(buffer, maxBuffer);// Relit un commentaire
    input.getline(buffer, maxBuffer);// Lit le mode de déplacement des vortex
    sbuffer = std::string(buffer,maxBuffer);
    ibuffer = std::stringstream(sbuffer);
    ibuffer >> isMobile;
    return std::make_tuple(vortices, isMobile, cartesianGrid, cloudOfPoints);
}




int main( int nargs, char* argv[] )
{
    #define P_AFFICHAGE 0
    #define P_CALCUL_1 1
    // MPI_Comm global, scomm;

    MPI_Comm global;
    int rank, nbp, provided;
    MPI_Init_thread(&nargs, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_dup(MPI_COMM_WORLD, &global);
    MPI_Comm_size(global, &nbp);
    MPI_Comm_rank(global, &rank);
    // int color = int(rank!=P_AFFICHAGE);     // affichage: color=0, calcul:color=1;
    // MPI_Comm_split(MPI_COMM_WORLD, color, rank, &scomm);

    int N_calcul = nbp-1;

    if (nargs==1)
    {
        std::cout << "Usage : vortexsimulator <nom fichier configuration>" << std::endl;
        return EXIT_FAILURE;
    }

    char const* filename;
    filename = argv[1];
    std::ifstream fich(filename);
    auto config = readConfigFile(fich);
    fich.close();
    Simulation::Vortices vortices = std::get<0>(config);
    int isMobile = std::get<1>(config);
    Numeric::CartesianGridOfSpeed grid     = std::get<2>(config);
    Geometry::CloudOfPoints cloud  = std::get<3>(config);
    grid.updateVelocityField(vortices);

    std::cout<<cloud.numberOfPoints()<<std::endl;



    double dt = 0.1;
    unsigned long nbOfVts=vortices.numberOfVortices();
    unsigned long nbOfPts=cloud.numberOfPoints();
    unsigned long subnbOfPts0 = nbOfPts/N_calcul;
    // unsigned long subnbOfPts = subnbOfPts0;
    // if(rank==nbp-1){
    //     subnbOfPts = nbOfPts - (N_calcul-1)*subnbOfPts0;
    // }
    int* recvcountsCloud = new int[nbp];    recvcountsCloud[0]=0;
    int* displsCloud = new int[nbp];    displsCloud[0]=0;
    for(int i=1;i<nbp;i++){
        if(i!=nbp-1){
            recvcountsCloud[i] = subnbOfPts0*2;
        }else{
            recvcountsCloud[i] = (nbOfPts - (N_calcul-1)*subnbOfPts0)*2;
        }
        
        displsCloud[i] = displsCloud[i-1]+recvcountsCloud[i-1];
        // std::cout<<"i = "<<i<<", displ = "<<displsCloud[i]<<", recvcount = "<<recvcountsCloud[i]<<std::endl;
    }
    // recvcountsCloud[N_calcul-1] = subnbOfPts*2;


    std::pair<std::size_t,std::size_t> gridSize = grid.cellGeometry();
    std::size_t gridHeight = gridSize.second;
    // std::size_t subHeight = gridHeight/N_calcul;
    std::size_t nbOfGridData = gridHeight * gridSize.first;
    // std::size_t subnbOfGridData = nbOfGridData/N_calcul;

    if(rank != P_AFFICHAGE){

        int signal;
        // MPI_Status status;
        MPI_Bcast(&signal, 1, MPI_INT, P_AFFICHAGE, MPI_COMM_WORLD);
        std::cout<<"Process["<<rank<<"] recv sig="<<signal<<" from Process["<<P_AFFICHAGE<<"]"<<std::endl;
        
        while(signal>=0){
            if(signal>=2){
                if(signal==2){
                    dt *=2;
                }else{
                    dt /=2;
                }
            }else{
                std::vector<unsigned long> sendInfo;
                Geometry::CloudOfPoints subCloud;
                std::vector<double> newCloudData(nbOfPts*2);
                if(signal){
                    // std::cout<<"rank = "<<rank<<", displ = "<<displsCloud[rank-1]<<", recvcount = "<<recvcountsCloud[rank-1]<<std::endl;

                    subCloud = sub_solve_RK4_vortices_new_cloud(dt, grid, cloud, rank, subnbOfPts0, recvcountsCloud[rank]);
                    // std::cout<<rank<<":"<<subCloud.numberOfPoints()<<":"<<recvcountsCloud[rank-1]<<":"<<displsCloud[rank-1]<<std::endl;
                    // MPI_Allgather(subCloud.data(), subnbOfPts*2, MPI_DOUBLE, newCloudData.data(), subnbOfPts*2,  MPI_DOUBLE, scomm);
                    // std::cout<<rank<<": Before Allgatherv()"<<std::endl;
                    MPI_Allgatherv(subCloud.data(), recvcountsCloud[rank], MPI_DOUBLE, newCloudData.data(), recvcountsCloud, displsCloud,  MPI_DOUBLE, global);
                    cloud = Geometry::CloudOfPoints (newCloudData);
                    // if(rank==1){
                    //     MPI_Send(cloud.data(), nbOfPts* 2, MPI_DOUBLE, P_AFFICHAGE, 400, MPI_COMM_WORLD);
                    // }
                    solve_RK4_movable_vortices_update_vortices(dt, grid, vortices, cloud);
                    if(rank==1){
                        MPI_Send(vortices.data(), nbOfVts* 3, MPI_DOUBLE, P_AFFICHAGE, 401, MPI_COMM_WORLD);
                        // std::cout<<"rank="<<rank<<", subPts="<<subnbOfPts<<", nbPts="<<nbOfPts<<", nbVts="<<nbOfVts<<std::endl;
                    }
                    // grid.subUpdateVelocitySubField(vortices, rank, subHeight);
                    // MPI_Allgather(grid.data()+subnbOfGridData*2*(rank-1), subnbOfGridData*2, MPI_DOUBLE, grid.data(), subnbOfGridData0*2,  MPI_DOUBLE, scomm);
                    grid.updateVelocityField(vortices);
                    if(rank==1){
                        MPI_Send(grid.data(), nbOfGridData*2, MPI_DOUBLE, P_AFFICHAGE, 402, MPI_COMM_WORLD);
                        // std::cout<<"rank="<<rank<<",send nbOfGridVal="<<nbOfGridData<<std::endl;
                    }
                    
                    
                }else{
                    // std::cout<<"rank="<<rank<<", subPts="<<subnbOfPts<<", nbPts="<<nbOfPts<<std::endl;
                    subCloud = sub_solve_RK4_vortices_new_cloud(dt, grid, cloud, rank, subnbOfPts0, recvcountsCloud[rank]);
                    // MPI_Allgather(subCloud.data(), subnbOfPts*2, MPI_DOUBLE, newCloudData.data(), subnbOfPts*2,  MPI_DOUBLE, scomm);
                    // MPI_Allgatherv(subCloud.data(), subnbOfPts*2, MPI_DOUBLE, newCloudData.data(), recvcountsCloud, displsCloud,  MPI_DOUBLE, scomm);
                    // std::cout<<rank<<": Before Allgatherv()"<<std::endl;
// 
                    MPI_Allgatherv(subCloud.data(), recvcountsCloud[rank], MPI_DOUBLE, newCloudData.data(), recvcountsCloud, displsCloud,  MPI_DOUBLE, global);
                    // std::cout<<rank<<": After Allgatherv()"<<std::endl;

                    cloud = Geometry::CloudOfPoints (newCloudData);
                    // if(rank==1){
                    //     std::cout<<rank<<": Before Send() cloud"<<std::endl;
                    //     MPI_Send(cloud.data(), nbOfPts* 2, MPI_DOUBLE, P_AFFICHAGE, 400, MPI_COMM_WORLD);
                    //     std::cout<<rank<<": after Send() cloud"<<std::endl;
                    // }
                    // MPI_Send(cloud.data(), nbOfPts* 2, MPI_DOUBLE, P_AFFICHAGE, 400, MPI_COMM_WORLD);
                }
 

            }
            // MPI_Recv(&signal,  1, MPI_INT, P_AFFICHAGE, 200, MPI_COMM_WORLD, &status);
            MPI_Bcast(&signal, 1, MPI_INT, P_AFFICHAGE, MPI_COMM_WORLD);
            // std::cout<<"Process["<<rank<<"] recv sig="<<signal<<" from Process["<<P_AFFICHAGE<<"]"<<std::endl;
        }

    }else if (rank==P_AFFICHAGE){
        bool animate=false;
        int Up = 2, Down = 3, Close=-1;
        std::size_t resx=800, resy=600;
        if (nargs>3)
        {
            resx = std::stoull(argv[2]);
            resy = std::stoull(argv[3]);
        }
        std::cout << "######## Vortex simultor ########" << std::endl << std::endl;
        std::cout << "Press P for play animation " << std::endl;
        std::cout << "Press S to stop animation" << std::endl;
        std::cout << "Press right cursor to advance step by step in time" << std::endl;
        std::cout << "Press down cursor to halve the time step" << std::endl;
        std::cout << "Press up cursor to double the time step" << std::endl;

        Graphisme::Screen myScreen( {resx,resy}, {grid.getLeftBottomVertex(), grid.getRightTopVertex()} );
        int t=0;float sum_fps=0;
        while (myScreen.isOpen())
        
        {
            auto start = std::chrono::system_clock::now();
            bool advance = false;
            // on inspecte tous les évènements de la fenêtre qui ont été émis depuis la précédente itération
            sf::Event event;
            while (myScreen.pollEvent(event))
            {

                // évènement "fermeture demandée" : on ferme la fenêtre
                if (event.type == sf::Event::Closed){
                    // MPI_Send(&Close, 1, MPI_INT, P_CALCUL_1, 200, MPI_COMM_WORLD);
                    // std::cout<<"Process["<<rank<<"] Bcast sig(Close)="<<Close<<std::endl;
                    MPI_Bcast(&Close, 1, MPI_INT, P_AFFICHAGE,MPI_COMM_WORLD);
                    myScreen.close();
                }
                if (event.type == sf::Event::Resized)
                {
                    // on met à jour la vue, avec la nouvelle taille de la fenêtre
                    myScreen.resize(event);
                }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::P)) animate = true;
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::S)) animate = false;
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up)) {
                    dt *= 2;
                    // std::cout<<"Process["<<rank<<"] Bcast sig(Up)="<<Up<<std::endl;
                    // MPI_Send(&Up, 1, MPI_INT, P_CALCUL_1, 200, MPI_COMM_WORLD);
                    MPI_Bcast(&Up, 1, MPI_INT, P_AFFICHAGE, MPI_COMM_WORLD);
                }
                
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down)) {                    
                    // std::cout<<"Process["<<rank<<"] Bcast sig(Down)="<<Down<<std::endl;
                    dt /= 2;
                    // MPI_Send(&Down, 1, MPI_INT, P_CALCUL_1, 200, MPI_COMM_WORLD);
                    MPI_Bcast(&Down, 1, MPI_INT, P_AFFICHAGE, MPI_COMM_WORLD);
                }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right)) advance = true;
            }
            if (animate | advance)
            {

                t++;
                // unsigned long nbOfPts, nbOfVts = vortices.numberOfVortices();
                MPI_Status status;
                // std::cout<<"Process["<<rank<<"] Bcast sig(isMobile)="<<isMobile<<std::endl;
                MPI_Bcast(&isMobile, 1, MPI_INT, P_AFFICHAGE, MPI_COMM_WORLD);
                if (isMobile)
                {
                    // unsigned long recvInfo[2];
                    // std::cout<<"Process["<<rank<<"] send sig(isMobile)="<<isMobile<<"to Process["<<P_CALCUL_1<<"]"<<std::endl;
                    // MPI_Send(&isMobile, 1, MPI_INT, P_CALCUL_1, 200, MPI_COMM_WORLD);
                    // MPI_Recv(&nbOfPts,  2, MPI_UNSIGNED_LONG, 1, 399, MPI_COMM_WORLD, &status);

                    // std::cout<<"Process["<<rank<<"] recv nb Of Pts="<<nbOfPts<<" from Process["<<P_CALCUL_1<<"]"<<std::endl;
                    // std::cout<<"Process["<<rank<<"] recv nb Of Vortices="<<recvInfo[0]<<" from Process["<<P_CALCUL_1<<"]"<<std::endl;
                    std::vector<double> newCloudData(nbOfPts*2);
                    // std::cout<<"Process["<<rank<<"] prepare recv Cloud Data from Process["<<P_CALCUL_1<<"]"<<std::endl;
                    // std::cout<<rank<<": Before Allgatherv()"<<std::endl;
                    MPI_Allgatherv(NULL, recvcountsCloud[rank], MPI_DOUBLE, newCloudData.data(), recvcountsCloud, displsCloud,  MPI_DOUBLE, global);

                    // MPI_Recv(newCloudData.data(), nbOfPts*2, MPI_DOUBLE, 1, 400, MPI_COMM_WORLD, &status);
                    // std::cout<<"Process["<<rank<<"] recv Cloud Data from Process["<<P_CALCUL_1<<"]"<<std::endl;
                    std::vector<double> newVorticeData(nbOfVts* 3);
                    MPI_Recv(newVorticeData.data(), nbOfVts*3, MPI_DOUBLE, 1, 401, MPI_COMM_WORLD, &status);
                    // std::cout<<"Process["<<rank<<"] recv Vortice Data from Process["<<P_CALCUL_1<<"]"<<std::endl;
                    for(unsigned long i=0; i<nbOfVts;i++){
                        vortices.setVortex(i, Geometry::Point<double>{newVorticeData[3*i], newVorticeData[3*i+1]}, newVorticeData[3*i+2]);
                    }
                    // cloud = Numeric::solve_RK4_movable_vortices(dt, grid, vortices, cloud);
                    Geometry::CloudOfPoints newCloud(newCloudData);
                    cloud = newCloud;

                    std::vector<double> newGridData(nbOfGridData*2);
                    MPI_Recv(newGridData.data(), nbOfGridData*2, MPI_DOUBLE, 1, 402, MPI_COMM_WORLD, &status);
                    // std::cout<<"P["<<rank<<"] recv grid data, numbre of values = "<<nbOfGridData<<std::endl; 
                    grid.SetVelocityFieldVal(newGridData);

                    
                    

                }
                else
                {
                    
                    // unsigned long nbOfPts=cloud.numberOfPoints();;
                    // unsigned long subnbOfPts = nbOfPts/N_calcul;
                    std::vector<double> newCloudData(nbOfPts*2);
                    // MPI_Recv(newCloudData.data(), nbOfPts*2, MPI_DOUBLE, 1, 400, MPI_COMM_WORLD, &status);
                    MPI_Allgatherv(NULL, recvcountsCloud[rank], MPI_DOUBLE, newCloudData.data(), recvcountsCloud, displsCloud,  MPI_DOUBLE, global);

                    Geometry::CloudOfPoints newCloud(newCloudData);
                    cloud = newCloud;
                }
            }
            myScreen.clear(sf::Color::Black);
            std::string strDt = std::string("Time step : ") + std::to_string(dt);
            myScreen.drawText(strDt, Geometry::Point<double>{50, double(myScreen.getGeometry().second-96)});
            myScreen.displayVelocityField(grid, vortices);
            myScreen.displayParticles(grid, vortices, cloud);
            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> diff = end - start;
            std::string str_fps = std::string("FPS : ") + std::to_string(1./diff.count());
            if(animate &&t<=100){
                // std::cout<<"t="<<t<<", fps= "<<std::to_string(1./diff.count())<<std::endl;
                sum_fps +=1./diff.count();

            }
            if(t==100)std::cout<<sum_fps/100<<std::endl;
            myScreen.drawText(str_fps, Geometry::Point<double>{300, double(myScreen.getGeometry().second-96)});
            myScreen.display();
            
            
        }

    }
   

    MPI_Finalize();
    return EXIT_SUCCESS;
    
 }