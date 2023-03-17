#include <iostream>
#include "runge_kutta.hpp"
#include "cartesian_grid_of_speed.hpp"
#include <omp.h>
using namespace Geometry;
// #include <mpi.h>

Geometry::CloudOfPoints
Numeric::solve_RK4_fixed_vortices( double dt, CartesianGridOfSpeed const& t_velocity, Geometry::CloudOfPoints const& t_points )
{
    constexpr double onesixth = 1./6.;
    using vector = Simulation::Vortices::vector;
    using point  = Simulation::Vortices::point;

    Geometry::CloudOfPoints newCloud(t_points.numberOfPoints());
    // On ne bouge que les points :
    #pragma omp parallel for
    for ( std::size_t iPoint=0; iPoint<t_points.numberOfPoints(); ++iPoint)
    {
        point  p = t_points[iPoint];
        vector v1 = t_velocity.computeVelocityFor(p);
        point p1 = p + 0.5*dt*v1;
        p1 = t_velocity.updatePosition(p1);
        vector v2 = t_velocity.computeVelocityFor(p1);
        point p2 = p + 0.5*dt*v2;
        p2 = t_velocity.updatePosition(p2);
        vector v3 = t_velocity.computeVelocityFor(p2);
        point p3 = p + dt*v3;
        p3 = t_velocity.updatePosition(p3);
        vector v4 = t_velocity.computeVelocityFor(p3);
        newCloud[iPoint] = t_velocity.updatePosition(p + onesixth*dt*(v1+2.*v2+2.*v3+v4));
    }
    return newCloud;
}

Geometry::CloudOfPoints
Numeric::solve_RK4_movable_vortices( double dt, CartesianGridOfSpeed& t_velocity, 
                                     Simulation::Vortices& t_vortices, 
                                     Geometry::CloudOfPoints const& t_points )
{
    constexpr double onesixth = 1./6.;
    using vector = Simulation::Vortices::vector;
    using point  = Simulation::Vortices::point;

    Geometry::CloudOfPoints newCloud(t_points.numberOfPoints());
    // On ne bouge que les points :
#pragma omp parallel for
    for ( std::size_t iPoint=0; iPoint<t_points.numberOfPoints(); ++iPoint)
    {
        point  p = t_points[iPoint];
        vector v1 = t_velocity.computeVelocityFor(p);
        point p1 = p + 0.5*dt*v1;
        p1 = t_velocity.updatePosition(p1);
        vector v2 = t_velocity.computeVelocityFor(p1);
        point p2 = p + 0.5*dt*v2;
        p2 = t_velocity.updatePosition(p2);
        vector v3 = t_velocity.computeVelocityFor(p2);
        point p3 = p + dt*v3;
        p3 = t_velocity.updatePosition(p3);
        vector v4 = t_velocity.computeVelocityFor(p3);
        newCloud[iPoint] = t_velocity.updatePosition(p + onesixth*dt*(v1+2.*v2+2.*v3+v4));
    }
    std::vector<point> newVortexCenter;
    newVortexCenter.reserve(t_vortices.numberOfVortices());
#pragma omp parallel for
    for (std::size_t iVortex=0; iVortex<t_vortices.numberOfVortices(); ++iVortex)
    {
        point p = t_vortices.getCenter(iVortex);
        vector v1 = t_vortices.computeSpeed(p);
        point p1 = p + 0.5*dt*v1;
        p1 = t_velocity.updatePosition(p1);
        vector v2 = t_vortices.computeSpeed(p1);
        point p2 = p + 0.5*dt*v2;
        p2 = t_velocity.updatePosition(p2);
        vector v3 = t_vortices.computeSpeed(p2);
        point p3 = p + dt*v3;
        p3 = t_velocity.updatePosition(p3);
        vector v4 = t_vortices.computeSpeed(p3);
        newVortexCenter.emplace_back(t_velocity.updatePosition(p + onesixth*dt*(v1+2.*v2+2.*v3+v4)));
    }
    // #pragma omp parallel for
    for (std::size_t iVortex=0; iVortex<t_vortices.numberOfVortices(); ++iVortex)
    {
        t_vortices.setVortex(iVortex, newVortexCenter[iVortex], 
                             t_vortices.getIntensity(iVortex));
    }
    t_velocity.updateVelocityField(t_vortices);
    return newCloud;

}



/*
    multi-processor
*/



// std::vector<Geometry::Point<double>>
// Numeric::sub_solve_RK4_vortices_new_cloud( double dt, CartesianGridOfSpeed& t_velocity, 
//                                      Geometry::CloudOfPoints const& t_points, int id, int N )
// {
//     constexpr double onesixth = 1./6.;
//     using vector = Simulation::Vortices::vector;
//     using point  = Simulation::Vortices::point;
//     unsigned long sub_nbOfPts = t_points.numberOfPoints()/N;
//     if(id==N){
//         sub_nbOfPts = t_points.numberOfPoints() - (N-1)*sub_nbOfPts;
//     }
//     std::vector<Geometry::Point<double>> newSubCloud(sub_nbOfPts);
//     for ( std::size_t idx=0, iPoint=sub_nbOfPts*id; idx<sub_nbOfPts; ++idx, ++iPoint)
//     {
//         point  p = t_points[iPoint];
//         vector v1 = t_velocity.computeVelocityFor(p);
//         point p1 = p + 0.5*dt*v1;
//         p1 = t_velocity.updatePosition(p1);
//         vector v2 = t_velocity.computeVelocityFor(p1);
//         point p2 = p + 0.5*dt*v2;
//         p2 = t_velocity.updatePosition(p2);
//         vector v3 = t_velocity.computeVelocityFor(p2);
//         point p3 = p + dt*v3;
//         p3 = t_velocity.updatePosition(p3);
//         vector v4 = t_velocity.computeVelocityFor(p3);
//         newSubCloud[idx] = t_velocity.updatePosition(p + onesixth*dt*(v1+2.*v2+2.*v3+v4));
//     }
//     return newSubCloud;
// }



// Geometry::CloudOfPoints
// Numeric::sub_solve_RK4_fixed_vortices( double dt, CartesianGridOfSpeed const& t_velocity, Geometry::CloudOfPoints const& t_points, int id, int N )
// {
// // id = rank of the process among all processes for calculation  (1 ~ N)
// // N  = number of processes of calculation
//     constexpr double onesixth = 1./6.;
//     using vector = Simulation::Vortices::vector;
//     using point  = Simulation::Vortices::point;

//     unsigned long sub_nbOfPts = t_points.numberOfPoints()/N;
//     if(id==N){
//         sub_nbOfPts = t_points.numberOfPoints() - (N-1)*sub_nbOfPts;
//     }

//     Geometry::CloudOfPoints newCloud(t_points.numberOfPoints());
//     // On ne bouge que les points :
//     // #pragma omp parallel for
//     for ( std::size_t idx = 0, iPoint=id*sub_nbOfPts; idx<t_points.numberOfPoints(); ++iPoint, ++idx)
//     {
//         point  p = t_points[iPoint];
//         vector v1 = t_velocity.computeVelocityFor(p);
//         point p1 = p + 0.5*dt*v1;
//         p1 = t_velocity.updatePosition(p1);
//         vector v2 = t_velocity.computeVelocityFor(p1);
//         point p2 = p + 0.5*dt*v2;
//         p2 = t_velocity.updatePosition(p2);
//         vector v3 = t_velocity.computeVelocityFor(p2);
//         point p3 = p + dt*v3;
//         p3 = t_velocity.updatePosition(p3);
//         vector v4 = t_velocity.computeVelocityFor(p3);
//         newCloud[iPoint] = t_velocity.updatePosition(p + onesixth*dt*(v1+2.*v2+2.*v3+v4));
//     }
//     return newCloud;
// }


// Geometry::CloudOfPoints
// Numeric::sub_solve_RK4_movable_vortices( double dt, CartesianGridOfSpeed& t_velocity, 
//                                      Simulation::Vortices& t_vortices, 
//                                      Geometry::CloudOfPoints const& t_points, int id, int N )
// // id = rank of the process among all processes for calculation  (1 ~ N)
// // N  = number of processes of calculation
// {
//     using vector = Simulation::Vortices::vector;
//     using point  = Simulation::Vortices::point;
//     std::vector<point> newSubCloud = sub_solve_RK4_vortices_new_cloud(dt, t_velocity, t_points, id, N );



//     std::vector<point> newVortexCenter;
//     newVortexCenter.reserve(t_vortices.numberOfVortices());
//     for (std::size_t iVortex=0; iVortex<t_vortices.numberOfVortices(); ++iVortex)
//     {
//         point p = t_vortices.getCenter(iVortex);
//         vector v1 = t_vortices.computeSpeed(p);
//         point p1 = p + 0.5*dt*v1;
//         p1 = t_velocity.updatePosition(p1);
//         vector v2 = t_vortices.computeSpeed(p1);
//         point p2 = p + 0.5*dt*v2;
//         p2 = t_velocity.updatePosition(p2);
//         vector v3 = t_vortices.computeSpeed(p2);
//         point p3 = p + dt*v3;
//         p3 = t_velocity.updatePosition(p3);
//         vector v4 = t_vortices.computeSpeed(p3);
//         newVortexCenter.emplace_back(t_velocity.updatePosition(p + onesixth*dt*(v1+2.*v2+2.*v3+v4)));
//     }
//     for (std::size_t iVortex=0; iVortex<t_vortices.numberOfVortices(); ++iVortex)
//     {
//         t_vortices.setVortex(iVortex, newVortexCenter[iVortex], 
//                              t_vortices.getIntensity(iVortex));
//     }
//     t_velocity.updateVelocityField(t_vortices);
//     return newCloud;

// }