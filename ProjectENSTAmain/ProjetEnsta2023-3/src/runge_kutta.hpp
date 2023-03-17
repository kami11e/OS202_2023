#ifndef _NUMERIC_RUNGE_KUTTA_HPP_
#define _NUMERIC_RUNGE_KUTTA_HPP_
#include <utility>
#include "vortex.hpp"
#include "cloud_of_points.hpp"
#include "cartesian_grid_of_speed.hpp"

namespace Numeric 
{

    Geometry::CloudOfPoints
    solve_RK4_fixed_vortices( double dt, CartesianGridOfSpeed const& speed, 
                              Geometry::CloudOfPoints const& t_points );

    Geometry::CloudOfPoints
    solve_RK4_movable_vortices( double dt, CartesianGridOfSpeed& t_velocity, 
                                Simulation::Vortices& t_vortices, 
                                Geometry::CloudOfPoints const& t_points );

    // Geometry::CloudOfPoints
    // Numeric::sub_solve_RK4_fixed_vortices( double dt, CartesianGridOfSpeed const& t_velocity, 
    //                                     Geometry::CloudOfPoints const& t_points, int id, int N );

    // Geometry::CloudOfPoints
    // sub_solve_RK4_movable_vortices( double dt, CartesianGridOfSpeed& t_velocity, 
    //                                  Simulation::Vortices& t_vortices, 
    //                                  Geometry::CloudOfPoints const& t_points, int id, int N );

    Geometry::CloudOfPoints
    sub_solve_RK4_vortices_new_cloud( double dt, CartesianGridOfSpeed& t_velocity, 
                                     Geometry::CloudOfPoints const& t_points, int id, unsigned long subnbOfPts0, unsigned long subnbOfPts );
    void
    solve_RK4_movable_vortices_update_vortices(double dt, CartesianGridOfSpeed& t_velocity, 
                                        Simulation::Vortices& t_vortices,
                                        Geometry::CloudOfPoints const& t_points);
}

#endif