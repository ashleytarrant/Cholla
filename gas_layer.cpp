#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "global.h"
#include "grid3D.h"
#include "mpi_routines.h"
#include "error_handling.h"
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>


void Grid3D::Gas_Layer_2D(Real n, Real vx, Real vy, Real vz, Real P)
{
    int i, j, id;
    int istart, jstart, kstart, iend, jend, kend; 

    Real x_pos, y_pos, z_pos;
    
    Real mu = 0.6;
    Real rho, T;

    rho = n*mu*MP / DENSITY_UNIT;
    vx  = vx * 1e5 / VELOCITY_UNIT;
    vy  = vy * 1e5 / VELOCITY_UNIT;
    vz  = vz * 1e5 / VELOCITY_UNIT;
    P   = P*KB / PRESSURE_UNIT;

    istart = H.n_ghost;
    iend   = H.nx-H.n_ghost;
    if (H.ny > 1) {
        jstart = H.n_ghost;
        jend   = H.ny-H.n_ghost;
    }
    else {
        jstart = 0;
        jend   = H.ny;
    }
    if (H.nz > 1) {
        kstart = H.n_ghost;
        kend   = H.nz-H.n_ghost;
    }
    else {
        kstart = 0;
        kend   = H.nz;
    } 

    // set the initial values of the box

    for (j=jstart; j<jend; j++) {
        for (i=istart; i<iend; i++) { 

            //get cell index
            id = i + j*H.nx;

            // get the centered x and y positions 
            Get_Position(i, j, H.n_ghost, &x_pos, &y_pos, &z_pos);

            // set constant initial states
            C.density[id]    = rho;
            C.momentum_x[id] = rho*vx;
            C.momentum_y[id] = rho*vy;
            C.momentum_z[id] = rho*vz;
            C.Energy[id]     = P/(gama-1.0) + 0.5*rho*(vx*vx + vy*vy + vz*vz);
            #ifdef DE
            C.GasEnergy[id]  = P/(gama-1.0);
            #endif

            Real n_layer, T_layer, P_layer, rho_layer, rho_layer_CU, P_layer_CU;
            Real z_layer, z_layerwidth, z_layer_max, z_layer_min; 

            // Define a specific region (layer) and change the conserved values accordingly
            n_layer = 8000.0 ; // 1/cm^3 
            rho_layer = n_layer * mu * MP ; // CGS units 
            rho_layer_CU = rho_layer / (MASS_UNIT/pow(LENGTH_UNIT,3.0)) ; // Code units Msun/kpc^3
            T_layer = 10000.0 ; // K 
            P_layer = n_layer * T_layer * KB ; // 
            P_layer_CU = P_layer / PRESSURE_UNIT ; 

            z_layer = 2.0  ; // layer starts at 2 kpc  
            z_layerwidth = 0.4 ; // 400 pc thick
            z_layer_max = z_layer + 0.5 * z_layerwidth ; 
            z_layer_min = z_layer - 0.5 * z_layerwidth ; 
            

           if (z_layer_min <= y_pos && y_pos <= z_layer_max)   
            {
                // Set up layer, changed the conserved variables 
                C.density[id]    = rho_layer_CU;
                C.momentum_x[id] = 0.0;
                C.momentum_y[id] = 0.0;
                C.momentum_z[id] = 0.0;
                C.Energy[id]     = P_layer_CU/(gama-1.0);
                #ifdef DE
                C.GasEnergy[id]  = P_layer_CU/(gama-1.0);
                #endif
            }
        }
    }

}

/*! \fn void Wind_Boundary()
* \brief Apply a wind boundary to the -x face */
void Grid3D::Wind_Boundary()
{
    int i, j, id;
    Real x_pos, y_pos, z_pos;
    Real n, rho, vx, vy, vz, P, T;

    // for now, hard code values in
    Real mu = 0.6;
    n = 0.1;  // 1/cm^-3
    vx = 0.0; // km/s
    vy = 1000.0;
    vz = 0.0;
    P = 1e1;

    // convert values to code units
    rho = n*mu*MP / DENSITY_UNIT;
    vx = vx*1e5 / VELOCITY_UNIT;
    vy = vy*1e5 / VELOCITY_UNIT;
    vz = vz*1e5 / VELOCITY_UNIT;
    P = P*KB / PRESSURE_UNIT;

    // set exact boundaries on the -x face
    for (i=0; i<H.nx; i++) {
        for (j=0; j<H.n_ghost; j++) {

            // get cell id
            id = i + j*H.nx;
            // get centered x, y, z positions
            Get_Position(i, j, H.n_ghost, &x_pos, &y_pos, &z_pos);

            // set the conserved quantities
            C.density[id] = rho;
            C.momentum_x[id] = vx*C.density[id];
            C.momentum_y[id] = vy*C.density[id];
            C.momentum_z[id] = vz*C.density[id];
            C.Energy[id] = P/(gama-1.0) + 0.5*C.density[id]*(vx*vx + vy*vy + vz*vz);
            #ifdef DE
            C.GasEnergy[id] = P/(gama-1.0);
            #endif //DE
        }
    }



}
