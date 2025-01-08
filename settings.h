// Copyright 2023 Will Snider

#pragma once
#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
#include <utility>
#include <random>
#include <chrono>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

namespace settings {
    // Freq Models: "Fornax (new, from Lelli et al.)", "Leo II (new, from Lelli et al.)", "Sculptor", "ToyModel"
    const string modelName = "Fornax (new, from Lelli et al.)";    // Set of characteristics inherited by the galaxy
    const int N = 25;                       // Number of bodies in the simulation
    const double Tmax = 1000;               // Runtime in My
    const double minstep = 0.005;            // Minimum length of a timestep (no higher than 0.01)
    const double initStep = minstep;        // Timestep stars are initialized with
    const double outputTime = 1;            // Time between outputs to the data file
    const double consoleWrites = floor(sqrt(25+pow(Tmax/10,2)));  // Number of console writes (for efficiency)                  
    const double bins = floor(sqrt(N));        // Groups for dispersion profile averages
    const double massPerc = 0.98;           // Mass percent radius defining maximum population range

    /* Switches
            dooubles: off = -1
            phi: azimuthal angle from +x axis CCW in xy plane 
    */
    const bool g_ratios = true;             // overrides provided model and uses standard toy model with given g ratios (fixed M & R_h)
        const double toy_mass = 2e8;         // standard dSph mass for toy model (2 * 10^8 solar masses)
        const double toy_hostR = 1e6;        // standard host distance for toy model (10^6 pc)
    const bool plCir = false;                // planar circular orbits in xy plane (only works for Newtonian)
    const double uniform_r = -1;           // stars initialized with same radius  
    const double uniform_phi = -1;          // stars initialized with same phi, use with plCir
    const double lin_dist_r = -1;           // stars radius dispersed evenly (to given max radius)
    const int lin_dist_phi = -1;             // stars phi dispersed evenly (in given number of groups)
    const bool CenterOfMass = false;         // tracks center of mass position 
    const bool blackHole = false;           // include central blackhole with mass mBlack
        const double mBlack = 1e6;          // mass of central blackhole
    const bool MOND = false;                // Modified Newtonian Dymanics switch
    const bool EFE = true;                 // External Field Effect switch
    const bool STVG = false;                 // experimental Scalar-Tensor-Vector correction based on Moffat & Toth 2023
    const bool trackTidalR = true;          // records pos mag of furthest bound star at every output

    // Dispersion Controls
    const int axis = 1;                     // Viewing Axis: 1 = x, 2 = y, 3 = z
    const double trunc_dist = 1e10;         // Distance to Truncate Dispersion & COM (none = 1e10)

    // Output Control
    const string simOutput = "results.txt";
    const string dispOutput = "disp_dat.txt";
    const bool run_dispersion = false;
    const bool pos_out = true, vel_out = false, acc_out = false;
        /* Format Options :
        *   0 = Animation [time N / ID x y z...]
        *   1 = Statistics [ID x y z / ...]
        *   2 = Animation w/ COM [time N COM / ID x y z...]
        *   3 = COM r and v [time x y z vx vy vz / ...]
        */
    const int format = 0;
}   // namespace set

namespace consts {
    const double years = 3600.0 * 24.0 * 365.25;    // seconds in a year
    const double AU = 1.495978707e11;               // AU in MKS
    const double G_mks = 6.672e-11;                 // gravitation constant (MKS)
    const double Msun = 1.98911e30;                 // solar mass in MKS(kg)
    const double pc = 3.0856776e16;                 // meters in a parsec
    const double G = 0.004498535262;                // pc**3/(My**2 M_sun)
    const double pi = acos(-1);
    const double tappr = 10.0;                      // about 1 million year
    const double a_mond = 3.88;                     // MOND acc. const in pc/My**2
    const double eta = 0.0005;                      // timestep adjustable parameter

    // random_device rd;
    // default_random_engine generator(rd());
    // uniform_real_distribution<double> distribution(0.0, 1.0);
}  // namespace consts

#endif