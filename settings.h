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
    const int N = 100;                       // Number of bodies in the simulation
    const double TmaxConst = 1000;           // Fixed runtime in My
    const double minstep = 0.005;            // Minimum length of a timestep (no higher than 0.01)
    const double initStep = minstep;        // Timestep stars are initialized with
    const double outputTime = 1;            // Time between outputs to the data file
    const double consoleWrites = floor(sqrt(25+pow(TmaxConst/10,2)));  // Number of console writes (for efficiency)                  
    const double bins = floor(sqrt(N));        // Groups for dispersion profile averages
    const double massPerc = 0.98;           // Mass percent radius defining maximum population range
    const double runs = 3;                  // Number of simulations to run sequentially (currently only works in conjunction w/ g_ratios)

    /* Switches
            doubles: off = -1
            phi: azimuthal angle from +x axis CCW in xy plane 
    */
    const bool MOND = true;                 // Modified Newtonian Dymanics switch
    const bool extField = true;              // External Field switch
        const bool diverging = true;         // divergence of ext field
    const bool STVG = false;                 // experimental Scalar-Tensor-Vector correction based on Moffat & Toth 2023
    const bool g_ratios = true;              // overrides provided model and uses standard toy model with given g ratios (fixed M & R_h)
        const double toy_mass = 2e8;         // standard dSph mass for toy model (2 * 10^8 solar masses)
        const double toy_hostR = 1e6;        // standard host distance for toy model (10^6 pc)
    const bool blackHole = false;            // include central blackhole with mass mBlack
        const double mBlack = 1e6;           // mass of central blackhole
    const bool constRuntime = false;         // fixed (true) or scale dependent (false) runtime
        const double crossings = 5;		 // number of crossing times used for scale dependent runtime
    /* Integrated Acceleration:     
    *   0 = gi + ge_star = g_total (orbiting)
    *   1 = g_total - ge_COM       (tidal)
    */ 
    const int integ_acc = 1;
    // Data Tracking
    const bool CenterOfMass = false;          // tracks center of mass
    const bool trackTidalR = false;          // records pos mag of furthest bound star at every output
    const bool trackSkews = true;            // records skew param at every console output
    const bool freezeStrays = true;          // fixes problematic stars in space with zero velocity
        const double freezeR = 5e4;          // Dist from host when stars freeze - 50 kpc
    // Initial Distribution
    const bool plCir = false;                // planar circular orbits in xy plane (only works for Newtonian)
    const double uniform_r = -1;             // stars initialized with same radius  
    const double uniform_phi = -1;           // stars initialized with same phi, use with plCir
    const double lin_dist_r = -1;            // stars radius dispersed evenly (to given max radius)
    const int lin_dist_phi = -1;             // stars phi dispersed evenly (in given number of groups)

    // Dispersion Controls
    const int axis = 1;                      // Viewing Axis: 1 = x, 2 = y, 3 = z
    const double trunc_dist = -1;            // Distance to Truncate Dispersion & COM 
    // Output Control
    const string simOutput = "results.txt";
    const string dispOutput = "disp_dat.txt";
    const string COM_Output = "centerOmass.txt";
    const string tidalOutput = "tidalR.txt";
    const string skewOutput = "skews.txt";
    const bool run_dispersion = false;
    const bool pos_out = true, vel_out = true, acc_out = false;
        /* Format Options :
        *   0 = Animation [time N r_half / ID x y z...]
        *   1 = Statistics [ID x y z / ...]
        *   2 = Animation w/ COM [time N COM / ID x y z...]
        *   3 = COM r and v [time x y z vx vy vz / ...]
        *   4 = Animation - COM vel [time N COM_v / ID x y z...]
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

namespace {
    double dotProduct(const vector<double>& a, const vector<double>& b) {
        double sum = 0;
        for (int i = 0; i < 3; i++) sum += a[i] * b[i];
        return sum;
    }

    vector<double> crossProduct(const vector<double>& a, const vector<double>& b) {
        return { (b[2] * a[1]) - (a[2] * b[1]),
                 (b[0] * a[2]) - (a[0] * b[2]),
                 (b[1] * a[0]) - (a[1] * b[0]) };
    }

    double calcMag(vector<double>& vec) {
        double mag = 0;
        for (double val : vec) { mag += pow(val, 2); }
        return sqrt(mag);
    }
}

#endif