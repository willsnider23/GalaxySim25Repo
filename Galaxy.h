// Copyright 2023 Will Snider

#pragma once
#ifndef GALAXY_HELPER_H
#define GALAXY_HELPER_H

#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <math.h>
#include <random>
#include "Star.h"
#include "settings.h"

using namespace std;
using namespace std::string_literals;
using PosMat = vector<vector<double>>;
using Model = unordered_map<string, double>;

class Galaxy {
    vector<Star> population;
    double baryonCloudMass;
    double r_half, a_s, r_tidal;
    double host_R, host_M;
    PosMat centerOfMass = { {}, {} };

public:
    Galaxy(Model& model);
    ~Galaxy() {}

// Getters
    double
    getMass() { return baryonCloudMass; }
    double
    getRHalf() { return r_half; }
    double
    getRTidal() { return r_tidal; }
    double
    getHostDist() { return host_R; }
    double 
    getHostMass() { return host_M; }
    Star&
    getStar(int ID);
    vector<Star>&
    getPopulation() { return population; }
    int
    getSize() { return population.size(); }
    PosMat&
    getCOM() { return centerOfMass; }
    bool
    isUniformTime();
    double
    calcAnisotropyFactor();
    double
    getGeff();
// Setters
    void
    setPopulation(vector<Star> pop) { population = pop; }
// Mutators
    void
    calcCOM();
    void
    wrangleStars(double time);
    void
    HITS(Star& s,double time);
};
 
#endif
