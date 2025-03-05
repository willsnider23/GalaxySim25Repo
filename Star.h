// Copyright 2023 Will Snider

#pragma once
#ifndef STAR_HELPER_H
#define STAR_HELPER_H

#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <math.h>
#include <random>
#include "settings.h"

using namespace std;
using namespace std::string_literals;
using PosMat = vector<vector<double>>;
using Model = unordered_map<string, double>;

class
Star {
    // position, velocity, acceleration, jerk, a2dot, and a3dot
    PosMat positionMatrix = {{}, {}, {}, {}, {}, {}};
    //    predicted r, v, a, j, and a2dot
    PosMat predMat = {{}, {}, {}, {}, {}};
    vector<vector<double>> log;
    double time, timestep;
    int ID, mass = 1;
    bool frozen = false;

public:
// Constuctor
    Star(int idx, double r_half, double cloudMass, double hostR, double hostM);
    Star(int idx, double r_half, double cloudMass, double hostR, double hostM, vector<double> initR, vector<double> initV);
// Getters
    int
    getID() const { return ID; }
    double
    getTime() const { return time;}
    double
    getTimestep() const { return timestep; }
    double
    getAfterTimestep() const { return time + timestep; }
    vector<double> 
    getPos() const { return positionMatrix[0];}
    vector<double> 
    getVel() { return positionMatrix[1]; }
    vector<double>
    getAcc() { return positionMatrix[2]; }
    vector<double>
    getJerk() { return positionMatrix[3]; }
    PosMat
    getPosMat() { return positionMatrix; }
    PosMat&
    getPredMat() { return predMat; }
    vector<double>&
    getLog(int idx) { return log[idx]; }
    vector<double>
    getPosSph();
    vector<double>
    getVelSph();
    bool
    isBound(double Geff, double host_M, double a_s);
    bool 
    isFrozen() const { return frozen ? true : false; }
// Setters
    void
    setTime(double t) { time = t; }
    void 
    setPos(vector<double>& pos) { positionMatrix[0] = pos; }
    void 
    setVel(vector<double>& vel) { positionMatrix[1] = vel; }
    void
    setPredMat(PosMat& preds) { predMat = preds; }
// Mutators
    void
    predict(double tstep);
    void
    a_and_adot(double r_half, double cloudMass, double hostR, double hostM, bool prediction);
    void
    correct(double step);
    void
    incrementTime();
    double
    updateStar(double rhalf, double gMass, double hostR, double hostM);
    void
    addLog();
};

#endif
// Hannha wuz here