// Copyright 2023 Will Snider

#include "Galaxy.h"
#include "settings.h"

vector<Star>
initialConditions(double r_half, double cloudMass, double hostR, double hostM) {
    vector<Star> stars;

    // Initialize N bodies, push into pop
    for (int idx = 0; idx < settings::N; idx++) {
        Star star(idx, r_half, cloudMass, hostR, hostM);
        stars.push_back(star);
    }        
    return stars;
}

// Constructor
Galaxy::Galaxy(Model& model) {
        baryonCloudMass = model["mass"];
        r_half = model["r_half"];
        a_s = r_half * sqrt(pow(2.0, (2.0 / 3.0)) - 1.0);
        r_tidal = model["r_tidal"];
        if (settings::EFE) { 
            host_R = model["R_mw"]; 
            host_M = model["M_mw"];
        } else {
            host_R = INFINITY;
            host_M = 0;
        }
        population = initialConditions(r_half, baryonCloudMass, host_R, host_M);
        calcCOM();
}

// Getters
Star&
Galaxy::getStar(int ID) {
    if (ID < static_cast<int>(population.size()))
        return population.at(ID);
    // else
        // throw std::invalid_argument("Not a valid star ID");
}

bool
Galaxy::isUniformTime() {
    double time = population[0].getTime();
    for (Star s : population) {
        if (s.getTime() > time + settings::minstep ||
            s.getTime() < time - settings::minstep)
            return false;
    }
    return true;
}

double
calcMag(vector<double> vec) {
    double mag = 0;
    for (double val : vec) { mag += pow(val, 2); }
    return sqrt(mag);
}

double
Galaxy::calcAnisotropyFactor() {
    double r_sum = 0, r2_sum = 0, t_sum = 0, t2_sum = 0;
    for (Star& s : population) {
        vector<double> vel_sph = s.getVelSph();
        double tangentSpeed = vel_sph[1]; // calcMag({ vel_sph[1], vel_sph[2] });
        r_sum += vel_sph[0];
        r2_sum += pow(vel_sph[0], 2);
        t_sum += tangentSpeed;
        t2_sum += pow(tangentSpeed, 2);
    }
    double r_avg = r_sum / population.size();
    double r2_avg = r2_sum / population.size();
    double t_avg = t_sum / population.size();
    double t2_avg = t2_sum / population.size();

    double sigma2_r = r2_avg - pow(r_avg, 2);
    double sigma2_t = t2_avg - pow(t_avg, 2);
    return 1 - (sigma2_t / sigma2_r);
}

// G_eff from external field mag at r_half for v_esc
double
Galaxy::getGeff() {
    double re = r_half - host_R;
    double gne = re * consts::G * host_M / pow(re, 3);
    double nu = 1.0 / (1.0 - exp(-sqrt(gne / consts::a_mond)));
    double ge = nu * gne;
    return consts::a_mond * consts::G / ge;
}

// Mutators
void
Galaxy::calcCOM() {
    vector<double> sumR = { 0, 0, 0 }, sumV = { 0, 0, 0 };
    for (Star& s : population)
    {
        for (int i = 0; i < 3; i++) {
            vector<double> pos = s.getPos();
            vector<double> vel = s.getVel();
            if (calcMag(pos) < settings::trunc_dist) {
                sumR[i] += s.getPos()[i];
                sumV[i] += s.getVel()[i];
            }
        }
    }
    centerOfMass[0] = { sumR[0] / settings::N, sumR[1] / settings::N, sumR[2] / settings::N };
    centerOfMass[1] = { sumV[0] / settings::N, sumV[1] / settings::N, sumV[2] / settings::N };
}

void
Galaxy::wrangleStars(double time) {
    for (Star& s : population) {
        this->HITS(s, time);
    }
}

void
Galaxy::HITS(Star& s, double time) {
    if (s.getAfterTimestep() > time) {
        s.predict(time - s.getTime());
        s.a_and_adot(r_half, baryonCloudMass, host_R, host_M, true);
        s.correct(time - s.getTime());
    }
}
