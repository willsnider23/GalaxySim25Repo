// Copyright 2023 Will Snider

#include "Galaxy.h"
#include "settings.h"

using PosMat = vector<vector<double>>;

vector<Star>
createPop(double r_half, double baryonCloudMass, double host_R, 
          double host_M, PosMat& centerOfMass) {
    vector<Star> stars;

    // Initialize N bodies, push into pop
    for (int idx = 0; idx < settings::N; idx++) {
        Star star(idx, r_half, baryonCloudMass, host_R, host_M, centerOfMass);
        stars.push_back(star);
    }        
    return stars;
}

double
calcTcross(double r_half, double baryonCloudMass) {
    double a_s = r_half * sqrt(pow(2.0, (2.0 / 3.0)) - 1.0);
    double coefficient = 128 * sqrt(6) / (9 * pow(consts::pi, (1.5)));
    return coefficient * pow(a_s, 1.5) / sqrt(baryonCloudMass * consts::G);
}

// Constructor
Galaxy::Galaxy(Model& model) {
        baryonCloudMass = model["mass"];
        r_half = model["r_half"];
        a_s = r_half * sqrt(pow(2.0, (2.0 / 3.0)) - 1.0);
        r_tidal = model["r_tidal"];
        if (settings::extField) { 
            host_R = model["R_mw"]; 
            host_M = model["M_mw"];
        } else {
            host_R = INFINITY;
            host_M = 0;
        }
        Tcross = calcTcross(r_half, baryonCloudMass);
        COMa_and_adot();
        population = createPop(r_half, baryonCloudMass, host_R, host_M, centerOfMass);
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
Galaxy::isUniformTime() const {
    double time = population[0].getTime();
    for (Star s : population) {
        if (s.getTime() > time + settings::minstep ||
            s.getTime() < time - settings::minstep)
            return false;
    }
    return true;
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
Galaxy::getGeff() const {
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
        vector<double> pos = s.getPos();
        vector<double> vel = s.getVel();
        for (int i = 0; i < 3; i++) {
            if (settings::trunc_dist == -1 || calcMag(pos) < settings::trunc_dist) {
                sumR[i] += s.getPos()[i];
                sumV[i] += s.getVel()[i];
            }
        }
    }
    centerOfMass[0] = { sumR[0] / settings::N, sumR[1] / settings::N, sumR[2] / settings::N };
    centerOfMass[1] = { sumV[0] / settings::N, sumV[1] / settings::N, sumV[2] / settings::N };
}

// Newtonian acceleration and jerk of the center of mass
// copied from Star external field submethod
void
Galaxy::COMa_and_adot() {
    centerOfMass[2] = { 0,0,0 };
    centerOfMass[3] = { 0,0,0 };
    // currently fixed at origin & static (have it follow calculated COM later?)
    vector<double> pos_e = { -host_R, 0, 0 };
    vector<double> vel(3, 0);
    double r_e = calcMag(pos_e);
    double r_e2 = pow(r_e, 2);
    double r_e3 = pow(r_e, 3);
    double factor1 = consts::G * host_M / r_e3;
    double re_dot_v = dotProduct(pos_e, vel);
    double factor2 = factor1 * 3.0 * re_dot_v / r_e2;

    /* note here we can use velocity of star, vel, since the host is not moving.
       If the host moves, this will have to be the relative velocity. */
    for (int i = 0; i < 3; i++) {
        centerOfMass[2][i] = -pos_e[i] * factor1;
        centerOfMass[3][i] = (-vel[i] * factor1) + (pos_e[i] * factor2);
    }
}

void 
Galaxy::calcSkewness() {
    // get list of all x-positions in population
    vector<double> posList;
    for (Star& s : population) {
        if (s.isBound(getGeff(), getMass(), getRHalf() * sqrt(pow(2.0, (2.0 / 3.0)) - 1.0)))
            posList.push_back(s.getPos()[0]);
    }
    sort(posList.begin(), posList.end());

    // finding 10th, 50th, and 90th percentile values for Kelly's coeff
    double P_ten = posList[floor(posList.size() * 0.1)];
    double P_fifty = posList[floor(posList.size() * 0.5)];
    double P_ninety = posList[floor(posList.size() * 0.9)];

    // Kelly's Coefficient of Skewness
    skewness = (P_ninety - 2*P_fifty + P_ten) / (P_ninety - P_ten);
}

void
Galaxy::wrangleStars(double time) {
    for (Star& s : population) {
        if (!s.isFrozen())
            this->HITS(s, time);
    }
}

void
Galaxy::HITS(Star& s, double time) {
    if (s.getAfterTimestep() > time) {
        s.predict(time - s.getTime());
        s.a_and_adot(r_half, baryonCloudMass, host_R, host_M, true, centerOfMass);
        s.correct(time - s.getTime());
    }
}
