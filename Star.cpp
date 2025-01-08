// Copyright 2023 Will Snider

#include "Star.h"
#include "settings.h"
#include <math.h>
#include <random>
#include <stdlib.h>
#include <time.h>

double 
LelliNu(double y) {
    return (1.0 / (1.0 - exp(-sqrt(y))));
}

double
LelliNu_Dot(double y, double nu, double g, double g_dot_gdot) {
    return (-(pow(nu, 2)*exp(-sqrt(y))) / (2.0 * sqrt(y))) * (g_dot_gdot) / (consts::a_mond * g);
}

// RNG [0-1]
double
RNG() {
    // static std::default_random_engine rng;
    // std::uniform_real_distribution<double> dist(0.0, 1.0);
    // dist(rng);

    return ((double)rand() / (RAND_MAX));
}

double
interiorMass(double cloudMass, double a_s, double r) {
    return cloudMass * pow(r, 3.0) / pow(pow(r, 2.0) + pow(a_s, 2.0), 1.5);
}

double
calcMag(vector<double>& vec) {
    double mag = 0;
    for (double val : vec) { mag += pow(val, 2); }
    return sqrt(mag);
}

vector<double> 
crossProduct(const vector<double>& a, const vector<double>& b)
{
    return { (b[2] * a[1]) - (a[2] * b[1]),
             (b[0] * a[2]) - (a[0] * b[2]),
             (b[1] * a[0]) - (a[1] * b[0]) };
}

double
dotProduct(const vector<double>& a, const vector<double>& b) {
    double sum = 0;
    for (int i = 0; i < 3; i++) sum += a[i] * b[i];
    return sum;
}

double
initial_r(double r_half) {
    vector<double> pos;
    double a_s = r_half * sqrt(pow(2.0, (2.0/3.0)) - 1.0);
    double r_mass_perc = r_half * settings::massPerc * sqrt((pow(2, (2.0 / 3.0)) - 1.0) / (1.0 - pow(settings::massPerc, (2.0 / 3.0))));

    double r, r_max = r_mass_perc;

    do {
        r = a_s / sqrt(pow(RNG(), (-2.0 / 3.0)) - 1); // from Bob & Alice paper
    } while (r > r_max);
    return r;
}

vector<double>
initial_pos(int idx, double r_half) {
    double r;
    if (settings::uniform_r != -1) r = settings::uniform_r;
    else if (settings::lin_dist_r != -1) r = ((double)idx + 1) * settings::lin_dist_r / settings::N;
    else r = initial_r(r_half);

    double phi_sp;
    if (settings::uniform_phi != -1) phi_sp = settings::uniform_phi;
    else if (settings::lin_dist_phi != -1) phi_sp = (double)idx * 2 * consts::pi / settings::lin_dist_phi;
    else phi_sp = 2.0 * consts::pi * RNG();

    double theta_sp;
    if (settings::plCir)  theta_sp = consts::pi / 2.0;
    else theta_sp = acos(1 - 2 * RNG());

    return { r * cos(phi_sp) * sin(theta_sp),
             r * sin(phi_sp) * sin(theta_sp),
             r * cos(theta_sp) };
}

double
calcEscapeVel(double a_newt, double r, double cloudMass, double a_s) {
    double vel_esc;
    if (settings::MOND) {
        double y = a_newt/consts::a_mond;
        double nu = LelliNu(y);
        double a = a_newt * nu;

        // Circular Speed
        double v_circ = sqrt(a * r);
        // Escape Velocity
        vel_esc = v_circ * sqrt(2.0 * (pow(r, 2.0) + pow(a_s, 2.0)) / pow(r, 2.0)); 
        return vel_esc;
    } else {
        // Newtonian Escape Velocity
        vel_esc = sqrt(2.0 * consts::G * cloudMass) *
            pow((pow(a_s, 2) + pow(r, 2)), -0.25);
    }
    return vel_esc;
}

double
distribution(double ratio_q) {
    return pow(ratio_q, 2) * pow((1.0 - pow(ratio_q, 2)), 3.5);
}

int
sign(double num) {
    if (num >= 0) return 1;
    return -1;
}

double
calcVelocityMag(double cloudMass, double a_s, double r, double a_newt) {
    // Rejection technique to determine vel magnitude as portion of esc vel
    double ratio_q = 0.0;  // sets up esc vel and distr func for Newtonian case
    double prob_g = 0.1;
    while (prob_g > distribution(ratio_q)) {
        ratio_q = RNG();
        prob_g = RNG() * 0.1;
    }
    double vel_esc = calcEscapeVel(a_newt, r, cloudMass, a_s);
    return ratio_q * vel_esc;
}

vector<double> 
initial_vel(double cloudMass, double r_half, vector<double> pos) {
    vector<double> vel;
    double a_s = r_half * sqrt(pow(2.0, (2.0/3.0)) - 1.0);
    double r = calcMag(pos);
    double mass_int = interiorMass(cloudMass, a_s, r);
    double a_newt = consts::G * mass_int / pow(r, 2.0);

    double vel_mag, v_theta, v_phi;
    if (settings::plCir) {
        vel_mag = sqrt(a_newt * r);
        vector<double> unitPos = { pos[0] / r, pos[1] / r, pos[2] / r };
        vector<double> v_hat;
        v_hat = crossProduct(unitPos, { 0, 0, 1 });
        if (pos[1] >= 0)
            v_theta = acos(v_hat[2]);
        else
            v_theta = acos(-v_hat[2]) + consts::pi;
        v_phi = atan(v_hat[1] / v_hat[0]);
    } else {
        vel_mag = calcVelocityMag(cloudMass, a_s, r, a_newt);
        // Randomize Angles
        v_theta = acos(1 - 2 * RNG());
        v_phi = 2 * consts::pi * RNG();
    }
    
    return {vel_mag * sin(v_theta) * cos(v_phi),
            vel_mag * sin(v_theta) * sin(v_phi),
            vel_mag * cos(v_theta)};
}

// Constructor
Star::Star(int idx, double r_half, double cloudMass, double hostR, double hostM) {
    this->ID = idx;
    this->time = 0;
    this->timestep = settings::initStep;
    this->positionMatrix[0] = initial_pos(idx, r_half);
    this->positionMatrix[1] = initial_vel(cloudMass, r_half, positionMatrix[0]);
    this->a_and_adot(r_half, cloudMass, hostR, hostM, false);
}

Star::Star(int idx, double r_half, double cloudMass, double hostR, double hostM, 
            vector<double> initR, vector<double> initV) {
    this->ID = idx;
    this->time = 0;
    this->timestep = settings::initStep;
    this->positionMatrix[0] = initR;
    this->positionMatrix[1] = initV;
    this->a_and_adot(r_half, cloudMass, hostR, hostM, false);
}

// Getters
// phi = aziumthal, theta = polar
vector<double>
Star::getPosSph() {
    double x = positionMatrix[0][0];
    double y = positionMatrix[0][1];
    double z = positionMatrix[0][2];
    double rho = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    double phi = (x >= 0 ? atan(y / x) : atan2(y, x));
    double theta = (z >= 0 ? atan(sqrt(pow(x, 2) + pow(y, 2)) / z) :
                             atan2(sqrt(pow(x, 2) + pow(y, 2)), z));
    return { rho, phi, theta };
}

vector<double>
Star::getVelSph() {
    vector<double> pos_sph = getPosSph();
    double rho = pos_sph[0], phi = pos_sph[1], theta = pos_sph[2];
    double vx = positionMatrix[1][0];
    double vy = positionMatrix[1][1];
    double vz = positionMatrix[1][2];
    double v_rho = (vx * sin(theta) * cos(phi)) + (vy * sin(theta) * sin(phi)) + (vz * cos(theta));
    double v_phi = (-vx * sin(phi)) + (vy * cos(phi));
    double v_theta = (vx * cos(theta) * cos(phi)) + (vy * cos(theta) * sin(phi)) - (vz * sin(theta));
    return {v_rho, v_phi, v_theta};
}

bool
Star::isBound(double G_eff, double host_M, double a_s) {
    double r = calcMag(positionMatrix[0]);
    double v = calcMag(positionMatrix[1]);
    double v_esc = sqrt(2 * G_eff * host_M) / pow(pow(r, 2) + pow(a_s, 2), 0.25); // from Bob & Alice with G_eff
    return v < v_esc;
}

// Setters
// Mutators

void 
Star::predict(double tstep) {
    double tstep2 = pow(tstep, 2);
    double tstep3 = pow(tstep, 3);
    vector<double> r_pred(3, 0), v_pred(3, 0);
    for (int i = 0; i < 3; i++) {
        r_pred[i] = positionMatrix[0][i] + positionMatrix[1][i]*tstep + 
                    positionMatrix[2][i]*tstep2/2.0 + 
                    positionMatrix[3][i]*tstep3/6.0;
        v_pred[i] = positionMatrix[1][i] + positionMatrix[2][i]*tstep + 
                    positionMatrix[3][i]*tstep2/2.0;
    }
    predMat[0] = r_pred;
    predMat[1] = v_pred;
}

// Newtonian external field
vector<vector<double>>
externalField(vector<double>& pos, vector<double>& vel,
    double hostR, double hostM) {
    //External acceleration and jerk
    vector<double> a_e(3, 0), j_e(3, 0);

    // position relative to the host galaxy at <hostR, 0, 0>
    vector<double> pos_e = { pos[0] - hostR, pos[1], pos[2] };
    double r_e = calcMag(pos_e);
    double r_e2 = pow(r_e, 2);
    double r_e3 = pow(r_e, 3);
    double factor1 = consts::G * hostM / r_e3;
    double re_dot_v = dotProduct(pos_e, vel);
    double factor2 = factor1 * 3.0 * re_dot_v / r_e2;

    /* components of gnedot EFE - 1; note here we can use velocity of star, vel, since the host
   is not moving. If the host moves, this will have to be the relative velocity. */
    for (int i = 0; i < 3; i++) {
        a_e[i] = -pos_e[i] * factor1;
        j_e[i] = (-vel[i] * factor1) + (pos_e[i] * factor2);
    }
    return { a_e, j_e };
}

// Returns EFE corrected internal acceleration and jerk
vector<vector<double>>
EFE(vector<double>& pos, vector<double>& vel, vector<vector<double>>& a_and_j, double hostR, double hostM) {
    // Newtonian external a  and a-dot
    vector<double> gne = a_and_j[1];
    vector<double> gne_dot = a_and_j[4];
    // Newtonian total a and a-dot (int + ext)
    vector<double> gntot = a_and_j[2];
    vector<double> gntot_dot = a_and_j[5];
    // Newtonian CoM external a and a-dot
    vector<double> origin(3, 0);
    vector<vector<double>> ext_fieldCoM = externalField(origin, origin, hostR, hostM);
    vector<double> gnCoM = ext_fieldCoM[0];
    vector<double> gnCoM_dot = ext_fieldCoM[1];

    double gne_mag = calcMag(gne);
    double gntot_mag = calcMag(gntot);
    double gnCoM_mag = calcMag(gnCoM);

    // y_e and nu_e for external field
    double y_e = gne_mag / consts::a_mond;
    double nu_e = LelliNu(y_e);
    double nu_e_dot = LelliNu_Dot(y_e, nu_e, gne_mag, dotProduct(gne, gne_dot));
    // y and nu for total field
    double y = gntot_mag / consts::a_mond;
    double nu = LelliNu(y);
    double nu_dot = LelliNu_Dot(y, nu, gntot_mag, dotProduct(gntot, gntot_dot));
    // y_c and nu_c for CoM
    double y_c = gnCoM_mag / consts::a_mond;
    double nu_c = LelliNu(y_c);
    double nu_c_dot = LelliNu_Dot(y_c, nu_c, gnCoM_mag, dotProduct(gnCoM, gnCoM_dot));

    // EFE MOND Corrected internal acc and jerk
    vector<double> gi(3), gi_dot(3);
    for (int i = 0; i < 3; i++) {
        gi[i] = nu * gntot[i] - nu_c * gnCoM[i];
        gi_dot[i] = nu * gntot_dot[i] + nu_dot * gntot[i]
            - nu_c * gnCoM_dot[i] - nu_c_dot * gnCoM[i];
    }
    return { gi, gi_dot };
}

PosMat
MONDCorrections(vector<double>& pos, vector<double>& vel, vector<vector<double>>& a_and_j, double hostR, double hostM) {
    // Switch notation here - use g's instead of a's for accelerations to
    // make it consistent with the MOND notes
    vector<double> gni = a_and_j[0];
    vector<double> gni_dot = a_and_j[3];

    PosMat corrections;
    // Isolated MOND case
    if (!settings::EFE) {
        double gni_mag = calcMag(gni);
        double y = gni_mag / consts::a_mond;
        double nu = LelliNu(y);
        double nu_dot = LelliNu_Dot(y, nu, gni_mag, dotProduct(gni, gni_dot));

        // MONDian Corrected acc and jerk
        vector<double> g(3), g_dot(3);
        for (int i = 0; i < 3; i++) {
            g[i] = gni[i] * nu;
            g_dot[i] = nu * gni_dot[i] + nu_dot * gni[i];
        }

        corrections.push_back(g);
        corrections.push_back(g_dot);
    } else { 
    // EFE case
        vector<vector<double>> g_EFE = EFE(pos, vel, a_and_j, hostR, hostM);

        corrections.push_back(g_EFE[0]);
        corrections.push_back(g_EFE[1]);
    }   
    return corrections;
}

PosMat
STVG(vector<double>& pos, vector<double>& vel, vector<double>& acc,
    vector<double>& jerk, double M, double M_dot) {
    // empirical factors from Moffat & Toth 2023
    double alpha_inf = 19;
    double E = 2.5 * pow(10, 4);
    double D = 6250;

    double r = calcMag(pos);
    double v = calcMag(vel);
    double alpha = alpha_inf * M / pow(sqrt(M) + E, 2);
    double r_0 = sqrt(M) / D;
    double acc_factor = 1 + alpha - alpha * (1 + (r / r_0)) * exp(-r / r_0);

    double A = 1 / pow(sqrt(M) + E, 2);
    double B = 1 + (r * D / sqrt(M));
    double C = exp(-r * D / sqrt(M));

    double term1 = (-1 / sqrt(M)) + (1 / (A * (sqrt(M) + E))) + (D * ((1 / B) - 1) * ((r / (2 * M)) - (v / M_dot)));
    double term2 = A - (sqrt(M) / sqrt(A)) + sqrt(M) * A * B * C * term1;
    double jerk_factor = alpha_inf * M_dot * term2;

    for (int i = 0; i < 3; i++) {
        acc[i] *= acc_factor;
        jerk[i] = (acc_factor * jerk[i]) + (acc[i] * jerk_factor);
    }
    
    return  { acc, jerk };
}

// bool prediction: 0 = positionMatrix, 1 = predictionMatrix
void
Star::a_and_adot(double r_half, double cloudMass, double hostR, double hostM, bool prediction) {
    // the integrated acceleration and jerk
    vector<double> a(3, 0);
    vector<double> j(3, 0);

    vector<double> pos;
    vector<double> vel;
    if (prediction) {
        pos = predMat[0];
        vel = predMat[1];
    } else {
        pos = positionMatrix[0];
        vel = positionMatrix[1];
    }

    // Distance from center, r, powers of r, & r dot v:
    double r = calcMag(pos);
    double r2 = pow(r, 2);
    double r3 = r2 * r;
    double r5 = r2 * r3;
    double rdotv = dotProduct(pos, vel);

    // Scalelength and internal mass functions
    double a_s = r_half * sqrt(pow(2.0, (2.0/3.0)) - 1.0);
    double int_mass = interiorMass(cloudMass, a_s, r);

    // If there's a centeral blackhole, add to int_mass
    if (settings::blackHole) {
        int_mass += settings::mBlack;
    }

    // Common factors needed for a and a-dot
    double factor_1 = int_mass / r3;
    double factor_2 = int_mass * 3.0 * rdotv / r5;
    double int_mass_dot = (3.0 * r * cloudMass * rdotv * pow(a_s, 2)) / 
                                        pow((r2 + pow(a_s, 2)), (2.5));
    double factor_3 = int_mass_dot / r3;

    // Newtonian internal, external, and total acceleration and jerk
    vector<double> a_i(3,0), a_e(3,0), a_t(3,0);
    vector<double> j_i(3,0), j_e(3,0), j_t(3,0);
    for (int i = 0; i < 3; i++) {
        a_i[i] = -(consts::G)*(pos[i] * factor_1);
        j_i[i] = (consts::G)*(-vel[i] * factor_1 + pos[i] * factor_2 - pos[i] * factor_3);
    }
    if (settings::EFE) {
        // External field on star
        vector<vector<double>> ext_field = externalField(pos, vel, hostR, hostM);
        a_e = ext_field[0];
        j_e = ext_field[1];

        // External field on CoM (origin for now)
        vector<double> origin(3, 0), a_CoM, j_CoM;
        vector<vector<double>> ext_fieldCoM = externalField(origin, origin, hostR, hostM);
        a_CoM = ext_fieldCoM[0];
        j_CoM = ext_fieldCoM[1];

        for (int i = 0; i < 3; i++) {
            a_t[i] = a_i[i] + a_e[i];
            j_t[i] = j_i[i] + j_e[i];
            a[i] = a_t[i] - a_CoM[i];
            j[i] = j_t[i] - j_CoM[i];
        }
    } else {
        a = a_i;
        j = j_i;
    }

    // Set of int, ext, and tot accelerations and jerks for easier transport
    vector<vector<double>> a_and_j = { a_i, a_e, a_t, j_i, j_e, j_t };
    if (settings::MOND) {
        PosMat corrections = MONDCorrections(pos, vel, a_and_j, hostR, hostM);
        a = corrections[0];
        j = corrections[1];
    } else if (settings::STVG) {
        PosMat corrections = STVG(pos, vel, a, j, int_mass, int_mass_dot);
        a = corrections[0];
        j = corrections[1];
    }

    if (prediction) {
       predMat[2] = a;
       predMat[3] = j;
    } else {
       positionMatrix[2] = a;
       positionMatrix[3] = j;
    }
}

void
Star::correct(double tstep) {
    double tstep2 = pow(tstep, 2);
    double tstep3 = tstep2 * tstep;
    double tstep4 = pow(tstep2, 2);
    double tstep5 = tstep2 * tstep3;
    positionMatrix[4] = {0, 0, 0};
    positionMatrix[5] = {0, 0, 0};

    for (int i = 0; i < 3; i++) {
        positionMatrix[4][i] = (-6.0 * (positionMatrix[2][i] - predMat[2][i]) + 
                    tstep * (4.0 * positionMatrix[3][i] + 2.0 * predMat[3][i])) 
                    / tstep2;
        positionMatrix[5][i] = (12.0 * (positionMatrix[2][i] - predMat[2][i]) + 
                            6.0 * tstep * (positionMatrix[3][i] + predMat[3][i]))
                                / tstep3;

        positionMatrix[0][i] = predMat[0][i] + tstep4*positionMatrix[4][i]/24 + 
                                               tstep5*positionMatrix[5][i]/120;
        positionMatrix[1][i] = predMat[1][i] + tstep3*positionMatrix[4][i]/6 + 
                                               tstep4*positionMatrix[5][i]/24;
        positionMatrix[2][i] = predMat[2][i];
    }
    predMat = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, 
            {0, 0, 0}, {0, 0, 0} };
}

double
newTimestep(PosMat& predMat, vector<double>& a3dotVec) {
    double newStep;
    
    double acc = calcMag(predMat[2]);
    double jerk = calcMag(predMat[3]);
    double a2dot = calcMag(predMat[4]);
    double a3dot = calcMag(a3dotVec);

    newStep = sqrt(consts::eta * (acc * a3dot + pow(jerk, 2)) / 
                (jerk * a3dot + pow(a2dot, 2)) );
    return newStep;
}

void
Star::incrementTime() {
    // Update time after integrating
    time += timestep;
    // Calc prediction for a2dot
    vector<double> a2dot_pred(3, 0);
    for (int i = 0; i < 3; i++)
        a2dot_pred[i] = positionMatrix[4][i] + timestep * positionMatrix[4][i];
    predMat[4] = a2dot_pred;
    double newStep = newTimestep(predMat, positionMatrix[5]);
    // Already being done in correct
    /* for (int j = 0; j < 4; j++) {
        positionMatrix[j] = predMat[j];
        predMat[j] = {0, 0, 0};
    } */
    // The new step can't be larger than 1.2*(the old step)
    newStep = min(1.2 * timestep, newStep);
    // Check to see if timestep less than minimum step
    timestep = max(newStep, settings::minstep);
}

double
Star::updateStar(double rhalf, double gMass, double hostR, double hostM) {
    // Integrate star forward one step
    this->predict(timestep);
    this->a_and_adot(rhalf, gMass, hostR, hostM, true);
    this->correct(timestep);
    // Update time to after timestep, give update star new timestep
    this->incrementTime();
    return this->getAfterTimestep();
}

void
Star::addLog() {
    vector<double> logEntry = {time};
    for (int i = 0; i < 3; i++)
        logEntry.push_back(positionMatrix[0][i]);
    /* for (int i = 0; i < 3; i++)
        logEntry.push_back(positionMatrix[1][i]); */
    log.push_back(logEntry);
}
