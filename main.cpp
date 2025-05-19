// Copyright 2025 Will Snider

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <utility>
#include <vector>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <math.h>
#include <omp.h>
#include <cmath>
#include <chrono>
#include <random>
#include "Star.h"
#include "Galaxy.h"
#include "settings.h"

using namespace std;
using namespace std::string_literals;
using PosMat = vector<vector<double>>;
using Model = unordered_map<string, double>;
using PairList = vector<vector<double>>;

class Timer
{
public:
    Timer() : beg_(clock_::now()) {}
    void reset() { beg_ = clock_::now(); }
    double elapsed() const {
        return std::chrono::duration_cast<second_>
            (clock_::now() - beg_).count();
    }

private:
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > second_;
    std::chrono::time_point<clock_> beg_;
};

void
getModelStats(Model& stats) {
    ifstream txtFile("modelData.txt");
    string skim, line;

    // Skip file lines until specified model is reached
    int lineNum = 0;
    while (getline(txtFile, skim)) {
        lineNum++;
        if (skim == settings::modelName) break;
    }
    // Hash all given parameters with their associated values
    while (getline(txtFile, line) && line != "-") {
        lineNum++;
        istringstream words(line);
        string name, eqs, value;
        words >> name >> eqs >> value;
        //cout << "Name: " << name << " eqs: " << eqs << " Value: " 
        //     << value << endl;
        if (value != "")
            stats[name] = stod(value);
    }
    if (stats.find("M_L") == stats.end() || stats["M_L"] == 0)
        cout << "Invalid Model, M/L ratio error" << endl;
    if (stats.find("lum") == stats.end() || stats["lum"] == 0)
        cout << "Invalid Model, luminosity error" << endl;
    if (stats.find("r_half") == stats.end() || stats["r_half"] == 0)
        cout << "Invalid Model, r_half error" << endl;

    stats["mass"] = stats["lum"] * stats["M_L"];
    // Tidal Radius Eq  (Banik 2022)
    stats["r_tidal"] = pow(2.0, 1.0 / 6.0) * stats["R_mw"] * pow(stats["mass"] / stats["M_mw"], 1.0 / 3.0) / 3.0;
}

void setModelStats(Model& stats, int itr) {
    double log_int_rat;
    vector<double> log_gi_a0 = { -2, -1, 0, 1 }; 
    vector<double> log_ge_a0 = { -2, -1, 0, 1 };

    if (settings::runs == 1) {
        cout << "Enter newtonian log(gi/a0) at r_half: ";
        cin >> log_int_rat;
    } else {
        cout << "Using log(gi/a0) = " << log_gi_a0[itr % log_gi_a0.size()] << endl;
        log_int_rat = log_gi_a0[itr % log_gi_a0.size()];
    }
    stats["log_int_rat"] = log_int_rat;
    cout << endl;

    stats["mass"] = settings::toy_mass;
    stats["r_half"] = sqrt(consts::G * settings::toy_mass / (2.0 * consts::a_mond * pow(10, log_int_rat)));

    if (settings::extField) {
        double log_ext_rat;

        if (settings::runs == 1) {
            cout << "Enter newtonian log(ge/a0) at CoM: ";
            cin >> log_ext_rat;
        } else {
            cout << "Using log(ge/a0) = " << log_ge_a0[itr / log_ge_a0.size()] << endl;
            log_ext_rat = log_ge_a0[itr / log_ge_a0.size()];
        }
        stats["log_ext_rat"] = log_ext_rat;
        cout << endl;

        stats["R_mw"] = settings::toy_hostR;
        // double ext_rat_N = pow(ext_rat, 2) / (1 + ext_rat);  -- used if ext_rat is dynamic/MOND ratio
        stats["M_mw"] = (pow(10, log_ext_rat) * consts::a_mond) * pow(settings::toy_hostR, 2) / consts::G;

        // Tidal Radius Eq  (Banik 2022)
        stats["r_tidal"] = pow(2.0, 1.0 / 6.0) * stats["R_mw"] * pow(stats["mass"] / stats["M_mw"], 1.0 / 3.0) / 3.0;
    }
}

void
printInitConds(Model& modelStats, double Tmax) {
    cout << "\nRunning with the following initial conditions..." << endl;
    if (!settings::g_ratios) cout << "Model Name: " << settings::modelName << endl;
    else {
        cout << "Acc Ratios:" << endl;
        cout << "\tgi/a0 = " << pow(10, modelStats["log_int_rat"]) << endl;
        cout << "\tge/a0 = " << pow(10, modelStats["log_ext_rat"]) << endl;
    }
    cout << "Dsph Mass: \t" << modelStats["mass"] << " solar masses" << endl;
    cout << "Pop. Size: \t" << settings::N << endl;
    cout << "R_Half: \t" << modelStats["r_half"] << " pc" << endl;
    cout << "Mass % Limit: \t" << settings::massPerc << endl;
    cout << "Runtime: \t" << (Tmax) / 1000 << " billion years" << endl;
    cout << "Switches:" << endl;
    cout << "\tPrinting: ";
        settings::pos_out ? cout << "\tpos" : cout << " ";
        settings::vel_out ? cout << "\tvel" : cout << " ";
        settings::acc_out ? cout << "\tacc" << endl : cout << " " << endl;
    if (settings::plCir) cout << "\tPlanar Circular Orbits: \tON" << endl;
    if (settings::uniform_r != -1) cout << "\tUniform R: \t" << settings::uniform_r << " pc" << endl;
    if (settings::uniform_phi != -1) cout << "\tUniform Phi: \t" << settings::uniform_phi << " rad" << endl;
    if (settings::lin_dist_r != -1) cout << "\tLinearly Distributed R: \t" << settings::lin_dist_r << " pc (max r)" << endl;
    if (settings::lin_dist_phi != -1) cout << "\tLinearly Distributed Phi: \t" << settings::lin_dist_phi << " groups" << endl;
    if (settings::CenterOfMass) cout << "\tCoM Tracking: \tON" << endl;
    if (settings::trackTidalR) cout << "\tTidal R Tracking: \tON" << endl;
    if (settings::trackSkews) cout << "\tSkew Tracking: \tON" << endl;
    if (settings::blackHole) cout << "\tBlackhole: \tON" << endl;
    if (settings::blackHole) cout << "\t\t\tmass: " << settings::mBlack << endl;
    if (settings::MOND) cout << "\tMOND: \t\tON" << endl;
    if (settings::extField) {
        cout << "\tExt Field: \t\tON" << endl;
        cout << "\t\tHost Mass: \t" << modelStats["M_mw"] << " solar masses" << endl;
        cout << "\t\tTidal Radius: \t" << modelStats["r_tidal"] << " pc" << endl;
    }
    cout << "Output Format: ";
    switch (settings::format) {
        case 0: 
            cout << "\tAnimation [time N / ID x y z...]" << endl;
            break;
        case 1:
            cout << "\tStatistics [ID x y z...]" << endl;
            break;
        case 2:
            cout << "\tAnimation w/ COM [time N COM / ID x y z...]" << endl;
            break;
    }
}

void
recordCOM(double time, Galaxy& g, vector<vector<double>>& out) {
    g.calcCOM();
    PosMat COM = g.getCOM();
    vector<double> record(7, 0);
    record[0] = time;
    for (int i = 1; i < 7; i++) record[i] = COM[i / 3][i % 3];
    out.push_back(record);
}

void
recordTides(double time, Galaxy& g, vector<vector<double>>& out) {
    if (time == 0) out.push_back({ g.getRTidal() });

    double r_tidal_exp = 0;
    vector<Star>& pop = g.getPopulation();
    for (int i = 0; i < settings::N; i++) {
        if (pop[i].isBound(g.getGeff(), g.getHostMass(), g.getRHalf() * sqrt(pow(2.0, (2.0 / 3.0)) - 1.0))) {
            vector<double> pos = pop[i].getPos();
            double r = calcMag(pos);
            if (r > r_tidal_exp) r_tidal_exp = r;
        }
    }
    out.push_back({ time, r_tidal_exp });
}

void
output(Galaxy& g, ofstream& outFile, double time) {
    // Header of output block, dependent on format setting
    if (settings::format != 1) {
        outFile << time << "\t" << settings::N << "\t" << g.getRHalf();
        if (settings::format == 2) {
            PosMat COM = g.getCOM();
            outFile << "\t" << COM[0][0] << " " << COM[0][1] << " " << COM[0][2];
        } else if (settings::format == 4) {
            PosMat COM = g.getCOM();
            outFile << "\t" << COM[1][0] << " " << COM[1][1] << " " << COM[1][2];
        }
        outFile << endl;
    }

    if (!g.isUniformTime()) // must integrate all stars to same time for output
        g.wrangleStars(time);

    vector<Star> pop = g.getPopulation();
    for (int i = 0; i < settings::N; i++) {
        outFile << pop[i].getID() << "\t";
        if (settings::pos_out) {
            vector<double> pos = pop[i].getPos();
            outFile << pos[0] << "\t" << pos[1] << "\t" << pos[2] << "\t";
        } 
        if (settings::vel_out) {
            vector<double> vel = pop[i].getVel();
            outFile << vel[0] << "\t" << vel[1] << "\t" << vel[2] << "\t";
        }
        if (settings::acc_out) {
            vector<double> acc = pop[i].getAcc();
            outFile << acc[0] << "\t" << acc[1] << "\t" << acc[2] << "\t";
        }
        if (pop[i].isBound(g.getGeff(), g.getMass(), g.getRHalf() * sqrt(pow(2.0, (2.0 / 3.0)) - 1.0)))
            outFile << "b" << "\t";
        else
            outFile << "u" << "\t";
        if (pop[i].isFrozen())
            outFile << "f" << "\t";
        else
            outFile << "m" << "\t";
        outFile << endl;
    }
}

Star&
minTimeStar(Galaxy& g) {
    vector<Star>& pop = g.getPopulation();
    double timeMin = pop[0].getAfterTimestep();
    size_t minIdx = 0;
    for (size_t i = 1; i < pop.size(); i++) {
        double t = pop[i].getAfterTimestep();
        if (settings::freezeStrays && pop[i].isFrozen()) continue;
        if (t < timeMin) {
            timeMin = t;
            minIdx = i;
        }
    }
    return g.getStar(minIdx);
}

PairList
makeProjection(ifstream& read) {
    
    string line, Time, skip, ID;
    string X, Y, Z, VX, VY, VZ; // , ax, ay, az;
    double x, y, z, vx, vy, vz;
    double r, perp_r, v_rad;

    PairList projection = {};
    getline(read, line);
    istringstream head(line);
    head >> Time >> skip;
        // time = stod(Time);
        for (int s = 0; s < settings::N; s++) {
            getline(read, line);
            istringstream vals(line);
            vals >> ID >> X >> Y >> Z >> VX >> VY >> VZ;
            // if (settings::acc_out) vals >> ax >> ay >> az;
            x = stod(X); y = stod(Y); z = stod(Z);
            vx = stod(VX); vy = stod(VY); vz = stod(VZ);
            r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
            if (settings::trunc_dist == -1 || r <= settings::trunc_dist) {
                if (settings::axis == 1) {
                    perp_r = sqrt(pow(y, 2) + pow(z, 2));
                    v_rad = vx;
                }
                else if (settings::axis == 2) {
                    perp_r = sqrt(pow(x, 2) + pow(z, 2));
                    v_rad = vy;
                }
                else {
                    perp_r = sqrt(pow(x, 2) + pow(y, 2));
                    v_rad = vz;
                }
                projection.push_back({ perp_r, v_rad });
            }
        }
    return projection;
}

PairList
binDispersion(PairList& projection) {
    PairList profile = {};
    for (double i = 0; i < settings::bins; i++) {
        double r_sum = 0, v_sum = 0, v2_sum = 0;
        for (double star = 0; star < settings::bins; star++) {
            double id = settings::bins * i + star;
            r_sum += projection[id][0];
            v_sum += projection[id][1];
            v2_sum += pow(projection[id][1], 2);
        }
        double r_avg = r_sum / settings::bins;
        double v_avg = v_sum / settings::bins;
        double v2_avg = v2_sum / settings::bins;
        double sigma = sqrt(v2_avg - pow(v_avg, 2));
        profile.push_back({ r_avg, sigma });
    }
    return profile;
}

void
outputDispProfile(const PairList& profile, int itr) {
    string dispFileName;
    if (settings::runs != 1) dispFileName = "Run_" + to_string(itr+1) + settings::dispOutput;
    else dispFileName = settings::skewOutput;
    ofstream out(dispFileName);
    for (int i = 0; i < settings::bins; i++) {
        out << profile[i][0] << "\t" << profile[i][1] << endl;
    }
}

double
calcBulkDispersion(PairList& profile) {
    double bulk = 0;
    for (int b = 0; b < settings::bins; b++) {
        bulk += profile[b][1];
    }
    return (bulk / settings::bins);
}

double
dispersion(int outputCount, int itr) {
    cout << "Running Dispersion Calculation" << endl;
    cout << "Number of data records = " << outputCount << endl;

    ifstream read(settings::simOutput);
    PairList projection, snapDisp;
    vector<PairList> dispTimeline = {};
    for (int i = 0; i < outputCount; i++) {
        projection = makeProjection(read);
        sort(projection.begin(), projection.end(),
            [](const std::vector<double>& a, const std::vector<double>& b) {
                return a[0] < b[0];
            });
        snapDisp = binDispersion(projection);
        dispTimeline.push_back(snapDisp);
    }
    read.close();
    cout << "Dispersion timeline created" << endl;

    PairList profile = {};
    for (int b = 0; b < settings::bins; b++) {
        double r_sum = 0, sig_sum = 0;
        for (int t = 0; t < outputCount; t++) {
            r_sum += dispTimeline[t][b][0];
            sig_sum += dispTimeline[t][b][1];
        }
        // Take average and convert disp units to km/s
        double r_avg = r_sum / outputCount;
        double sig_avg = sig_sum * 0.978 / outputCount;
        profile.push_back({ r_avg, sig_avg });
        //dispOut << r_avg << ", " << sig_avg
    }
    outputDispProfile(profile, itr);
    double bulk = calcBulkDispersion(profile);
    return bulk;
}

// Submethod for Haghi dipsersion prediction
double 
F(double a_e, double x) {
    double y = log10(a_e / consts::a_mond);
    double A = 5.3 / (10.56 + pow((y + 2), 3.22));
    double B = pow(10, -(1.65 * y + 0.0065));
    double C = 3.788 * y + 0.006;

    return (-A / 4.0) * (log(exp(-x / A) + B) + C);
}

// Based on Haghi et. al. 2019
double
dispPredHaghi(double M, double r_half, double hostR, double hostM) {
    double x = log10((consts::G * M / (2 * pow(r_half, 2))) / consts::a_mond);
    double sigmaM = pow((4.0 / 81.0) * consts::G * M * consts::a_mond, 0.25) * pow(1 + 0.56 * exp(3.02*x),0.184);
    if (!settings::extField) {
        return sigmaM;
    }  else {
        double a_e = (consts::G * hostM) / pow(hostR, 2);
        double sigmaEF = sqrt((consts::G * M * consts::a_mond) / (4 * hostR * a_e)); // pow(10, (log10(sigmaM) + F(a_e, x)));
        return sigmaEF;
    }
}

void
dataDump(int itr, Galaxy& dsph, int outputCount, PairList& skews, 
         PairList& tideOutput, vector<vector<double>>& COMrecord) {
    if (settings::trackSkews) {
        string skewFileName;
        if (settings::runs != 1) skewFileName = "Run_" + to_string(itr+1) + settings::skewOutput;
		else skewFileName = settings::skewOutput;
        ofstream skewFile(skewFileName);
        for (int i = 0; i < skews.size(); i++) skewFile << skews[i][0] << "\t" << skews[i][1] << endl;
        skewFile.close();
    }
    if (settings::trackTidalR) {
        string tideFileName;
        if (settings::runs != 1) tideFileName = "Run_" + to_string(itr+1) + settings::tidalOutput;
        else tideFileName = settings::tidalOutput;
        ofstream tidal(tideFileName);
        tidal << tideOutput[0][0] << endl;
        for (int i = 1; i < tideOutput.size(); i++) tidal << tideOutput[i][0] << "\t" << tideOutput[i][1] << endl;
    }
    if (settings::CenterOfMass) {
        string COMFileName;
        if (settings::runs != 1) COMFileName = "Run_" + to_string(itr+1) + settings::COM_Output;
        else COMFileName = settings::COM_Output;
        ofstream COMout(COMFileName);
        for (int i = 0; i < COMrecord.size(); i++) {
            for (int j = 0; j < 7; j++) COMout << COMrecord[i][j] << "\t";
            COMout << endl;
        }
    }

    if (settings::run_dispersion) {
        if (settings::pos_out && settings::vel_out) {
            double bulkDisp = dispersion(outputCount, itr);
            cout << "Bulk Dispersion: " << bulkDisp << endl;
            if (settings::MOND) {
                double HaghiPred = dispPredHaghi(dsph.getMass(), dsph.getRHalf(), dsph.getHostDist(), dsph.getHostMass());
                cout << "Haghi Disp. Prediction: " << HaghiPred << endl;
                if (!settings::extField) {
                    double sigmaIsoM = pow((4 * consts::G * dsph.getMass() * consts::a_mond) / 81, 0.25); // (Milgrom 1994; McGaugh& Milgrom 2013)
                    cout << "Milgram Isolated MOND disp prediction: " << sigmaIsoM << endl;
                }
            }
        }
        else
            cout << "Unable to calculate dispersion profile without recorded positions and velocities" << endl;
    }
}

// Reassign population with a controllable states for testing
/*
    vector<Star> testPop = {};
    for (int id = 0; id < N; id++) {
        testPop.push_back(test);
    }
    vector<double> initR = { 500, 0, 0 };
    vector<double> initV = { 0, -6.70711, 0 };  // -6.70711
    Star test(0, N, modelStats["r_half"], modelStats["lum"]*modelStats["M_L"], initR, initV, false);
    dsph.setPopulation({ test });
*/

int main(int argc, char* argv[]) {
    for (int itr = 0; itr < settings::runs; itr++) {
        // Initialize galaxy with given parameters
        Model modelStats;
        if (!settings::g_ratios)
            getModelStats(modelStats);
        else
            setModelStats(modelStats, itr);
        //cout << "mass: " << modelStats["M_mw"] << " radius: " << modelStats["R_mw"];
        
        Galaxy dsph(modelStats);
        cout << "\nSuccessfully initialized galaxy!\n" << endl;

        // Determine runtime
        double Tmax;
        if (settings::constRuntime) Tmax = settings::TmaxConst;
        else Tmax = settings::crossings * dsph.getTcross();
        printInitConds(modelStats, Tmax);

        // Set up output file and data storage
        string outFileName;
        if (settings::runs != 1) outFileName = "Run_" + to_string(itr + 1) + settings::simOutput;
        else outFileName = settings::simOutput;
        ofstream outFile(outFileName);
        PairList tideOutput = {};
        PairList skews = {};
        vector<vector<double>> COMrecord = {};
        
        // Realtime timer and simulation time
        Timer timer;
        double time = 0;
        output(dsph, outFile, time);  // Output initial positions
        if (settings::trackTidalR) recordTides(time, dsph, tideOutput);
        // cout << "Initial anisotropy coefficient: " << dsph.calcAnisotropyFactor() << endl;
        int outputCount = 1; // tracks outputs for calculating output times
        // Begin integration loop
        while (time < Tmax) {
            // Find star with smallest post-timestep time
            Star& minStar = minTimeStar(dsph);
            // Integrate star to next time step, returns time after timestep
            time = minStar.updateStar(dsph.getRHalf(), dsph.getMass(), dsph.getHostDist(), dsph.getHostMass(), dsph.getCOM());
            // Check if output time
            if (time >= outputCount * settings::outputTime) {
                if (settings::CenterOfMass) recordCOM(time, dsph, COMrecord);
                if (settings::trackTidalR) recordTides(time, dsph, tideOutput);
                output(dsph, outFile, time);
                if (settings::trackSkews) {
                    dsph.calcSkewness();
                    vector<double> pair = { time, dsph.getSkewness() };
                    skews.push_back(pair);
                }                
                // Check if time to write to console
                int timePerWrite = Tmax / settings::consoleWrites;
                if (outputCount % timePerWrite == 0) {
                    cout << "Star " << minStar.getID() << " caused output at t = " << time
                        << " My   (" << outputCount / timePerWrite << "/" << settings::consoleWrites << ")"
                        << "\tElasped Time : " << timer.elapsed() << " s" << endl;
                    timer.reset();
                }
                outputCount++;
                // if (outputCount % 50 == 0) cout << "Anisotropy coefficient: " << dsph.calcAnisotropyFactor() << endl;
            }
        }
        cout << "Complete! :)" << endl;
        outFile.close();

        // Output all remaining stored data
        dataDump(itr, dsph, outputCount, skews, tideOutput, COMrecord);
    }
}
/////////////////////////   End of Real Code    ////////////////////////////////
/*
    // Create file of all known locations and then create a datafile based on
    // those locations
    double populationTimeMax = 0;
    while (getMaxTime(dsph) < Tmax) {
        for (Star& s : population) {
            double step = s.getTimestep();
            s.predict(step);
            PosMat predictionMat = s.getPredMat();
            s.a_and_adot(dsph.getRHalf(), dsph.getMass(), predictionMat, MOND);
            s.setPredMat(predictionMat);
            s.correct(step);
            s.incrementTime();
        }
    }
*/

/*  (argc > 4 ? (True == argv[4]) : true),   // Optional mask flag
        (argc > 5 ? std::stoi(argv[5]) : 75),    // Optional percentMatch
        (argc > 6 ? std::stoi(argv[6]) : 32));   // Optional tolerance
*/

/*
double
lerp(double a, double b, double t) {
    return a + (b - a) * t;
}

vector<double>
linearInterpolation(double time, Star& s) {
    int counter = 0;
    double logTime = 0;
    vector<double>& currLog = s.getLog(counter);
    do {
        counter++;
        currLog = s.getLog(counter);
        logTime = currLog[0];
    } while (logTime < time && counter < static_cast<int>(currLog.size()));
    if (counter < static_cast<int>(currLog.size())) {
        cout << "Whoops invalid log index";
        return {0, 0, 0};
    }
    vector<double> prevLog = s.getLog(counter - 1);
    double t = (time - prevLog[0]) / (logTime - prevLog[0]);
    double x = lerp(prevLog[1], currLog[1], t);
    double y = lerp(prevLog[2], currLog[2], t);
    double z = lerp(prevLog[3], currLog[3], t);
    return {x, y, z};
}
*/

/*
    // Output initial positions
    for (int i = 0; i < g.getSize(); i++) {
        Star s = g.getStar(i);
        vector<double> initPos = s.getLog(0);
        outFile << s.getID() << " " << initPos[1] << " " << initPos[2] << " "
                << initPos[3] << endl;
    }
    vector<Star> pop = g.getPopulation();
    while (time < Tmax) {
        time += outputTime;
        outFile << time << "\t" << N << endl;
        for (int i = 0; i < g.getSize(); i++) {
            Star s = g.getStar(i);
            vector<double> estPos = linearInterpolation(time, s);
            outFile << s.getID() << " " << estPos[0] << " " << estPos[1]
                    << " " << estPos[2] << endl;
        }
    }
*/

// End of source code
