/**
 * Kalman Filter 
 * Skeleton code for teaching 
 * A3M33MKR
 * Czech Technical University 
 * Faculty of Electrical Engineering
 * Intelligent and Mobile Robotics group
 *
 * Authors: Zdeněk Kasl, Karel Košnar kosnar@labe.felk.cvut.cz
 *
 * this code is inspired by tutorial on Kalman Filter by Greg Czerniak 
 *
 * Licence: MIT (see LICENSE file)
 **/
/*TODO set up the convarience matrices better and re-measure the noice constants (process noise/measurenment noice)*/

#include<cmath>
#include<cassert>
#include<cstdlib>
#include<fstream>
#include<iostream>
#include <Eigen/Dense>

#include "gui/gui.h"
#include "systemSimulator/system.h"

using namespace imr;
using namespace gui;
using namespace Eigen;
using namespace std;

#define PRINT(name) Print(#name, (name))

void Print(const char *name, Matrix4f m);

Point KalmanFilter(Point measuredPosition);

IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");

Matrix4f A, B, Sigma, R, Q, E, C;
//Matrix2f Q;
//Matrix<float, 2, 4> C;
Vector4f mu, u;
float dt = 0.1;

//noice constants
double q = 100;
double r = 0.1;

//step constant
int step = 10;

void initMatrices() {
    A << 1, 0, dt, 0,
            0, 1, 0, dt,
            0, 0, 1, 0,
            0, 0, 0, 1;

    B << 0, 0, 0, 0,
            0, -0.5 * dt * dt, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, -dt;

    u << 9.8, 9.8, 9.8, 9.8;

    mu << -300, 100, 0, 0;

    Sigma << 0.1, 0, 0, 0,
            0, 0.1, 0, 0,
            0, 0, 0.1, 0,
            0, 0, 0, 0.1;

// /Noice consts.
//-------------------------------------------------------------------
    C <<    1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;

    Q <<    q, 0, 0, 0,
            0, q, 0, 0,
            0, 0, q, 0,
            0, 0, 0, q;

    R <<    r, 0, 0, 0,
            0, r, 0, 0,
            0, 0, r, 0,
            0, 0, 0, r;
//-------------------------------------------------------------------
    E <<    1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;
}


void help(char **argv) {
    std::cout << "\nUsage of the program " << argv[0] + 2 << std::endl
              << "Parameter [-h or -H] displays this message." << std::endl
              << " Parametr [-n or -N] number of simulation steps."
              << std::endl;
}

int main(int argc, char **argv) {
    int nSteps = 1000;
    char *dataFile;

    // Parse all console parameters
    for (int i = 0; i < argc; i++) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                //> HELP
                case 'H' :
                case 'h' :
                    help(argv);
                    break;
                case 'N':
                case 'n':
                    assert(i + 1 < argc);
                    assert(atoi(argv[i + 1]) > 1);
                    step = atoi(argv[i + 1]);
                    break;
                default :
                    std::cout << "Parameter \033[1;31m" << argv[i] << "\033[0m is not valid!\n"
                              << "Use parameter -h or -H for help." << std::endl;
                    break;
            }
        }
    }
    // All parameters parsed

    Gui gui;
    System system;

    //> comment line below in order to let the program
    //> continue right away
    gui.startInteractor();

    initMatrices();

    Point measurement;
    Point truth;
    Point kfPosition;


    for (int i = 1; i < nSteps; i++) {
        system.makeStep();
        truth = system.getTruthPosition();
        measurement = system.getMeasurement();
        kfPosition = KalmanFilter(measurement);
        gui.setPoints(truth, measurement, kfPosition);

        //> comment line below in order to let the program
        //> continue right away
        //if (i % step == 0) gui.startInteractor();

    }
    gui.startInteractor();
    return EXIT_SUCCESS;
}


/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//                               implement Kalman Filter here
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


Point KalmanFilter(const Point measuredPosition) {
    Point ret;
    Vector4f mu_actual, z, y;
    Matrix4f Sigma_actual, K, S, S_inv;

    z << measuredPosition.x,
            measuredPosition.y,
            0,
            0;

    PRINT(A);
    PRINT(C);
    //-------------------------------------------------------------------------------
    mu_actual = A * mu + B * u;

    //Error estimation
    //-------------------------------------------------------------------------------
    Sigma_actual = A * Sigma * A.transpose() + R;
    PRINT(Sigma_actual);

    //K gain
    //-------------------------------------------------------------------------------
    S = C * Sigma_actual * C.transpose() + Q;
    PRINT(S);

    K = Sigma_actual * C.transpose() * S.inverse();
    PRINT(K);

    //improvement
    //-------------------------------------------------------------------------------
    y = z - C * mu_actual;
    mu = mu_actual + K * y;
    Sigma = (E - K * C) * Sigma_actual;
    PRINT(Sigma);

    ret.x = mu(0);
    ret.y = mu(1);

    //Simulation of system
    //ret.x = mu_actual(0);
    //ret.y = mu_actual(1);
    //mu = mu_actual;

    cout << "X: " << ret.x << "    Y: " << ret.y <<
             "    Vx: " << mu(2) << "    Vy: " << mu(3) << endl;
    cout << "-------------------------------------------------------------------" << endl;

    return ret;
}

void Print(const char *name, Matrix4f m){
    cout << name << ": " << endl;
    cout << m.format(CleanFmt) << endl;
    cout << "-------------------------------------------------------------------" << endl;
}

/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



