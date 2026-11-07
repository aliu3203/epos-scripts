#ifndef STAR200_H
#define STAR200_H
// STAR Au+Au 200 GeV, HEPData ins2928164 (auto-extracted from the CSV tables).
// Order is most-peripheral-first as published: 70-80,60-70,50-60,40-50,30-40,20-30,10-20,5-10,0-5
const int NSTAR = 9;
const double STAR_cent[NSTAR] = {75,65,55,45,35,25,15,7.5,2.5};
const double STAR_v2[NSTAR]     = {0.0236, 0.03597, 0.04761, 0.0542, 0.05392, 0.04685, 0.03447, 0.02377, 0.01445};
const double STAR_v2_stat[NSTAR]= {0.0002, 0.00011, 6e-05, 4e-05, 3e-05, 3e-05, 3e-05, 7e-05, 0.00017};
const double STAR_v2_sys[NSTAR] = {0.0004, 0.00043, 0.00047, 0.00044, 0.00037, 0.00032, 0.00029, 0.00046, 0.00075};
const double STAR_g112OS[NSTAR]     = {-0.00012, -7e-06, 1.2e-05, -3e-06, 2e-06, -3e-06, -5e-06, -6e-06, 2e-06};
const double STAR_g112OS_stat[NSTAR]= {0.00015, 4.5e-05, 1.8e-05, 9e-06, 6e-06, 4e-06, 4e-06, 9e-06, 2e-05};
const double STAR_g112OS_sys[NSTAR] = {0.00016, 1.4e-05, 4e-06, 3e-06, 2e-06, 5e-06, 5e-06, 6e-06, 6e-06};
const double STAR_g112SS[NSTAR]     = {-0.00047, -0.00051, -0.000385, -0.000271, -0.000184, -0.000115, -6.1e-05, -4.1e-05, -2e-06};
const double STAR_g112SS_stat[NSTAR]= {0.00015, 4e-05, 1.8e-05, 9e-06, 6e-06, 4e-06, 4e-06, 9e-06, 2e-05};
const double STAR_g112SS_sys[NSTAR] = {4e-05, 6e-05, 1.1e-05, 3e-06, 3e-06, 3e-06, 5e-06, 0, 2.7e-05};
const double STAR_g132OS[NSTAR]     = {0.00021, 5.8e-05, 0.000111, 9.07e-05, 5.49e-05, 2.9e-05, 1.78e-05, -1e-05, 4e-06};
const double STAR_g132OS_stat[NSTAR]= {0.00015, 4.5e-05, 1.8e-05, 9.4e-06, 5.8e-06, 4e-06, 4.3e-06, 9e-06, 2e-05};
const double STAR_g132OS_sys[NSTAR] = {3e-05, 1.5e-05, 6e-06, 7e-07, 1.2e-06, 7e-06, 8e-07, 0, 4e-06};
const double STAR_g132SS[NSTAR]     = {0.00012, 4.1e-05, 5e-06, -3.5e-06, -4.1e-05, -3.6e-05, -2.9e-05, -2.1e-05, -5.4e-05};
const double STAR_g132SS_stat[NSTAR]= {0.00015, 4.5e-05, 1.8e-05, 9.4e-06, 6e-06, 4e-06, 4e-06, 9e-06, 2e-05};
const double STAR_g132SS_sys[NSTAR] = {3e-05, 1.1e-05, 6e-06, 1.2e-06, 8e-06, 2e-06, 4e-06, 2e-06, 6e-06};
const double STAR_dOS[NSTAR]     = {0.003334, 0.003034, 0.002213, 0.0015379, 0.001065, 0.0007462, 0.0005421, 0.0004478, 0.0003871};
const double STAR_dOS_stat[NSTAR]= {6e-06, 3e-06, 1.3e-06, 8e-07, 5e-07, 3e-07, 2e-07, 3e-07, 2e-07};
const double STAR_dOS_sys[NSTAR] = {3.5e-05, 5.8e-05, 4.38e-05, 2.86e-05, 2.15e-05, 1.84e-05, 2.39e-05, 3.49e-05, 4.39e-05};
const double STAR_dSS[NSTAR]     = {0.000674, 6.8e-05, -0.0001724, -0.0002541, -0.0002653, -0.0002283, -0.0001646, -0.0001219, -9.04e-05};
const double STAR_dSS_stat[NSTAR]= {6e-06, 3e-06, 1.3e-06, 8e-07, 5e-07, 3e-07, 2e-07, 3e-07, 2e-07};
const double STAR_dSS_sys[NSTAR] = {2.3e-05, 1.2e-05, 6.6e-06, 7.1e-06, 9.7e-06, 1.56e-05, 3.05e-05, 4.93e-05, 6.52e-05};
const double STAR_k112[NSTAR]     = {5, 4.7, 3.49, 2.759, 2.597, 2.45, 2.3, 2.6, 0.7};
const double STAR_k112_stat[NSTAR]= {3, 0.6, 0.23, 0.138, 0.115, 0.14, 0.3, 0.9, 4};
const double STAR_k112_sys[NSTAR] = {3, 0.5, 0.14, 0.017, 0.015, 0.09, 0.3, 0.7, 1.7};
const double STAR_k132[NSTAR]     = {1.6, 0.17, 0.94, 0.971, 1.34, 1.42, 1.91, 0.8, 8};
const double STAR_k132_stat[NSTAR]= {3.4, 0.6, 0.23, 0.138, 0.12, 0.14, 0.25, 0.9, 4};
const double STAR_k132_sys[NSTAR] = {0.5, 0.19, 0.05, 0.017, 0.15, 0.1, 0.04, 0.5, 0};
#endif
