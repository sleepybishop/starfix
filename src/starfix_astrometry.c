#include "starfix_astrometry.h"

#include <math.h>
#include <stdio.h>

#define J2000 51544.5
#define PI 3.14159265358979323846
#define TWOPI 6.28318530717958647692
#define degrad(x) ((x) * PI / 180.0)
#define raddeg(x) ((x) * 180.0 / PI)
#define arcsecrad(x) ((x) * 4.84813681109535993e-6)

double starfix_year_to_mjd(double year) { return J2000 + (year - 2000.0) * 365.25; }

double starfix_utc_to_mjd(int year, int month, int day, int hour, int minute, double second) {
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;
    long jdn = (long)day + (153L * m + 2L) / 5L + 365L * y + y / 4L - y / 100L + y / 400L - 32045L;
    double mjd_day = (double)(jdn - 2400001L);
    return mjd_day + ((double)hour + (double)minute / 60.0 + second / 3600.0) / 24.0;
}

static void precession_matrix(double equinox, double P[3][3]) {
    if (fabs(equinox - J2000) < 1.0) {
        P[0][0] = P[1][1] = P[2][2] = 1.0;
        P[0][1] = P[0][2] = P[1][0] = P[1][2] = P[2][0] = P[2][1] = 0.0;
        return;
    }

    double T = (equinox - J2000) / 36525.0;
    /* coefficients: USNO/HMNAO, rutherford appleton laboratory, ISSN 0737-6421 */
    double zeta_A = arcsecrad(
        2.650545 +
        T * (2306.083227 + T * (0.2988499 + T * (0.01801828 + T * (-5.971e-6 - 3.173e-7 * T)))));
    double z_A = arcsecrad(
        -2.650545 +
        T * (2306.077181 + T * (1.0927348 + T * (0.01826837 + T * (-28.596e-6 - 2.904e-7 * T)))));
    double theta_A = arcsecrad(
        T * (2004.191903 - T * (0.4294934 + T * (0.04182264 + T * (7.089e-6 + 1.274e-7 * T)))));

    double cz = cos(zeta_A), sz = sin(zeta_A);
    double cza = cos(z_A), sza = sin(z_A);
    double ct = cos(theta_A), st = sin(theta_A);

    P[0][0] = cz * ct * cza - sz * sza;
    P[0][1] = -sz * ct * cza - cz * sza;
    P[0][2] = -st * cza;

    P[1][0] = cz * ct * sza + sz * cza;
    P[1][1] = -sz * ct * sza + cz * cza;
    P[1][2] = -st * sza;

    P[2][0] = cz * st;
    P[2][1] = -sz * st;
    P[2][2] = ct;
}
static void obliquity(double mj, double* eps) {
    double t = (mj - J2000) / 36525.0;
    *eps = degrad(23.4392911 + t * (-46.8150 + t * (-0.00059 + t * (0.001813))) / 3600.0);
}
/* IAU 1980 nutation series; stern */

#define NUT_SCALE 1e4
#define NUT_SERIES 106
#define NUT_MAXMUL 4
#define SECPERCIRC (3600. * 360.)

/* Delaunay arguments, in arc seconds; they differ slightly from ELP82B */
static double delaunay[5][4] = {
    {485866.733, 1717915922.633, 31.310, 0.064},  /* M', moon mean anom */
    {1287099.804, 129596581.224, -0.577, -0.012}, /* M, sun mean anom */
    {335778.877, 1739527263.137, -13.257, 0.011}, /* F, moon arg lat */
    {1072261.307, 1602961601.328, -6.891, 0.019}, /* D, elong moon sun */
    {450160.280, -6962890.539, 7.455, 0.008},     /* Om, moon l asc node */
};

/* multipliers for Delaunay arguments */
static short multarg[NUT_SERIES][5] = {
    /* bounds:  -2..3, -2..2, -2/0/2/4, -4..4, 0..2 */
    {0, 0, 0, 0, 1},   {0, 0, 0, 0, 2},   {-2, 0, 2, 0, 1},  {2, 0, -2, 0, 0},  {-2, 0, 2, 0, 2},
    {1, -1, 0, -1, 0}, {0, -2, 2, -2, 1}, {2, 0, -2, 0, 1},  {0, 0, 2, -2, 2},  {0, 1, 0, 0, 0},
    {0, 1, 2, -2, 2},  {0, -1, 2, -2, 2}, {0, 0, 2, -2, 1},  {2, 0, 0, -2, 0},  {0, 0, 2, -2, 0},
    {0, 2, 0, 0, 0},   {0, 1, 0, 0, 1},   {0, 2, 2, -2, 2},  {0, -1, 0, 0, 1},  {-2, 0, 0, 2, 1},
    {0, -1, 2, -2, 1}, {2, 0, 0, -2, 1},  {0, 1, 2, -2, 1},  {1, 0, 0, -1, 0},  {2, 1, 0, -2, 0},
    {0, 0, -2, 2, 1},  {0, 1, -2, 2, 0},  {0, 1, 0, 0, 2},   {-1, 0, 0, 1, 1},  {0, 1, 2, -2, 0},
    {0, 0, 2, 0, 2},   {1, 0, 0, 0, 0},   {0, 0, 2, 0, 1},   {1, 0, 2, 0, 2},   {1, 0, 0, -2, 0},
    {-1, 0, 2, 0, 2},  {0, 0, 0, 2, 0},   {1, 0, 0, 0, 1},   {-1, 0, 0, 0, 1},  {-1, 0, 2, 2, 2},
    {1, 0, 2, 0, 1},   {0, 0, 2, 2, 2},   {2, 0, 0, 0, 0},   {1, 0, 2, -2, 2},  {2, 0, 2, 0, 2},
    {0, 0, 2, 0, 0},   {-1, 0, 2, 0, 1},  {-1, 0, 0, 2, 1},  {1, 0, 0, -2, 1},  {-1, 0, 2, 2, 1},
    {1, 1, 0, -2, 0},  {0, 1, 2, 0, 2},   {0, -1, 2, 0, 2},  {1, 0, 2, 2, 2},   {1, 0, 0, 2, 0},
    {2, 0, 2, -2, 2},  {0, 0, 0, 2, 1},   {0, 0, 2, 2, 1},   {1, 0, 2, -2, 1},  {0, 0, 0, -2, 1},
    {1, -1, 0, 0, 0},  {2, 0, 2, 0, 1},   {0, 1, 0, -2, 0},  {1, 0, -2, 0, 0},  {0, 0, 0, 1, 0},
    {1, 1, 0, 0, 0},   {1, 0, 2, 0, 0},   {1, -1, 2, 0, 2},  {-1, -1, 2, 2, 2}, {-2, 0, 0, 0, 1},
    {3, 0, 2, 0, 2},   {0, -1, 2, 2, 2},  {1, 1, 2, 0, 2},   {-1, 0, 2, -2, 1}, {2, 0, 0, 0, 1},
    {1, 0, 0, 0, 2},   {3, 0, 0, 0, 0},   {0, 0, 2, 1, 2},   {-1, 0, 0, 0, 2},  {1, 0, 0, -4, 0},
    {-2, 0, 2, 2, 2},  {-1, 0, 2, 4, 2},  {2, 0, 0, -4, 0},  {1, 1, 2, -2, 2},  {1, 0, 2, 2, 1},
    {-2, 0, 2, 4, 2},  {-1, 0, 4, 0, 2},  {1, -1, 0, -2, 0}, {2, 0, 2, -2, 1},  {2, 0, 2, 2, 2},
    {1, 0, 0, 2, 1},   {0, 0, 4, -2, 2},  {3, 0, 2, -2, 2},  {1, 0, 2, -2, 0},  {0, 1, 2, 0, 1},
    {-1, -1, 0, 2, 1}, {0, 0, -2, 0, 1},  {0, 0, 2, -1, 2},  {0, 1, 0, 2, 0},   {1, 0, -2, -2, 0},
    {0, -1, 2, 0, 1},  {1, 1, 0, -2, 1},  {1, 0, -2, 2, 0},  {2, 0, 0, 2, 0},   {0, 0, 2, 4, 2},
    {0, 1, 0, 1, 0}};

/* 1/NUT_SCALE arcsec: {index, dpsi, dpsi per T/10, deps, deps per T/10} */
static long ampsecul[][5] = {{0, -171996, -1742, 92025, 89},
                             {1, 2062, 2, -895, 5},
                             {8, -13187, -16, 5736, -31},
                             {9, 1426, -34, 54, -1},
                             {10, -517, 12, 224, -6},
                             {11, 217, -5, -95, 3},
                             {12, 129, 1, -70, 0},
                             {15, 17, -1, 0, 0},
                             {17, -16, 1, 7, 0},
                             {30, -2274, -2, 977, -5},
                             {31, 712, 1, -7, 0},
                             {32, -386, -4, 200, 0},
                             {33, -301, 0, 129, -1},
                             {37, 63, 1, -33, 0},
                             {38, -58, -1, 32, 0},
                             /* termination */
                             {
                                 -1,
                             }};

/* {dpsi, deps} in 1/NUT_SCALE arcsec; secular entries are zeroed */
static short ampconst[NUT_SERIES][2] = {
    {0, 0},    {0, 0},    {46, -24},  {11, 0},    {-3, 1},   {-3, 0},  {-2, 1},   {1, 0},
    {0, 0},    {0, 0},    {0, 0},     {0, 0},     {0, 0},    {48, 1},  {-22, 0},  {0, 0},
    {-15, 9},  {0, 0},    {-12, 6},   {-6, 3},    {-5, 3},   {4, -2},  {4, -2},   {-4, 0},
    {1, 0},    {1, 0},    {-1, 0},    {1, 0},     {1, 0},    {-1, 0},  {0, 0},    {0, 0},
    {0, 0},    {0, 0},    {-158, -1}, {123, -53}, {63, -2},  {0, 0},   {0, 0},    {-59, 26},
    {-51, 27}, {-38, 16}, {29, -1},   {29, -12},  {-31, 13}, {26, -1}, {21, -10}, {16, -8},
    {-13, 7},  {-10, 5},  {-7, 0},    {7, -3},    {-7, 3},   {-8, 3},  {6, 0},    {6, -3},
    {-6, 3},   {-7, 3},   {6, -3},    {-5, 3},    {5, 0},    {-5, 3},  {-4, 0},   {4, 0},
    {-4, 0},   {-3, 0},   {3, 0},     {-3, 1},    {-3, 1},   {-2, 1},  {-3, 1},   {-3, 1},
    {2, -1},   {-2, 1},   {2, -1},    {-2, 1},    {2, 0},    {2, -1},  {1, -1},   {-1, 0},
    {1, -1},   {-2, 1},   {-1, 0},    {1, -1},    {-1, 1},   {-1, 1},  {1, 0},    {1, 0},
    {1, -1},   {-1, 0},   {-1, 0},    {1, 0},     {1, 0},    {-1, 0},  {1, 0},    {1, 0},
    {-1, 0},   {-1, 0},   {-1, 0},    {-1, 0},    {-1, 0},   {-1, 0},  {-1, 0},   {1, 0},
    {-1, 0},   {1, 0}};

/* compute nutation in obliquity (deps) and longitude (dpsi), in radians */
static void nutation(double mj, double* deps, double* dpsi) {
    double T, T2, T3, T10; /* jul cent since J2000 */
    int i, isecul;         /* index in term table */
    double delcache[5][2 * NUT_MAXMUL + 1];
    /* delaunay multiples in radians, indexed by multiple + NUT_MAXMUL */

    T = (mj - J2000) / 36525.;
    T2 = T * T;
    T3 = T2 * T;
    T10 = T / 10.;

    /* calculate delaunay args and place in cache */
    for (i = 0; i < 5; ++i) {
        double x;
        short j;

        x = delaunay[i][0] + delaunay[i][1] * T + delaunay[i][2] * T2 + delaunay[i][3] * T3;

        /* convert to radians */
        x /= SECPERCIRC;
        x -= floor(x);
        x *= 2. * PI;

        /* fill cache table */
        for (j = 0; j <= 2 * NUT_MAXMUL; ++j) delcache[i][j] = (j - NUT_MAXMUL) * x;
    }

    /* find dpsi and deps */
    double sum_dpsi = 0., sum_deps = 0.;
    for (i = isecul = 0; i < NUT_SERIES; ++i) {
        double arg = 0., ampsin, ampcos;
        short j;

        if (ampconst[i][0] || ampconst[i][1]) {
            /* take non-secular terms from simple array */
            ampsin = ampconst[i][0];
            ampcos = ampconst[i][1];
        } else {
            /* secular terms from different array */
            ampsin = (double)ampsecul[isecul][1] + (double)ampsecul[isecul][2] * T10;
            ampcos = (double)ampsecul[isecul][3] + (double)ampsecul[isecul][4] * T10;
            ++isecul;
        }

        for (j = 0; j < 5; ++j) arg += delcache[j][NUT_MAXMUL + multarg[i][j]];

        sum_dpsi += ampsin * sin(arg);

        sum_deps += ampcos * cos(arg);
    }

    *dpsi = degrad(sum_dpsi / 3600. / NUT_SCALE);
    *deps = degrad(sum_deps / 3600. / NUT_SCALE);
}

void starfix_astrometry_frame(double mj, double P[3][3], double aberration[3]) {
    double prec[3][3], nut[3][3];
    double epsilon, dpsi, deps;

    precession_matrix(mj, prec);
    obliquity(mj, &epsilon);
    nutation(mj, &deps, &dpsi);

    double se = sin(epsilon), ce = cos(epsilon);
    double sp = sin(dpsi), cp = cos(dpsi);
    double sede = sin(epsilon + deps), cede = cos(epsilon + deps);

    nut[0][0] = cp;
    nut[0][1] = -sp * ce;
    nut[0][2] = -sp * se;

    nut[1][0] = cede * sp;
    nut[1][1] = cede * cp * ce + sede * se;
    nut[1][2] = cede * cp * se - sede * ce;

    nut[2][0] = sede * sp;
    nut[2][1] = sede * cp * ce - cede * se;
    nut[2][2] = sede * cp * se + cede * ce;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            P[i][j] = nut[i][0] * prec[0][j] + nut[i][1] * prec[1][j] + nut[i][2] * prec[2][j];
        }
    }

    double L = TWOPI * (0.27908 + 100.00214 * (mj - J2000) / 36525.0);
    double cL = cos(L);
    aberration[0] = -0.994e-4 * sin(L);
    aberration[1] = 0.912e-4 * cL;
    aberration[2] = 0.395e-4 * cL;
}

double starfix_gast_deg(double mjd_utc) {
    double d_ut1 = mjd_utc - J2000;
    double t = d_ut1 / 36525.0;
    double gmst =
        280.46061837 + 360.98564736629 * d_ut1 + 0.000387933 * t * t - (t * t * t) / 38710000.0;
    double eps = 0.0, dpsi = 0.0, deps = 0.0;
    obliquity(mjd_utc, &eps);
    nutation(mjd_utc, &deps, &dpsi);
    double eq_equinox_deg = raddeg(dpsi * cos(eps + deps));
    double gast = fmod(gmst + eq_equinox_deg, 360.0);
    if (gast < 0.0) gast += 360.0;
    return gast;
}

static void unrefractLT15(double pr, double tr, double aa, double* ta);
static void unrefractGE15(double pr, double tr, double aa, double* ta);

void starfix_unrefract(double pr, double tr, double aa, double* ta) {
#define LTLIM 14.5
#define GELIM 15.5

    double aadeg = raddeg(aa);

    if (aadeg < LTLIM)
        unrefractLT15(pr, tr, aa, ta);
    else if (aadeg >= GELIM)
        unrefractGE15(pr, tr, aa, ta);
    else {
        /* smooth blend -- important for inverse */
        double taLT, taGE, p;

        unrefractLT15(pr, tr, aa, &taLT);
        unrefractGE15(pr, tr, aa, &taGE);
        p = (aadeg - LTLIM) / (GELIM - LTLIM);
        *ta = taLT + (taGE - taLT) * p;
    }
}

static void unrefractGE15(double pr, double tr, double aa, double* ta) {
    double r;

    r = 7.888888e-5 * pr / ((273 + tr) * tan(aa));
    *ta = aa - r;
}

static void unrefractLT15(double pr, double tr, double aa, double* ta) {
    double aadeg = raddeg(aa);
    double r, a, b;

    a = ((2e-5 * aadeg + 1.96e-2) * aadeg + 1.594e-1) * pr;
    b = (273 + tr) * ((8.45e-2 * aadeg + 5.05e-1) * aadeg + 1);
    r = degrad(a / b);

    *ta = (aa < 0 && r < 0) ? aa : aa - r; /* 0 below ~5 degs */
}

void starfix_refract(double pr, double tr, double ta, double* aa) {
    if (ta <= 0.0) {
        *aa = ta;
        return;
    }
    double est_aa = ta;
    for (int iter = 0; iter < 3; iter++) {
        double est_ta = est_aa;
        starfix_unrefract(pr, tr, est_aa, &est_ta);
        est_aa += (ta - est_ta);
    }
    *aa = est_aa;
}
