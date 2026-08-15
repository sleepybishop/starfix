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

/* Precession computation to or from J2000
 * If dir < 0 go from equinox to J2000
 * If dir > 0 go from J2000 to equinox
 */
static void do_precess(double equinox, double* alpha, double* delta, int dir) {
    double T;
    double zeta_A, z_A, theta_A;
    double A, B, C;

    /* Precession progresses about 1 arc second in .047 years = 17 days
     * Threshold must be at most 1.7 days for 0.1 arcsec accuracy */
    double threshhold = 1;

    if (fabs(equinox - J2000) < threshhold) return;

    T = (equinox - J2000) / 36525.0;

    /* Compute coefficients in arcseconds from the Naval
     * Observatory, Her Majesty's Nautical Amanac Office;
     * Rutherford Appleton Laboratory, ISSN 0737-6421. */

    zeta_A = arcsecrad(
        2.650545 +
        T * (2306.083227 + T * (0.2988499 + T * (0.01801828 + T * (-5.971e-6 - 3.173e-7 * T)))));
    z_A = arcsecrad(
        -2.650545 +
        T * (2306.077181 + T * (1.0927348 + T * (0.01826837 + T * (-28.596e-6 - 2.904e-7 * T)))));
    theta_A = arcsecrad(
        T * (2004.191903 - T * (0.4294934 + T * (0.04182264 + T * (7.089e-6 + 1.274e-7 * T)))));

    if (dir < 0) {
        /* If dir is negative, go from equinox to 2000.0 */
        A = sin(*alpha - z_A) * cos(*delta);
        B = cos(*alpha - z_A) * cos(theta_A) * cos(*delta) + sin(theta_A) * sin(*delta);
        C = -cos(*alpha - z_A) * sin(theta_A) * cos(*delta) + cos(theta_A) * sin(*delta);
        *alpha = atan2(A, B) - zeta_A;
    } else {
        /* If dir is positive, go from 2000.0 to equinox */
        A = sin(*alpha + zeta_A) * cos(*delta);
        B = cos(*alpha + zeta_A) * cos(theta_A) * cos(*delta) - sin(theta_A) * sin(*delta);
        C = cos(*alpha + zeta_A) * sin(theta_A) * cos(*delta) + cos(theta_A) * sin(*delta);
        *alpha = atan2(A, B) + z_A;
    }
    while (*alpha < 0.0) {
        *alpha += TWOPI;
    }
    while (*alpha >= TWOPI) {
        *alpha -= TWOPI;
    }
    *delta = asin(C);
}

/* Corrects ra and dec, both in radians, for precession from epoch 1 to
 * epoch 2.  The epochs are given by their modified JDs, from_mjd and
 * to_mjd, respectively.  N.B. ra and dec are modifed IN PLACE.
 */
void starfix_precess(double to_mjd, double* ra, double* dec) { do_precess(to_mjd, ra, dec, 1); }
static void obliquity(double mj, double* eps) {
    double t = (mj - J2000) / 36525.0;
    *eps = degrad(23.4392911 + t * (-46.8150 + t * (-0.00059 + t * (0.001813))) / 3600.0);
}
/* nutation (in IAU (1980) expression) and abberation; stern
 * on an HP PA processor, this reproduces the Almanac nutation values
 * (given to 0.001") EXACTLY over 750 days (1995 and 1996)
 */

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

/* amplitudes which  have secular terms; in 1/NUT_SCALE arc seconds
 * {index, constant dPSI, T/10 in dPSI, constant in dEPS, T/10 in dEPS}
 */
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

/* amplitudes which only have constant terms; same unit as above
 * {dPSI, dEPS}
 * indexes which are already in ampsecul[][] are zeroed
 */
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

/* given the modified JD, mj, find the nutation in obliquity, *deps, and
 * the nutation in longitude, *dpsi, each in radians.
 */
void nutation(double mj, double* deps, /* on input:  precision parameter in arc seconds */
              double* dpsi) {
    static double lastmj = -10000, lastdeps, lastdpsi;
    double T, T2, T3, T10; /* jul cent since J2000 */
    double prec;           /* series precis in arc sec */
    int i, isecul;         /* index in term table */
    static double delcache[5][2 * NUT_MAXMUL + 1];
    /* cache for multiples of delaunay args
     * [M',M,F,D,Om][-min*x, .. , 0, .., max*x]
     * make static to have unfilled fields cleared on init
     */

    if (fabs(mj - lastmj) < 1e-9) {
        *deps = lastdeps;
        *dpsi = lastdpsi;
        return;
    }

    prec = 0.0;

#if 0 /* this is if deps should contain a precision value */
	prec =* deps;
	if (prec < 0.0 || prec > 1.0)	/* accept only sane value */
		prec = 1.0;
#endif

    /* augment for abundance of small terms */
    prec *= NUT_SCALE / 10;

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
    lastdpsi = lastdeps = 0.;
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

        if (fabs(ampsin) >= prec) lastdpsi += ampsin * sin(arg);

        if (fabs(ampcos) >= prec) lastdeps += ampcos * cos(arg);
    }

    /* convert to radians.
     */
    lastdpsi = degrad(lastdpsi / 3600. / NUT_SCALE);
    lastdeps = degrad(lastdeps / 3600. / NUT_SCALE);

    lastmj = mj;
    *deps = lastdeps;
    *dpsi = lastdpsi;
}

/* given the modified JD, mj, correct, IN PLACE, the right ascension *ra
 * and declination *dec (both in radians) for nutation.
 */
void starfix_nutation(double mj, double* ra, double* dec) {
    static double lastmj = -10000;
    static double a[3][3]; /* rotation matrix */
    double xold, yold, zold, x, y, z;

    if (fabs(mj - lastmj) > 1e-9) {
        double epsilon, dpsi, deps;
        double se, ce, sp, cp, sede, cede;

        obliquity(mj, &epsilon);
        nutation(mj, &deps, &dpsi);

        /* the rotation matrix a applies the nutation correction to
         * a vector of equatoreal coordinates Xeq to Xeq' by 3 subsequent
         * rotations:  R1 - from equatoreal to ecliptic system by
         * rotation of angle epsilon about x, R2 - rotate ecliptic
         * system by -dpsi about its z, R3 - from ecliptic to equatoreal
         * by rotation of angle -(epsilon + deps)
         *
         *	Xeq' = A * Xeq = R3 * R2 * R1 * Xeq
         *
         *		[ 1       0          0    ]
         * R1 =	[ 0   cos(eps)   sin(eps) ]
         *		[ 0  - sin(eps)  cos(eps) ]
         *
         *		[ cos(dpsi)  - sin(dpsi)  0 ]
         * R2 =	[ sin(dpsi)   cos(dpsi)   0 ]
         *		[      0           0      1 ]
         *
         *		[ 1         0                 0         ]
         * R3 =	[ 0  cos(eps + deps)  - sin(eps + deps) ]
         *		[ 0  sin(eps + deps)   cos(eps + deps)  ]
         *
         * for efficiency, here is a explicitely:
         */

        se = sin(epsilon);
        ce = cos(epsilon);
        sp = sin(dpsi);
        cp = cos(dpsi);
        sede = sin(epsilon + deps);
        cede = cos(epsilon + deps);

        a[0][0] = cp;
        a[0][1] = -sp * ce;
        a[0][2] = -sp * se;

        a[1][0] = cede * sp;
        a[1][1] = cede * cp * ce + sede * se;
        a[1][2] = cede * cp * se - sede * ce;

        a[2][0] = sede * sp;
        a[2][1] = sede * cp * ce - cede * se;
        a[2][2] = sede * cp * se + cede * ce;

        lastmj = mj;
    }

    double cd = cos(*dec);
    xold = cos(*ra) * cd;
    yold = sin(*ra) * cd;
    zold = sin(*dec);
    x = a[0][0] * xold + a[0][1] * yold + a[0][2] * zold;
    y = a[1][0] * xold + a[1][1] * yold + a[1][2] * zold;
    z = a[2][0] * xold + a[2][1] * yold + a[2][2] * zold;
    double r = sqrt(x * x + y * y + z * z);
    *dec = asin(z / r);
    *ra = atan2(y, x);            /* radius should be 1.0 */
    if (*ra < 0.) *ra += 2. * PI; /* make positive for display */
}

void starfix_aberration(double mj, double* ra, double* dec) {
    double x, y, z;
    double vx, vy, vz;
    double L;
    double cL;
    double r;

    double cd = cos(*dec);
    x = cos(*ra) * cd;
    y = sin(*ra) * cd;
    z = sin(*dec);

    L = TWOPI * (0.27908 + 100.00214 * (mj - J2000) / 36525.0);
    cL = cos(L);
    vx = -0.994e-4 * sin(L);
    vy = 0.912e-4 * cL;
    vz = 0.395e-4 * cL;
    x += vx;
    y += vy;
    z += vz;

    r = sqrt(x * x + y * y + z * z);
    *dec = asin(z / r);
    *ra = atan2(y, x);
    if (*ra < 0) *ra += TWOPI;
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

/* correct the true altitude, ta, for refraction to the apparent altitude, aa,
 * each in radians, given the local atmospheric pressure, pr, in mbars, and
 * the temperature, tr, in degrees C.
 */
void starfix_refract(double pr, double tr, double ta, double* aa) {
#define MAXRERR degrad(0.1 / 3600.) /* desired accuracy, rads */

    double d, t, t0, a;

    if (isnan(ta)) {
        *aa = ta;
        return;
    }

    /* first guess of error is to go backwards.
     * make use that we know delta-apparent is always < delta-true.
     */
    starfix_unrefract(pr, tr, ta, &t);
    d = 0.8 * (ta - t);
    t0 = t;
    a = ta;

    /* use secant method to discover a value that unrefracts to ta.
     * max=7 ave=2.4 loops in hundreds of test cases.
     */
    while (1) {
        a += d;
        starfix_unrefract(pr, tr, a, &t);
        if (fabs(ta - t) <= MAXRERR) break;
        d *= -(ta - t) / (t0 - t);
        t0 = t;
    }

    *aa = a;

#undef MAXRERR
}
