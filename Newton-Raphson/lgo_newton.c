/*
 * lgo_newton.c — Newton-Raphson solver for LGO apogee radius
 *
 * Drop-in replacement for the grid-search find_hindex() used in
 * LGO_C2.c / LGO_GUI.c.
 *
 * ── Algorithm ───────────────────────────────────────────────────────────
 *  For each target m/n value:
 *   1. Coarse log-spaced scan (N_SCAN=30) to bracket the root.
 *   2. Newton-Raphson step with central-difference derivative.
 *   3. If Newton step leaves the bracket, bisect (Brent-style fallback).
 *   4. Stop when |m_n(rA) - target| < TOL.
 *
 * ── Integration ─────────────────────────────────────────────────────────
 *  Either:
 *   (a) #include "lgo_newton.c" at the top of LGO_C2.c / LGO_GUI.c
 *   (b) Add to the same compilation unit alongside the existing code.
 *
 *  Then replace every find_hindex() call site (see BEFORE/AFTER below).
 *
 * ── Build (standalone test) ──────────────────────────────────────────────
 *   gcc -O2 -DLGO_NEWTON_TEST -o lgo_newton_test lgo_newton.c -lm
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

/* ── Constants (must match LGO_C2.c / LGO_GUI.c) ───────────────────────── */
#ifndef MU
#  define MU    398600.4
#  define OME   7.2921159e-5
#  define RMEAN 6371.0
#  define RE    6378.0
#  define J2    0.0010827
#  define PI    3.14159265358979323846
#  define sq(x) ((x)*(x))
#endif

/* ── Solver parameters ──────────────────────────────────────────────────── */
#define NEWTON_TOL      1e-9    /* convergence: |m_n - target|               */
#define NEWTON_MAXITER  80      /* max iterations per target                 */
#define NEWTON_NSCAN    30      /* coarse scan points for bracket finding     */
#define NEWTON_RA_LO    (RE + 300.0)          /* [km] lower bound            */
#define NEWTON_RA_HI    (RMEAN + 60000.0)     /* [km] upper bound            */

/* ══════════════════════════════════════════════════════════════════════════
 *  eval_mn_at_rA
 *
 *  Evaluates m/n repetition factor at a given apogee radius.
 *
 *  rA       [km]  apogee geocentric radius
 *  i_rad    [rad] inclination
 *  phi_rad  [rad] apogee sub-satellite latitude  (= i_rad for special LGO)
 *  cosg     [-]   cos(gamma) = cos(i)/cos(phi)   (= 1 for special LGO)
 *  special  [0/1] 1 = special LGO formula,  0 = general LGO formula
 *
 *  Returns m/n [-], or NAN for physically invalid inputs.
 * ══════════════════════════════════════════════════════════════════════════ */
static double eval_mn_at_rA(double rA,
                             double i_rad, double phi_rad,
                             double cosg,  int special)
{
    double ci  = cos(i_rad),   ci2  = ci*ci;
    double cp  = cos(phi_rad), cp2  = cp*cp;
    double rA3 = rA*rA*rA;

    /* ── LGO perigee radius ─────────────────────────────────────────────── */
    double den = special ? (2.0*MU - sq(OME)*rA3*ci2)
                         : (2.0*MU*cosg - sq(OME)*rA3*cp2);

    if (den <= 0.0 || !isfinite(den)) return NAN;

    double rp = sq(OME)*rA3*rA*ci2 / den;    /* [km] perigee geocentric r   */

    if (rp <= 0.0 || rp >= rA) return NAN;

    /* ── Orbital elements ───────────────────────────────────────────────── */
    double a        = (rA + rp) / 2.0;       /* [km] semimajor axis         */
    double e        = (rA - rp) / (rA + rp); /* [-]  eccentricity           */
    double one_m_e2 = 1.0 - e*e;

    if (a <= 0.0 || one_m_e2 <= 0.0) return NAN;

    double a2   = a*a;
    double a3h  = pow(a, 1.5);               /* a^(3/2)  [km^(3/2)]         */
    double a7h  = pow(a, 3.5);               /* a^(7/2)  [km^(7/2)]         */
    double e2sq = sq(one_m_e2);              /* (1-e^2)^2                   */
    double sqe  = sqrt(one_m_e2);            /* (1-e^2)^(1/2)               */

    /* ── J2 nodal correction K ──────────────────────────────────────────── */
    double K = 1.0 - (1.5*J2*RE*RE / (a2*e2sq)) *
        (  (10.0*ci2 - 2.0)/4.0
         - (3.0*ci2  - 1.0)/4.0 * sqe
         + sqrt(MU)*ci / (OME*a7h*e2sq) );

    double mn = 1.0 / (OME * a3h/sqrt(MU) * K);   /* [-] m/n factor        */

    return (isfinite(mn) && mn > 0.0) ? mn : NAN;
}

/* ══════════════════════════════════════════════════════════════════════════
 *  newton_solve_rA
 *
 *  Finds rA [km] such that m_n(rA) == mn_target using Newton-Raphson with
 *  bisection fallback.
 *
 *  mn_target  [-]    target m/n repetition factor
 *  i_rad      [rad]  inclination
 *  phi_rad    [rad]  apogee sub-satellite latitude  (= i_rad for special)
 *  special    [0/1]  1 = special LGO,  0 = general LGO
 *
 *  rA_out     [km]   solution apogee radius  (NAN on failure)
 *  Va_out     [km/s] apogee velocity at solution
 *  rP_out     [km]   perigee radius at solution
 *
 *  Returns true on convergence, false if no solution found.
 * ══════════════════════════════════════════════════════════════════════════ */
bool newton_solve_rA(double mn_target,
                     double i_rad, double phi_rad,
                     int    special,
                     double *rA_out, double *Va_out, double *rP_out)
{
    *rA_out = NAN;
    if (Va_out) *Va_out = NAN;
    if (rP_out) *rP_out = NAN;

    if (mn_target <= 0.0 || !isfinite(mn_target)) return false;

    /* cos(gamma) = cos(i)/cos(phi)  →  1 for special case (phi=i)          */
    double cosg = special ? 1.0 : cos(i_rad)/cos(phi_rad);

    /* ── Step 1: coarse log-spaced scan to find a bracket ──────────────── */
    double log_lo  = log(NEWTON_RA_LO);
    double log_hi  = log(NEWTON_RA_HI);
    double dlog    = (log_hi - log_lo) / (NEWTON_NSCAN - 1);

    double rA_lo = NAN, rA_hi = NAN;
    double mn_lo = NAN, mn_hi;
    double prev_rA = NAN;

    for (int s = 0; s < NEWTON_NSCAN; s++) {
        double rA = exp(log_lo + s*dlog);
        double mn = eval_mn_at_rA(rA, i_rad, phi_rad, cosg, special);

        if (isfinite(mn) && isfinite(mn_lo)) {
            double f_prev = mn_lo - mn_target;
            double f_curr = mn   - mn_target;
            if (f_prev * f_curr <= 0.0) {      /* sign change → bracket      */
                rA_lo = prev_rA;
                rA_hi = rA;
                break;
            }
        }
        if (isfinite(mn)) { mn_lo = mn; prev_rA = rA; }
    }

    if (!isfinite(rA_lo)) return false;        /* no bracket → no solution   */

    /* ── Step 2: Newton-Raphson with bisection fallback ────────────────── */
    double rA = (rA_lo + rA_hi) / 2.0;        /* start at bracket midpoint  */

    for (int iter = 0; iter < NEWTON_MAXITER; iter++) {
        double mn_val = eval_mn_at_rA(rA, i_rad, phi_rad, cosg, special);

        if (!isfinite(mn_val)) {               /* invalid → bisect           */
            rA = (rA_lo + rA_hi) / 2.0;
            continue;
        }

        double f = mn_val - mn_target;

        if (fabs(f) < NEWTON_TOL) break;       /* ── converged ─────────── */

        /* Keep bracket tight for fallback */
        if (f > 0.0) rA_lo = rA;
        else         rA_hi = rA;

        /* Central-difference derivative  dm_n/drA */
        double h    = rA * 1e-5;
        double mn_p = eval_mn_at_rA(rA+h, i_rad, phi_rad, cosg, special);
        double mn_m = eval_mn_at_rA(rA-h, i_rad, phi_rad, cosg, special);
        double df;

        if      (isfinite(mn_p) && isfinite(mn_m)) df = (mn_p - mn_m) / (2.0*h);
        else if (isfinite(mn_p))                    df = (mn_p - mn_val) / h;
        else if (isfinite(mn_m))                    df = (mn_val - mn_m) / h;
        else { rA = (rA_lo + rA_hi) / 2.0; continue; }

        if (fabs(df) < 1e-20) {                /* flat region → bisect       */
            rA = (rA_lo + rA_hi) / 2.0;
            continue;
        }

        double rA_newton = rA - f / df;

        /* Accept Newton step only if it stays inside the bracket            */
        rA = (rA_newton > rA_lo && rA_newton < rA_hi)
             ? rA_newton
             : (rA_lo + rA_hi) / 2.0;         /* bisect fallback            */
    }

    /* ── Step 3: verify and extract orbit parameters ────────────────────── */
    double mn_final = eval_mn_at_rA(rA, i_rad, phi_rad, cosg, special);
    if (!isfinite(mn_final) || fabs(mn_final - mn_target) > 1e-4) return false;

    /* Recompute rP from the solution rA */
    double ci2 = sq(cos(i_rad));
    double cp2 = sq(cos(phi_rad));
    double rA3 = rA*rA*rA;
    double den = special ? (2.0*MU - sq(OME)*rA3*ci2)
                         : (2.0*MU*cosg - sq(OME)*rA3*cp2);
    double rp = sq(OME)*rA3*rA*ci2 / den;

    *rA_out = rA;
    if (Va_out) *Va_out = sqrt(2.0*MU*rp / (rA*(rA+rp)));
    if (rP_out) *rP_out = rp;
    return true;
}

/* ══════════════════════════════════════════════════════════════════════════
 *  newton_solve_all
 *
 *  Convenience wrapper: solves for all n_targets values in mn0[].
 *
 *  rA_sol[t], Va_sol[t], rP_sol[t] receive results;
 *  set to NAN if the t-th target has no solution.
 *  Returns number of targets successfully solved.
 * ══════════════════════════════════════════════════════════════════════════ */
int newton_solve_all(const double *mn0,   int n_targets,
                     double i_rad,  double phi_rad,  int special,
                     double *rA_sol, double *Va_sol, double *rP_sol)
{
    int n_solved = 0;
    for (int t = 0; t < n_targets; t++) {
        rA_sol[t] = Va_sol[t] = rP_sol[t] = NAN;
        if (newton_solve_rA(mn0[t], i_rad, phi_rad, special,
                            &rA_sol[t], &Va_sol[t], &rP_sol[t]))
            n_solved++;
    }
    return n_solved;
}


/* ══════════════════════════════════════════════════════════════════════════
 *  HOW TO REPLACE find_hindex() IN LGO_C2.c / LGO_GUI.c
 *
 *  ── BEFORE (special LGO, inside the Z==0 branch) ─────────────────────
 *
 *    // ... grid arrays mn_k[N_H], Va_k[N_H] already filled ...
 *    int *hindex = malloc(n_mn0 * sizeof(int));
 *    find_hindex(mn_k, N_H, mn0, n_mn0, hindex);
 *
 *    for (int t = 0; t < n_mn0; t++) {
 *        int h = hindex[t];
 *        if (h < 0 || !isfinite(Va_k[h])) continue;
 *        DATA[J].H   = H[h];
 *        DATA[J].Va  = Va_k[h];
 *        DATA[J].mn  = mn0[t];
 *        DATA[J].i_rad = i_rad;
 *        DATA[J].phi   = i_rad;
 *        DATA[J].has_phi = 0;
 *        J++;
 *    }
 *    free(hindex);
 *
 *  ── AFTER (drop-in replacement) ──────────────────────────────────────
 *
 *    double rA_s[16], Va_s[16], rP_s[16];   // max 16 targets
 *    newton_solve_all(mn0, n_mn0, i_rad, i_rad, 1,
 *                     rA_s, Va_s, rP_s);
 *
 *    for (int t = 0; t < n_mn0; t++) {
 *        if (!isfinite(rA_s[t])) continue;
 *        DATA[J].H     = rA_s[t] - RMEAN;   // altitude [km]
 *        DATA[J].Va    = Va_s[t];
 *        DATA[J].mn    = mn0[t];
 *        DATA[J].i_rad = i_rad;
 *        DATA[J].phi   = i_rad;
 *        DATA[J].has_phi = 0;
 *        J++;
 *    }
 *
 *  ── General LGO (Z==1 branch, same pattern) ──────────────────────────
 *
 *    newton_solve_all(mn0, n_mn0, i_rad, phi, 0,
 *                     rA_s, Va_s, rP_s);
 *    // phi already computed from omega or passed directly
 *
 *  NOTE: the grid arrays (mn_k, Va_k) are still needed for the gnuplot
 *  CURVES (the dash-dot m/n line and dashed Va line). Only the marker
 *  identification step moves to Newton. find_hindex() and hindex[] are
 *  no longer needed.
 * ══════════════════════════════════════════════════════════════════════════ */


/* ══════════════════════════════════════════════════════════════════════════
 *  Standalone test  (gcc -O2 -DLGO_NEWTON_TEST -o test lgo_newton.c -lm)
 * ══════════════════════════════════════════════════════════════════════════ */
#ifdef LGO_NEWTON_TEST
int main(void)
{
    printf("── Newton solver test ────────────────────────────────────────────\n");
    printf("%-5s  %-8s  %-8s  %-10s  %-10s  %-10s  %-8s  %s\n",
           "i[°]","mn_tgt","mn_sol","H_A[km]","Va[km/s]","H_P[km]","e[-]","Case");
    printf("─────────────────────────────────────────────────────────────────\n");

    /* Test cases */
    struct { double i_deg; double phi_deg; double mn; int sp; } cases[] = {
        {  0.0,  0.0, 1.0, 1 },
        {  0.0,  0.0, 1.5, 1 },
        {  0.0,  0.0, 2.0, 1 },
        { 63.4, 63.4, 1.0, 1 },
        { 63.4, 63.4, 1.5, 1 },
        { 63.4, 63.4, 2.0, 1 },
        { 63.4, 53.0, 1.0, 0 },   /* general LGO */
        { 63.4, 30.0, 1.0, 0 },
    };
    int n = sizeof(cases)/sizeof(cases[0]);

    for (int t = 0; t < n; t++) {
        double i_r   = cases[t].i_deg   * PI / 180.0;
        double phi_r = cases[t].phi_deg * PI / 180.0;
        double rA, Va, rP;
        bool ok = newton_solve_rA(cases[t].mn, i_r, phi_r, cases[t].sp,
                                  &rA, &Va, &rP);
        if (ok) {
            double Ha = rA - RMEAN;
            double Hp = rP - RMEAN;
            double e  = (rA - rP)/(rA + rP);
            double cosg = cases[t].sp ? 1.0 : cos(i_r)/cos(phi_r);
            double mn_check = eval_mn_at_rA(rA, i_r, phi_r, cosg, cases[t].sp);
            printf("%-5.1f  %-8.4f  %-8.6f  %-10.1f  %-10.5f  %-10.1f  %-8.4f  %s\n",
                   cases[t].i_deg, cases[t].mn, mn_check,
                   Ha, Va, Hp, e,
                   cases[t].sp ? "special" : "general");
        } else {
            printf("%-5.1f  %-8.4f  %-8s  %-10s  %-10s  %-10s  %-8s  no solution\n",
                   cases[t].i_deg, cases[t].mn,
                   "—","—","—","—","—");
        }
    }

    printf("─────────────────────────────────────────────────────────────────\n");
    return 0;
}
#endif /* LGO_NEWTON_TEST */
