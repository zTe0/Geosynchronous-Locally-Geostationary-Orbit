/*
 * LGO_C2.c — Locally Geostationary Orbit (LGO) Analysis
 *
 * C equivalent of:  function DATA = LGO_C2(I, m_n_0, om, Z)
 *
 * Inputs (see main() for example call):
 *   I      [deg]   inclination array
 *   m_n_0  [-]     target m/n repetition factors (grid method)
 *   om     [deg]   argument of perigee (used when Z=1, general case)
 *   Z      [-]     0 = special LGO (phi=i), 1 = general LGO (phi from om)
 *
 * Output:
 *   LGOData[]  struct array — one entry per identified LGO orbit
 *              fields: H [km], Va [km/s], mn [-], i_rad [rad],
 *                      phi [rad], has_phi (0=special, 1=general)
 *
 * Build:  gcc -O2 -o lgo2 LGO_C2.c -lm (linux) 
 * or gcc LGO_C2win.c -o LGO_C2.exe (Windows from cmd)
 * Run:    ./lgo2 (linux)
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <sys/stat.h>

#ifdef _WIN32
    #include <windows.h>
    #include <direct.h>
    #define MKDIR(dir) _mkdir(dir)
#else
    #define MKDIR(dir) mkdir(dir, 0777)
#endif

/* ── Physical constants ─────────────────────────────────────────────────── */
#define MU    398600.4       /* [km^3/s^2] gravitational parameter            */
#define OME   7.2921159e-5   /* [s^-1]     Earth rotation rate                */
#define RMEAN 6371.0         /* [km]       Earth mean radius                  */
#define RE    6378.0         /* [km]       Earth equatorial radius             */
#define J2    0.0010827      /* [-]        second zonal harmonic               */
#define PI    3.14159265358979323846

/* ── Altitude grid ──────────────────────────────────────────────────────── */
#define N_H    60001   /* H = 0, 1, …, 60000 [km]                            */

/* ── Background fixed parameters ────────────────────────────────────────── */
#define N_HP    3      /* fixed perigee altitudes for background Va curves    */
#define N_I_BG  2      /* fixed inclinations for background Vi curves         */

/* ── Flat 2-D array accessors ───────────────────────────────────────────── */
#define VA_BG(h,k)  va_bg[(h)*N_HP   +(k)]
#define VI_BG(h,j)  vi_bg[(h)*N_I_BG +(j)]

/* ── gem12 colour palette ───────────────────────────────────────────────── */
static const char *gem12[12] = {
    "#0072BD","#D95319","#EDB120","#7E2F8E",
    "#77AC30","#4DBEEE","#A2142F","#FF1493",
    "#00CED1","#FF8C00","#228B22","#8B0000"
};

/* ── DATA output struct ─────────────────────────────────────────────────── */
typedef struct {
    double H;       /* [km]   apogee altitude at identified LGO orbit        */
    double Va;      /* [km/s] apogee velocity                                */
    double mn;      /* [-]    m/n target value used in grid method           */
    double i_rad;   /* [deg]  inclination                                    */
    double phi;     /* [deg]  latitude of apogee sub-satellite point         */
    double omega;      /* [deg]  argument of perigee                             */
    int    has_phi; /* 0 = special case (phi=i),  1 = general case           */
} LGOData;

static inline double sq(double x) { return x * x; }

/* ── compute_K_mn ───────────────────────────────────────────────────────── *
 *  Computes nodal correction factor K [-] and m/n repetition factor [-].   *
 *  Returns false if orbit parameters are invalid (a<=0 or 1-e^2<=0).       */
static bool compute_K_mn(double a,      /* [km]  semimajor axis              */
                         double e,      /* [-]   eccentricity                */
                         double i_rad,  /* [rad] inclination                 */
                         double *K_out, double *mn_out)
{
    double one_m_e2 = 1.0 - e*e;
    if (a <= 0.0 || one_m_e2 <= 0.0) return false;

    double ci   = cos(i_rad);
    double a2   = a*a;
    double a3h  = pow(a, 1.5);    /* [km^(3/2)] a^(3/2) */
    double a7h  = pow(a, 3.5);    /* [km^(7/2)] a^(7/2) */
    double e2sq = sq(one_m_e2);   /* [-] (1-e^2)^2       */
    double sqe  = sqrt(one_m_e2); /* [-] (1-e^2)^(1/2)   */

    double K = 1.0 - (1.5*J2*RE*RE / (a2*(2.0 - e*e))) *
        (  (10.0*sq(ci)-2.0)/4.0
         - (3.0*sq(ci)-1.0)/4.0 * sqe
         + sqrt(MU)*ci / (OME*a7h*e2sq) );

    double mn = 1.0 / (OME * a3h/sqrt(MU) * K);  /* [-] m/n repetition factor */

    if (!isfinite(K) || !isfinite(mn) || mn <= 0.0) return false;
    *K_out  = K;
    *mn_out = mn;
    return true;
}

/* ── find_hindex ────────────────────────────────────────────────────────── *
 *  Grid method: for each target mn0[t], scan mn_arr[0..n_h-1] and store    *
 *  in hindex[t] the index h with smallest |mn_arr[h] - mn0[t]|.            *
 *  Sets hindex[t]=-1 when no valid entry is found.                          *
 *  Implements: df = M_N - m_n_0;  [~,index] = min(abs(df),[],1)            */
static void find_hindex(const double *mn_arr, int n_h,
                        const double *mn0,    int n_mn0,
                        int *hindex)
{
    const double TOL = 0.005; // Strict tolerance for m/n matching
    for (int t = 0; t < n_mn0; t++) {
        double best = 1e30;
        hindex[t] = -1;
        for (int h = 0; h < n_h; h++) {
            double v = mn_arr[h];
            if (!isfinite(v) || v <= 0.0) continue;
            double d = fabs(v - mn0[t]);
            if (d < best) { best = d; hindex[t] = h; }
        }
         // If the best available H results in an m/n error > TOL, reject it.
        // This stops points from "sticking" to the circular curve incorrectly.
        if (hindex[t] != -1 && fabs(mn_arr[hindex[t]] - mn0[t]) > TOL) {
            hindex[t] = -1;
        } else {
            hindex[t] = hindex[t];
        }
    }
}

/* ══════════════════════════════════════════════════════════════════════════ *
 *  LGO_C2 — main computation function                                       *
 *  Returns heap-allocated DATA array; caller must free().                   *
 * ══════════════════════════════════════════════════════════════════════════ */
LGOData *LGO_C2(const double *I,    int n_I,    /* [deg] inclinations        */
                const double *mn0,  int n_mn0,  /* [-]   m/n targets         */
                const double *om,   int n_om,   /* [deg] argument of perigee */
                int    Z,                       /* 0=special, 1=general      */
                int   *n_data_out)
{
    /* ── Altitude grid ──────────────────────────────────────────────────── */
    double *H  = malloc(N_H * sizeof(double));  /* [km]   altitude           */
    double *r  = malloc(N_H * sizeof(double));  /* [km]   geocentric radius  */
    double *vc = malloc(N_H * sizeof(double));  /* [km/s] circular velocity  */

    if (!H || !r || !vc) { fprintf(stderr,"malloc failed\n"); exit(1); }

    for (int h = 0; h < N_H; h++) {
        H[h]  = (double)h;
        r[h]  = H[h] + RMEAN;
        vc[h] = sqrt(MU / r[h]);
    }

    /* ── Background Va: apogee velocity for fixed perigee altitudes ─────── */
    const double Hp_bg[N_HP]    = {500.0, 5000.0, 15000.0};  /* [km] */
    double *va_bg = malloc(N_H * N_HP   * sizeof(double));    /* [km/s] */

    for (int k = 0; k < N_HP; k++) {
        double rp = Hp_bg[k] + RMEAN;              /* [km] perigee geocentric radius */
        for (int h = 0; h < N_H; h++) {
            double rA = r[h];
            VA_BG(h,k) = sqrt(2.0*MU*rp / (rA*(rA+rp)));
        }
    }

    /* ── Background Vi: geocentric-tip velocity for fixed inclinations ───── */
    const double i_bg_deg[N_I_BG] = {0.0, 63.4};             /* [deg] */
    double *vi_bg = malloc(N_H * N_I_BG * sizeof(double));    /* [km/s] */

    for (int j = 0; j < N_I_BG; j++) {
        double ci = cos(i_bg_deg[j] * PI / 180.0);
        for (int h = 0; h < N_H; h++)
            VI_BG(h,j) = OME * r[h] * ci;
    }

    /* ── Write background gnuplot data files ────────────────────────────── */
    FILE *f;

    f = fopen("data/lgo2_vc.dat","w");
    for (int h=0; h<N_H; h++) fprintf(f,"%.0f %.8f\n", H[h], vc[h]);
    fclose(f);

    for (int k=0; k<N_HP; k++) {
        char fn[64]; sprintf(fn,"data/lgo2_va_bg%d.dat",k);
        f = fopen(fn,"w");
        for (int h=0; h<N_H; h++)
            if (vc[h] >= VA_BG(h,k))
                fprintf(f,"%.0f %.8f\n", H[h], VA_BG(h,k));
        fclose(f);
    }

    for (int j=0; j<N_I_BG; j++) {
        char fn[64]; sprintf(fn,"data/lgo2_vi_bg%d.dat",j);
        f = fopen(fn,"w");
        for (int h=0; h<N_H; h++) {
            bool ok = (j==0) ? (H[h]>=15000.0 && H[h]<=45000.0)
                             : (H[h]>=19000.0);
            if (ok) fprintf(f,"%.0f %.8f\n", H[h], VI_BG(h,j));
        }
        fclose(f);
    }

    /* ── DATA output buffer ─────────────────────────────────────────────── */
    int max_data = n_I * n_mn0 * (n_om > 0 ? n_om+1 : 1) + 32;
    LGOData *DATA = calloc(max_data, sizeof(LGOData));
    int J = 0;   /* DATA entry counter */

    /* ── gnuplot plot command — built incrementally ─────────────────────── */
    int   pbufsz = 16384;
    char *plot_cmd = malloc(pbufsz);
    int   plen = 0;

    plen += snprintf(plot_cmd+plen, pbufsz-plen,
        "plot \\\n"
        "  'data/lgo2_vc.dat'     u 1:2 axes x1y1 w l lt 1 t 'V_c(H_c)', \\\n"
        "  'data/lgo2_va_bg0.dat' u 1:2 axes x1y1 w l lt 2 t 'V_A(H_P=500 km)', \\\n"
        "  'data/lgo2_va_bg1.dat' u 1:2 axes x1y1 w l lt 3 t 'V_A(H_P=5000 km)', \\\n"
        "  'data/lgo2_va_bg2.dat' u 1:2 axes x1y1 w l lt 4 t 'V_A(H_P=15000 km)', \\\n"
        "  'data/lgo2_vi_bg0.dat' u 1:2 axes x1y1 w l lt 5 t 'V_i(i=0 deg)', \\\n"
        "  'data/lgo2_vi_bg1.dat' u 1:2 axes x1y1 w l lt 6 t 'V_i(i=63.4 deg)'");

    int gp_color = 0;   /* gem12 index for LGO-specific curves */

    /* ── Main inclination loop ──────────────────────────────────────────── */
    for (int k = 0; k < n_I; k++) {

        double i_rad = I[k] * PI / 180.0;  /* [rad] inclination */
        double ci    = cos(i_rad);
        double ci2   = ci*ci;

        if (Z == 0) {
            /* ── Special LGO: phi = i ─────────────────────────────────── */
            double *rP   = malloc(N_H * sizeof(double));  /* [km]   special LGO perigee */
            double *mn_k = malloc(N_H * sizeof(double));  /* [-]    m/n factor          */
            double *Va_k = malloc(N_H * sizeof(double));  /* [km/s] apogee velocity     */
            bool   *vld  = malloc(N_H * sizeof(bool));

            double omega = 270;

            for (int h = 0; h < N_H; h++) {
                double rAh = r[h];
                double rA3 = rAh*rAh*rAh;
                double den = 2.0*MU - sq(OME)*rA3*ci2;

                vld[h]=false; rP[h]=mn_k[h]=Va_k[h]=NAN;
                if (den <= 0.0 || fabs(den) < 1e-6) continue;

                double rp = sq(OME)*rA3*rAh*ci2 / den;   /* [km] special LGO perigee */
                if (rp <= 0.0 || rp > rAh) continue;

                double a = (rAh+rp)/2.0;                  /* [km] semimajor axis      */
                double e = (rAh-rp)/(rAh+rp);             /* [-]  eccentricity        */
                double K, mn_val;
                if (!compute_K_mn(a, e, i_rad, &K, &mn_val)) continue;

                rP[h]   = rp;
                mn_k[h] = mn_val;
                Va_k[h] = sqrt(2.0*MU*rp / (rAh*(rAh+rp)));  /* [km/s] */
                vld[h]  = true;
            }

            /* Grid method: find H closest to each mn0 target */
            int *hindex = malloc(n_mn0 * sizeof(int));
            find_hindex(mn_k, N_H, mn0, n_mn0, hindex);

            /* Write data files */
            char fn_mrk[64], fn_mn[64], fn_va[64];
            sprintf(fn_mrk,"data/lgo2_mrk_sp%d.dat", k);
            sprintf(fn_mn, "data/lgo2_mn_sp%d.dat",  k);
            sprintf(fn_va, "data/lgo2_va_sp%d.dat",  k);

            FILE *fmrk=fopen(fn_mrk,"w"), *fmn=fopen(fn_mn,"w"), *fva=fopen(fn_va,"w");
            for (int h=0; h<N_H; h++) {
                if (!vld[h]) continue;
                fprintf(fmn,"%.0f %.8f\n", H[h], mn_k[h]);
                fprintf(fva,"%.0f %.8f\n", H[h], Va_k[h]);
            }
            for (int t=0; t<n_mn0; t++) {
                int h = hindex[t];
                if (h>=0 && isfinite(Va_k[h]))
                    fprintf(fmrk,"%.0f %.8f %.2f\n", H[h], Va_k[h], mn0[t]);
            }
            fclose(fmrk); fclose(fmn); fclose(fva);

            /* Append gnuplot commands:
             *   markers  → left  axis (x1y1)
             *   m/n, Va  → right axis (x1y2)  as in MATLAB yyaxis right  */
            int c = gp_color % 12;
            plen += snprintf(plot_cmd+plen, pbufsz-plen,
                ", \\\n"
                "  '%s' u 1:2 axes x1y1 w p lc rgb '%s' pt 7 ps 1.4"
                "  t 'V_A markers (i=%.4g deg / m_n targets)', \\\n"
                "  '%s' u 1:2 axes x1y2 w l lc rgb '%s' dt 4 lw 2"
                "  t 'm/n (i=%.4g deg / special LGO)', \\\n"
                "  '%s' u 1:2 axes x1y2 w l lc rgb '%s' dt 2 lw 2"
                "  t 'V_A (i=%.4g deg / special LGO)'",
                fn_mrk, gem12[c], I[k],
                fn_mn,  gem12[c], I[k],
                fn_va,  gem12[c], I[k]);
            gp_color++;

            /* Populate DATA */
            for (int t=0; t<n_mn0; t++) {
                int h = hindex[t];
                if (h<0 || !isfinite(Va_k[h])) continue;
                DATA[J].H       = H[h];
                DATA[J].Va      = Va_k[h];
                DATA[J].mn      = mn0[t];
                DATA[J].i_rad   = i_rad;
                DATA[J].phi     = i_rad;   /* phi = i in special case */
                DATA[J].omega   = omega;
                DATA[J].has_phi = 0;
                J++;
            }

            free(rP); free(mn_k); free(Va_k); free(vld); free(hindex);

        } else {
            /* ── General LGO: phi computed from om, Z=1 ──────────────── */
            for (int o = 0; o < n_om; o++) {
                // If inclination is 0, only calculate for the first omega to avoid duplicates
                if (I[k] == 0 && o > 0) continue;

                double phi = asin(sin((om[o] * PI / 180.0) + PI) * sin(i_rad));
                double omega = om[o];

                if (fabs(phi) > fabs(i_rad)) continue; 

                double cp  = cos(phi);
                double sp  = sin(phi);
                double tp  = sp / cp; 
                double cosg = sqrt(1.0 - sq(sin(i_rad)) + sq(tp) * sq(cos(i_rad)));

                double *mn_g = malloc(N_H * sizeof(double));
                double *Va_g = malloc(N_H * sizeof(double));
                bool   *vld  = malloc(N_H * sizeof(bool));

                for (int h = 0; h < N_H; h++) {
                    double rAh = r[h];
                    double rA3 = pow(rAh, 3);
                    double den = 2.0 * MU * cosg * cosg - sq(OME) * rA3 * sq(cp);

                    vld[h] = false; 
                    mn_g[h] = Va_g[h] = NAN;
                    
                    if (den <= 0.0 || fabs(den) < 1e-6) continue;

                    double rp = (sq(OME) * pow(rAh, 4) * sq(cp)) / den;
                    if (rp <= 0.0 || rp > rAh) continue;

                    double a = (rAh + rp) / 2.0;
                    double e = (rAh - rp) / (rAh + rp);
                    double K, mn_val;
                    if (!compute_K_mn(a, e, i_rad, &K, &mn_val)) continue;

                    mn_g[h] = mn_val;
                    Va_g[h] = sqrt((2.0 * MU * rp) / (rAh * (rAh + rp)));
                    vld[h]  = true;
                }
                
                int *hindex = malloc(n_mn0 * sizeof(int));
                find_hindex(mn_g, N_H, mn0, n_mn0, hindex);

                char fn_mrk[64], fn_mn[64], fn_va[64];
                // FIX: Use 'o' instead of 'jj'
                sprintf(fn_mrk,"data/lgo2_mrk_gen_%d_%d.dat", k, o);
                sprintf(fn_mn, "data/lgo2_mn_gen_%d_%d.dat",  k, o);
                sprintf(fn_va, "data/lgo2_va_gen_%d_%d.dat",  k, o);

                FILE *fmrk=fopen(fn_mrk,"w"), *fmn=fopen(fn_mn,"w"), *fva=fopen(fn_va,"w");
                for (int h=0; h<N_H; h++) {
                    if (!vld[h]) continue;
                    fprintf(fmn,"%.0f %.8f\n", H[h], mn_g[h]);
                    fprintf(fva,"%.0f %.8f\n", H[h], Va_g[h]);
                }
                for (int t=0; t<n_mn0; t++) {
                    int h = hindex[t];
                    if (h>=0 && isfinite(Va_g[h]))
                        fprintf(fmrk,"%.0f %.8f %.2f\n", H[h], Va_g[h], mn0[t]);
                }
                fclose(fmrk); fclose(fmn); fclose(fva);

                int c = gp_color % 12;
                // FIX: Use 'o' instead of 'jj' inside snprintf for plot_cmd
                plen += snprintf(plot_cmd+plen, pbufsz-plen,
                    ", \\\n"
                    "  '%s' u 1:2 axes x1y1 w p lc rgb '%s' pt 7 ps 1.4"
                    "  t 'V_A markers (i=%.4g deg / phi=%.4g deg)', \\\n"
                    "  '%s' u 1:2 axes x1y2 w l lc rgb '%s' dt 4 lw 2"
                    "  t 'm/n (i=%.4g deg / phi=%.4g deg)', \\\n"
                    "  '%s' u 1:2 axes x1y2 w l lc rgb '%s' dt 2 lw 2"
                    "  t 'V_A (i=%.4g deg / phi=%.4g deg)'",
                    fn_mrk, gem12[c], I[k], phi*180.0/PI,
                    fn_mn,  gem12[c], I[k], phi*180.0/PI,
                    fn_va,  gem12[c], I[k], phi*180.0/PI);
                gp_color++;

                /* Populate DATA */
                for (int t=0; t<n_mn0; t++) {
                    int h = hindex[t];
                    if (h<0 || !isfinite(Va_g[h])) continue;
                    DATA[J].H       = H[h];
                    DATA[J].Va      = Va_g[h];
                    DATA[J].mn      = mn0[t];
                    DATA[J].i_rad   = i_rad;
                    DATA[J].phi     = phi;
                    DATA[J].omega   = omega;
                    DATA[J].has_phi = 1;
                    J++;
                }

                free(mn_g); free(Va_g); free(vld); free(hindex);
            }
        }
    }

    /* ── Gnuplot ─────────────────────────────────────────────────────────── */
    FILE *gp = popen("gnuplot -persistent","w");
    if (!gp) { fprintf(stderr,"Cannot open gnuplot\n"); goto done; }

    fprintf(gp,"set terminal qt size 800,750 enhanced font 'Sans,10'\n");
    fprintf(gp,"set title 'Orbital Velocity and Repetition Factor vs Altitude'\n");
    fprintf(gp,"set xlabel 'H_A [km]'\n");
    fprintf(gp,"set ylabel 'V [km/s]'\n");
    fprintf(gp,"set y2label 'm/n [-]'\n");
    fprintf(gp,"set ytics  nomirror\n");
    fprintf(gp,"set y2tics\n");
    fprintf(gp,"set yrange  [0:8]\n");
    fprintf(gp,"set y2range [0:8]\n");
    fprintf(gp,"set xrange  [0:60000]\n");
    fprintf(gp,"set grid\n");
    fprintf(gp,"set key outside right top font ',8' spacing 1.2\n");
    fprintf(gp,"set format x '%%g'\n");
    for (int c=0; c<12; c++)
        fprintf(gp,"set linetype %d lc rgb '%s' lw 2\n", c+1, gem12[c]);

    /* Close plot command and send */
    plen += snprintf(plot_cmd+plen, pbufsz-plen, "\n");
    fprintf(gp,"%s", plot_cmd);
    pclose(gp);

done:
    /* ── Print DATA table to stdout ─────────────────────────────────────── */
    printf("\n── DATA output ───────────────────────────────────────────────────────────────────\n");
    printf("%-4s  %-10s  %-10s  %-7s  %-10s  %-10s  %-10s  %s\n",
           "Idx","H [km]","Va [km/s]","m/n","i [deg]","phi [deg]","om [deg]","Case");
    for (int j=0; j<J; j++)
        printf("%-4d  %-10.1f  %-10.5f  %-7.3f  %-10.4f  %-10.4f  %-10.1f  %s\n",
               j+1,
               DATA[j].H,
               DATA[j].Va,
               DATA[j].mn,
               DATA[j].i_rad * 180.0/PI,
               DATA[j].phi   * 180.0/PI,
               DATA[j].omega,
               DATA[j].has_phi ? "general" : "special");
    printf("──────────────────────────────────────────────────────────────────────────────────\n");
    printf("Total DATA entries: %d\n\n", J);

    FILE *f_table = fopen("data/table_results.txt", "w");
    if (f_table == NULL) {
        printf("Error: Could not open or create the file. Make sure the folder exists!\n");
        return NULL;
    }

    fprintf(f_table,"\n── DATA output ─────────────────────────────────────────────────────\n");
    fprintf(f_table,"%-4s  %-10s  %-10s  %-7s  %-10s  %-10s  %-10s  %s\n",
           "Idx","H [km]","Va [km/s]","m/n","i [deg]","phi [deg]","om [deg]","Case");
    for (int j=0; j<J; j++)
        fprintf(f_table,"%-4d  %-10.1f  %-10.5f  %-7.3f  %-10.4f  %-10.4f  %-10.1f  %s\n",
               j+1,
               DATA[j].H,
               DATA[j].Va,
               DATA[j].mn,
               DATA[j].i_rad * 180.0/PI,
               DATA[j].phi   * 180.0/PI,
               DATA[j].omega,
               DATA[j].has_phi ? "general" : "special");
    fprintf(f_table,"─────────────────────────────────────────────────────────────────\n");
    fprintf(f_table,"Total DATA entries: %d\n\n", J);
    fclose(f_table);

    *n_data_out = J;
    free(H); free(r); free(vc);
    free(va_bg); free(vi_bg);
    free(plot_cmd);
    return DATA;
}

/* ══════════════════════════════════════════════════════════════════════════ *
 *  main — example call                                                       *
 * ══════════════════════════════════════════════════════════════════════════ */
int main(void)
{
    double I[]   = {0.0, 63.4};                       /* [deg] inclinations   */
    double mn0[] = {1.0, 1.5, 2.0};                   /* [-]   m/n targets    */
    double om[]  = {135.0, 90.0};                     /* [deg] arg. perigee   */
    int    Z     = 1;                                 /* 0=special, 1=general */

    int n_I   = sizeof(I)   / sizeof(I[0]);
    int n_mn0 = sizeof(mn0) / sizeof(mn0[0]);
    int n_om  = sizeof(om)  / sizeof(om[0]);

    int n_data;

    // Set console output to UTF-8 for Windows
    #ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    #endif

    const char *dir = "data";
    // Create the directory (ignores if it already exists)
    MKDIR(dir);

    LGOData *DATA = LGO_C2(I, n_I, mn0, n_mn0, om, n_om, Z, &n_data);

    printf("Press Enter to exit...");
    getchar();
    free(DATA);
    
    return 0;
}
