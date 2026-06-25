/*
 * LGO_GUI.c — GTK3 GUI front-end for Geosynchronous LGO Analysis
 *
 * Self-contained: includes all computation from LGO_C2 plus the GUI layer.
 *
 * ── Build (Linux) ───────────────────────────────────────────────────────
 *   sudo apt install libgtk-3-dev gnuplot
 *   sudo add-apt-repository ppa:kisak/kisak-mesa
 *   sudo apt update && sudo apt upgrade -y
 *   gcc -O2 -Wno-deprecated-declarations $(pkg-config --cflags gtk+-3.0) -o lgo_gui LGO_GUI.c $(pkg-config --libs gtk+-3.0) -lm
 *
 * ── Build (Windows / MSYS2 MinGW64 shell) ───────────────────────────────
 *   pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-pkg-config mingw-w64-x86_64-gtk3 mingw-w64-x86_64-gnuplot
 *   gcc -O2 $(pkg-config --cflags gtk+-3.0) -o lgo_gui.exe LGO_GUI.c $(pkg-config --libs gtk+-3.0) -lm -mwindows
 *   cmd terminal:
 *   bash -c "gcc -O2 -Wno-deprecated-declarations $(pkg-config --cflags gtk+-3.0) -o lgo_gui.exe LGO_GUI.c $(pkg-config --libs gtk+-3.0) -lm -mwindows"
 *
 * ── Run ─────────────────────────────────────────────────────────────────
 *   ./lgo_gui          (Linux)
 *   lgo_gui.exe        (Windows)
 */

#include <gtk/gtk.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#ifdef _WIN32
  #include <windows.h>
  #include <direct.h>
  #define MKDIR(p)  _mkdir(p)
  #define POPEN     _popen
  #define PCLOSE    _pclose
#else
  #include <sys/stat.h>
  #define MKDIR(p)  mkdir(p, 0777)
  #define POPEN     popen
  #define PCLOSE    pclose
#endif

/* ══════════════════════════════════════════════════════════════════════════
 *  PHYSICS / COMPUTATION
 * ══════════════════════════════════════════════════════════════════════════ */

#define MU    398600.4       /* [km^3/s^2] gravitational parameter            */
#define OME   7.2921159e-5   /* [s^-1]     Earth rotation rate                */
#define RMEAN 6371.0         /* [km]       Earth mean radius                  */
#define RE    6378.0         /* [km]       Earth equatorial radius             */
#define J2    0.0010827      /* [-]        second zonal harmonic               */
#define PI    3.14159265358979323846

#define N_H    60001   /* altitude grid: H = 0, 1, …, 60000 [km]             */
#define N_HP   3       /* fixed perigee altitudes for background Va curves    */
#define N_I_BG 2       /* fixed inclinations for background Vi curves         */

#define VA_BG(h,k)  va_bg[(h)*N_HP   +(k)]
#define VI_BG(h,j)  vi_bg[(h)*N_I_BG +(j)]

static const char *gem12[12] = {
    "#0072BD","#D95319","#EDB120","#7E2F8E",
    "#77AC30","#4DBEEE","#A2142F","#FF1493",
    "#00CED1","#FF8C00","#228B22","#8B0000"
};

typedef struct {
    double H;       /* [km]   apogee altitude                                */
    double Va;      /* [km/s] apogee velocity                                */
    double mn;      /* [-]    m/n repetition factor (target value)           */
    double i_rad;   /* [rad]  inclination                                    */
    double phi;     /* [rad]  latitude of apogee sub-satellite point         */
    double omega;   /* [deg]  argument of perigee (NAN if phi-direct input)  */
    double a;
    double e;
    int    has_phi; /* 0 = special case (phi=i),  1 = general case           */
} LGOData;

static inline double sq(double x) { return x * x; }

/* ── compute_K_mn ───────────────────────────────────────────────────────── */
static bool compute_K_mn(double a, double e, double i_rad,
                          double *K_out, double *mn_out)
{
    double one_m_e2 = 1.0 - e*e;
    if (a <= 0.0 || one_m_e2 <= 0.0) return false;

    double ci   = cos(i_rad);
    double a2   = a*a;
    double a3h  = pow(a, 1.5);    /* [km^(3/2)] a^(3/2) */
    double a7h  = pow(a, 3.5);    /* [km^(7/2)] a^(7/2) */
    double e2sq = sq(one_m_e2);
    double sqe  = sqrt(one_m_e2);

    double K = 1.0 - (1.5*J2*RE*RE / (a2*e2sq)) *
        (  (10.0*sq(ci)-2.0)/4.0
         - (3.0*sq(ci)-1.0)/4.0 * sqe
         + sqrt(MU)*ci / (OME*a7h*e2sq) );

    double mn = 1.0 / (OME * a3h/sqrt(MU) * K);  /* [-] m/n repetition factor */

    if (!isfinite(K) || !isfinite(mn) || mn <= 0.0) return false;
    *K_out  = K;
    *mn_out = mn;
    return true;
}

/* ── find_hindex ────────────────────────────────────────────────────────── */
static void find_hindex(const double *mn_arr, int n_h,
                        const double *mn0, int n_mn0, int *hindex)
{
    const double TOL = 0.005;   /* strict tolerance for m/n matching         */
    for (int t = 0; t < n_mn0; t++) {
        double best = 1e30;
        hindex[t] = -1;
        for (int h = 0; h < n_h; h++) {
            double v = mn_arr[h];
            if (!isfinite(v) || v <= 0.0) continue;
            double d = fabs(v - mn0[t]);
            if (d < best) { best = d; hindex[t] = h; }
        }
        if (hindex[t] != -1 && fabs(mn_arr[hindex[t]] - mn0[t]) > TOL)
            hindex[t] = -1;
    }
}

/* ── LGO_C2 ─────────────────────────────────────────────────────────────── *
 *  angles[]: omega [deg] when use_om=1, phi [deg] when use_om=0             *
 *  Returns allocated DATA array (caller must free).                         */
static LGOData *LGO_C2(const double *I,      int n_I,
                        const double *mn0,    int n_mn0,
                        const double *angles, int n_angles,
                        int Z, int use_om,
                        int *n_data_out)
{
    /* ── Altitude grid ──────────────────────────────────────────────────── */
    double *H  = malloc(N_H * sizeof(double));  /* [km]   altitude           */
    double *r  = malloc(N_H * sizeof(double));  /* [km]   geocentric radius  */
    double *vc = malloc(N_H * sizeof(double));  /* [km/s] circular velocity  */

    for (int h = 0; h < N_H; h++) {
        H[h]  = (double)h;
        r[h]  = H[h] + RMEAN;
        vc[h] = sqrt(MU / r[h]);
    }

    /* ── Background Va for fixed perigee altitudes ──────────────────────── */
    const double Hp_bg[N_HP] = {500.0, 5000.0, 15000.0};     /* [km] */
    double *va_bg = malloc(N_H * N_HP   * sizeof(double));
    for (int k = 0; k < N_HP; k++) {
        double rp = Hp_bg[k] + RMEAN;
        for (int h = 0; h < N_H; h++) {
            double rA = r[h];
            VA_BG(h,k) = sqrt(2.0*MU*rp / (rA*(rA+rp)));
        }
    }

    /* ── Background Vi for fixed inclinations ───────────────────────────── */
    const double i_bg_deg[N_I_BG] = {0.0, 63.4};             /* [deg] */
    double *vi_bg = malloc(N_H * N_I_BG * sizeof(double));
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
    int max_data = n_I * n_mn0 * (n_angles > 0 ? n_angles+1 : 1) + 32;
    LGOData *DATA = calloc(max_data, sizeof(LGOData));
    int J = 0;

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

    int gp_color = 0;

    /* ── Main inclination loop ──────────────────────────────────────────── */
    for (int k = 0; k < n_I; k++) {

        double i_rad = I[k] * PI / 180.0;
        double ci    = cos(i_rad);
        double ci2   = ci*ci;
        double *A = malloc(N_H * sizeof(double));  /* [km] Semi-major axis   */
        double *E = malloc(N_H * sizeof(double));  /* [-] Eccentricity       */

        if (Z == 0) {
            /* ── Special LGO: phi = i ─────────────────────────────────── */
            double *rP   = malloc(N_H * sizeof(double));  /* [km]   perigee   */
            double *mn_k = malloc(N_H * sizeof(double));  /* [-]    m/n       */
            double *Va_k = malloc(N_H * sizeof(double));  /* [km/s] Va        */
            bool   *vld  = malloc(N_H * sizeof(bool));

            double omega = 270;

            for (int h = 0; h < N_H; h++) {
                double rAh = r[h];
                double rA3 = rAh*rAh*rAh;
                double den = 2.0*MU - sq(OME)*rA3*ci2;
                vld[h]=false; rP[h]=mn_k[h]=Va_k[h]=NAN;
                if (den <= 0.0 || fabs(den) < 1e-6) continue;
                double rp = sq(OME)*rA3*rAh*ci2 / den;
                if (rp <= 0.0 || rp > rAh) continue;
                double a = (rAh+rp)/2.0;
                A[h] = a;
                double e = (rAh-rp)/(rAh+rp);
                E[h] = e;
                double K, mn_val;
                if (!compute_K_mn(a, e, i_rad, &K, &mn_val)) continue;
                rP[h]=rp; mn_k[h]=mn_val;
                Va_k[h] = sqrt(2.0*MU*rp / (rAh*(rAh+rp)));
                vld[h]=true;
            }

            int *hindex = malloc(n_mn0 * sizeof(int));
            find_hindex(mn_k, N_H, mn0, n_mn0, hindex);

            char fn_mrk[64], fn_mn[64], fn_va[64];
            sprintf(fn_mrk,"data/lgo2_mrk_sp%d.dat",k);
            sprintf(fn_mn, "data/lgo2_mn_sp%d.dat", k);
            sprintf(fn_va, "data/lgo2_va_sp%d.dat", k);

            FILE *fmrk=fopen(fn_mrk,"w"), *fmn=fopen(fn_mn,"w"), *fva=fopen(fn_va,"w");
            for (int h=0; h<N_H; h++) {
                if (!vld[h]) continue;
                fprintf(fmn,"%.0f %.8f\n", H[h], mn_k[h]);
                fprintf(fva,"%.0f %.8f\n", H[h], Va_k[h]);
            }
            for (int t=0; t<n_mn0; t++) {
                int h=hindex[t];
                if (h>=0 && isfinite(Va_k[h]))
                    fprintf(fmrk,"%.0f %.8f %.2f\n", H[h], Va_k[h], mn0[t]);
            }
            fclose(fmrk); fclose(fmn); fclose(fva);

            int c = gp_color % 12;
            plen += snprintf(plot_cmd+plen, pbufsz-plen,
                ", \\\n"
                "  '%s' u 1:2 axes x1y1 w p lc rgb '%s' pt 7 ps 1.4"
                "  t 'V_A markers (i=%.4g deg)', \\\n"
                "  '%s' u 1:2 axes x1y2 w l lc rgb '%s' dt 4 lw 2"
                "  t 'm/n (i=%.4g deg / special LGO)', \\\n"
                "  '%s' u 1:2 axes x1y2 w l lc rgb '%s' dt 2 lw 2"
                "  t 'V_A (i=%.4g deg / special LGO)'",
                fn_mrk,gem12[c],I[k], fn_mn,gem12[c],I[k], fn_va,gem12[c],I[k]);
            gp_color++;

            for (int t=0; t<n_mn0; t++) {
                int h=hindex[t];
                if (h<0 || !isfinite(Va_k[h])) continue;
                DATA[J].H=H[h]; 
                DATA[J].Va=Va_k[h]; 
                DATA[J].mn=mn0[t];
                DATA[J].i_rad=i_rad; 
                DATA[J].phi=i_rad;
                DATA[J].omega=omega;
                DATA[J].a = A[h];
                DATA[J].e = E[h]; 
                DATA[J].has_phi=0;
                J++;
            }
            free(rP); free(mn_k); free(Va_k); free(vld); free(hindex);

        } else {
            /* ── General LGO: phi from om or phi direct ───────────────── */
            for (int o = 0; o < n_angles; o++) {
                if (I[k] == 0.0 && o > 0) continue;  /* skip duplicates at i=0 */

                double phi, omega_deg;
                if (use_om) {
                    /* angles[] are omega [deg] → compute phi */
                    omega_deg = angles[o];
                    phi = asin(sin(angles[o]*PI/180.0 + PI) * sin(i_rad));
                } else {
                    /* angles[] are phi [deg] → use directly */
                    phi = angles[o] * PI / 180.0;
                    omega_deg = (asin(sin(phi)/sin(i_rad)) - PI) * 180/PI;
                }

                if (fabs(phi) > fabs(i_rad)) continue;  /* |phi| > |i| invalid */

                double cp   = cos(phi);
                double sp   = sin(phi);
                double tp   = (fabs(cp) > 1e-9) ? sp/cp : 0.0;  /* tan(phi) */
                double cosg = sqrt(fmax(0.0, 1.0 - sq(sin(i_rad)) + sq(tp)*ci2));

                double *mn_g = malloc(N_H * sizeof(double));
                double *Va_g = malloc(N_H * sizeof(double));
                bool   *vld  = malloc(N_H * sizeof(bool));

                for (int h=0; h<N_H; h++) {
                    double rAh = r[h];
                    double rA3 = rAh*rAh*rAh;
                    double den = 2.0*MU*cosg*cosg - sq(OME)*rA3*sq(cp);
                    vld[h]=false; mn_g[h]=Va_g[h]=NAN;
                    if (den <= 0.0 || fabs(den) < 1e-6) continue;
                    double rp = sq(OME)*rA3*rAh*sq(cp) / den;
                    if (rp <= 0.0 || rp > rAh) continue;
                    double a=(rAh+rp)/2.0, e=(rAh-rp)/(rAh+rp);
                    A[h] = a;
                    E[h] = e;
                    double K, mn_val;
                    if (!compute_K_mn(a, e, i_rad, &K, &mn_val)) continue;
                    mn_g[h]=mn_val;
                    Va_g[h]=sqrt(2.0*MU*rp / (rAh*(rAh+rp)));
                    vld[h]=true;
                }

                int *hindex = malloc(n_mn0 * sizeof(int));
                find_hindex(mn_g, N_H, mn0, n_mn0, hindex);

                char fn_mrk[64], fn_mn[64], fn_va[64];
                sprintf(fn_mrk,"data/lgo2_mrk_gen_%d_%d.dat",k,o);
                sprintf(fn_mn, "data/lgo2_mn_gen_%d_%d.dat", k,o);
                sprintf(fn_va, "data/lgo2_va_gen_%d_%d.dat", k,o);

                FILE *fmrk=fopen(fn_mrk,"w"), *fmn=fopen(fn_mn,"w"), *fva=fopen(fn_va,"w");
                for (int h=0; h<N_H; h++) {
                    if (!vld[h]) continue;
                    fprintf(fmn,"%.0f %.8f\n", H[h], mn_g[h]);
                    fprintf(fva,"%.0f %.8f\n", H[h], Va_g[h]);
                }
                for (int t=0; t<n_mn0; t++) {
                    int h=hindex[t];
                    if (h>=0 && isfinite(Va_g[h]))
                        fprintf(fmrk,"%.0f %.8f %.2f\n", H[h], Va_g[h], mn0[t]);
                }
                fclose(fmrk); fclose(fmn); fclose(fva);

                int c = gp_color % 12;
                plen += snprintf(plot_cmd+plen, pbufsz-plen,
                    ", \\\n"
                    "  '%s' u 1:2 axes x1y1 w p lc rgb '%s' pt 7 ps 1.4"
                    "  t 'V_A markers (i=%.4g / phi=%.4g deg)', \\\n"
                    "  '%s' u 1:2 axes x1y2 w l lc rgb '%s' dt 4 lw 2"
                    "  t 'm/n (i=%.4g / phi=%.4g deg)', \\\n"
                    "  '%s' u 1:2 axes x1y2 w l lc rgb '%s' dt 2 lw 2"
                    "  t 'V_A (i=%.4g / phi=%.4g deg)'",
                    fn_mrk,gem12[c],I[k],phi*180.0/PI,
                    fn_mn, gem12[c],I[k],phi*180.0/PI,
                    fn_va, gem12[c],I[k],phi*180.0/PI);
                gp_color++;

                for (int t=0; t<n_mn0; t++) {
                    int h=hindex[t];
                    if (h<0 || !isfinite(Va_g[h])) continue;
                    DATA[J].H=H[h]; 
                    DATA[J].Va=Va_g[h]; 
                    DATA[J].mn=mn0[t];
                    DATA[J].i_rad=i_rad; 
                    DATA[J].phi=phi;
                    DATA[J].omega=omega_deg;
                    DATA[J].a = A[h];
                    DATA[J].e = E[h]; 
                    DATA[J].has_phi=1;
                    J++;
                }
                free(mn_g); free(Va_g); free(vld); free(hindex);
            }
        }
    }

    /* ── Launch gnuplot ─────────────────────────────────────────────────── */
    FILE *gp = POPEN("gnuplot -persistent","w");
    if (gp) {
        fprintf(gp,"set terminal qt size 800,750 enhanced font 'Sans,10'\n");
        fprintf(gp,"set title 'Orbital Velocity and Repetition Factor vs Altitude'\n");
        fprintf(gp,"set xlabel 'H_A [km]'\n");
        fprintf(gp,"set ylabel 'V [km/s]'\n");
        fprintf(gp,"set y2label 'm/n [-]'\n");
        fprintf(gp,"set ytics nomirror\n");
        fprintf(gp,"set y2tics\n");
        fprintf(gp,"set yrange  [0:8]\n");
        fprintf(gp,"set y2range [0:8]\n");
        fprintf(gp,"set xrange  [0:60000]\n");
        fprintf(gp,"set grid\n");
        fprintf(gp,"set key outside right top font ',8' spacing 1.2\n");
        fprintf(gp,"set format x '%%g'\n");
        for (int c=0; c<12; c++)
            fprintf(gp,"set linetype %d lc rgb '%s' lw 2\n", c+1, gem12[c]);
        plen += snprintf(plot_cmd+plen, pbufsz-plen, "\n");
        fprintf(gp,"%s", plot_cmd);
        PCLOSE(gp);
    }

    // /* ── Print DATA table to stdout ─────────────────────────────────────── */
    // printf("\n── DATA output ───────────────────────────────────────────────────────────────────\n");
    // printf("%-4s  %-10s  %-10s  %-7s  %-10s  %-10s  %-10s  %s\n",
    //        "Idx","H [km]","Va [km/s]","m/n","i [deg]","phi [deg]","om [deg]","Case");
    // for (int j=0; j<J; j++)
    //     printf("%-4d  %-10.1f  %-10.5f  %-7.3f  %-10.4f  %-10.4f  %-10.1f  %s\n",
    //            j+1,
    //            DATA[j].H,
    //            DATA[j].Va,
    //            DATA[j].mn,
    //            DATA[j].i_rad * 180.0/PI,
    //            DATA[j].phi   * 180.0/PI,
    //            DATA[j].omega,
    //            DATA[j].has_phi ? "general" : "special");
    // printf("──────────────────────────────────────────────────────────────────────────────────\n");
    // printf("Total DATA entries: %d\n\n", J);

    free(H); free(r); free(vc);
    free(va_bg); free(vi_bg);
    free(plot_cmd);

    *n_data_out = J;
    return DATA;
}

/* ── format_data_table ──────────────────────────────────────────────────── *
 *  Formats DATA array as a plain-text table for display in the GUI.         *
 *  Returns newly allocated string; caller must g_free().                   */
static gchar *format_data_table(const LGOData *DATA, int n, int use_om)
{
    GString *s = g_string_new(NULL);

    /* Header */
    // if (use_om)
        g_string_append_printf(s,
            "%-4s %-7s %-9s %-5s %-7s %-10s %-8s %-9s %-7s %s\n",
            "Idx","H[km]","Va[km/s]","m/n","i[deg]","phi[deg]","om[deg]","a[km]","e","Case");
    // else
        // g_string_append_printf(s,
            // "%-4s  %-10s  %-10s  %-7s  %-10s  %-10s  %s\n",
            // "Idx","H [km]","Va [km/s]","m/n [-]","i [deg]","phi [deg]","Case");

    g_string_append_printf(s,
        "───────────────────────────────────────────────────────────────────────────────────\n");

    for (int j = 0; j < n; j++) {
        // if (use_om)
            g_string_append_printf(s,
                "%-4d %-7.0f %-9.5f %-5.1f %-7.1f %-10.4f %-8.1f %-9.1f %-7.4f %s\n",
                j+1, 
                DATA[j].H, 
                DATA[j].Va, 
                DATA[j].mn,
                DATA[j].i_rad*180.0/PI, 
                DATA[j].phi*180.0/PI,
                DATA[j].omega,
                DATA[j].a,
                DATA[j].e,
                DATA[j].has_phi ? "general" : "special");
        // else
            // g_string_append_printf(s,
                // "%-4d  %-10.1f  %-10.5f  %-7.3f  %-10.4f  %-10.4f  %s\n",
                // j+1, DATA[j].H, DATA[j].Va, DATA[j].mn,
                // DATA[j].i_rad*180.0/PI, DATA[j].phi*180.0/PI,
                // DATA[j].has_phi ? "general" : "special");
    }

    g_string_append_printf(s,
        "───────────────────────────────────────────────────────────────────────────────────\n"
        "Total: %d orbit%s identified\n", n, n==1?"":"s");

    return g_string_free(s, FALSE);   /* transfers ownership of the char*   */
}

/* ══════════════════════════════════════════════════════════════════════════
 *  GUI
 * ══════════════════════════════════════════════════════════════════════════ */

typedef struct {
    /* Input widgets */
    GtkWidget *entry_I;          /* inclinations [deg]                       */
    GtkWidget *entry_mn0;        /* m/n targets [-]                          */
    GtkWidget *radio_special;    /* Z = 0 radio                              */
    GtkWidget *radio_general;    /* Z = 1 radio                              */
    GtkWidget *radio_om;         /* use omega as input                       */
    GtkWidget *radio_phi;        /* use phi as input                         */
    GtkWidget *entry_angles;     /* omega or phi values                      */
    GtkWidget *lbl_angle_unit;   /* "ω [deg]:" / "φ [deg]:"                 */
    GtkWidget *general_frame;    /* container greyed when special selected   */
    /* Output widgets */
    GtkWidget *textview;
    GtkWidget *lbl_status;
    GtkWidget *btn_run;
} AppWidgets;

/* ── parse_doubles ──────────────────────────────────────────────────────── *
 *  Splits a comma/space-separated string into a double array.               *
 *  Returns number of values parsed (≤ max_n).                               */
static int parse_doubles(const char *str, double *out, int max_n)
{
    int n = 0;
    char buf[512];
    g_strlcpy(buf, str, sizeof(buf));
    char *tok = strtok(buf, ", \t");
    while (tok && n < max_n) {
        out[n++] = atof(tok);
        tok = strtok(NULL, ", \t");
    }
    return n;
}

/* ── on_case_toggled ────────────────────────────────────────────────────── */
static void on_case_toggled(GtkToggleButton *btn, gpointer user_data)
{
    AppWidgets *w = user_data;
    gboolean general = gtk_toggle_button_get_active(
                           GTK_TOGGLE_BUTTON(w->radio_general));
    gtk_widget_set_sensitive(w->general_frame, general);
}

/* ── on_mode_toggled ────────────────────────────────────────────────────── */
static void on_mode_toggled(GtkToggleButton *btn, gpointer user_data)
{
    AppWidgets *w = user_data;
    gboolean om_selected = gtk_toggle_button_get_active(
                               GTK_TOGGLE_BUTTON(w->radio_om));
    gtk_label_set_text(GTK_LABEL(w->lbl_angle_unit),
                       om_selected ? "ω values [deg]:" : "φ values [deg]:");
    /* Update placeholder hint */
    gtk_entry_set_placeholder_text(GTK_ENTRY(w->entry_angles),
        om_selected ? "e.g.  135, 90" : "e.g.  53, 30");
}

/* ── on_run_clicked ─────────────────────────────────────────────────────── */
static void on_run_clicked(GtkButton *btn, gpointer user_data)
{
    AppWidgets *w = user_data;

    /* Disable button during run */
    gtk_widget_set_sensitive(w->btn_run, FALSE);
    gtk_label_set_text(GTK_LABEL(w->lbl_status), "Running…");
    while (gtk_events_pending()) gtk_main_iteration();

    /* ── Parse inputs ─────────────────────────────────────────────────── */
    double I_arr[16], mn0_arr[16], ang_arr[16];
    int n_I   = parse_doubles(gtk_entry_get_text(GTK_ENTRY(w->entry_I)),   I_arr,   16);
    int n_mn0 = parse_doubles(gtk_entry_get_text(GTK_ENTRY(w->entry_mn0)), mn0_arr, 16);
    int n_ang = parse_doubles(gtk_entry_get_text(GTK_ENTRY(w->entry_angles)), ang_arr, 16);

    int Z      = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w->radio_general)) ? 1 : 0;
    int use_om = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w->radio_om))      ? 1 : 0;

    /* ── Validate ─────────────────────────────────────────────────────── */
    if (n_I < 1 || n_mn0 < 1) {
        gtk_label_set_text(GTK_LABEL(w->lbl_status),
                           "⚠ Provide at least one inclination and one m/n target.");
        gtk_widget_set_sensitive(w->btn_run, TRUE);
        return;
    }
    if (Z == 1 && n_ang < 1) {
        gtk_label_set_text(GTK_LABEL(w->lbl_status),
                           "⚠ General case requires at least one angle value.");
        gtk_widget_set_sensitive(w->btn_run, TRUE);
        return;
    }

    /* ── Create data dir ──────────────────────────────────────────────── */
    MKDIR("data");

    /* ── Run computation ──────────────────────────────────────────────── */
    int n_data = 0;
    LGOData *DATA = LGO_C2(I_arr, n_I, mn0_arr, n_mn0,
                            ang_arr, n_ang, Z, use_om, &n_data);

    /* ── Display results ──────────────────────────────────────────────── */
    gchar *table = format_data_table(DATA, n_data, use_om);
    GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(w->textview));
    gtk_text_buffer_set_text(buf, table, -1);
    g_free(table);

    /* Also write to file */
    FILE *ft = fopen("data/table_results.txt","w");
    if (ft) {
        gchar *t2 = format_data_table(DATA, n_data, use_om);
        fputs(t2, ft);
        g_free(t2);
        fclose(ft);
    }

    free(DATA);

    /* ── Status ───────────────────────────────────────────────────────── */
    gchar status[128];
    g_snprintf(status, sizeof(status),
               "✓  %d orbit%s identified — gnuplot window launched",
               n_data, n_data==1?"":"s");
    gtk_label_set_text(GTK_LABEL(w->lbl_status), status);
    gtk_widget_set_sensitive(w->btn_run, TRUE);
}

/* ── activate — build the UI ────────────────────────────────────────────── */
static void activate(GtkApplication *app, gpointer user_data)
{
    AppWidgets *w = user_data;

    /* ── CSS styling ──────────────────────────────────────────────────── */
    GtkCssProvider *css = gtk_css_provider_new();
    gtk_css_provider_load_from_data(css,
        "button.run {"
        "  background: #0072BD; color: white;"
        "  font-weight: bold; border-radius: 4px;"
        "  padding: 5px 24px; border: none; }"
        "button.run:hover { background: #005a9e; }"
        // "textview text { font-family: monospace; font-size: 9pt; }"
        "frame > label { font-weight: bold; }",
        -1, NULL);
    gtk_style_context_add_provider_for_screen(
        gdk_screen_get_default(),
        GTK_STYLE_PROVIDER(css),
        GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);

    /* ── Main window ──────────────────────────────────────────────────── */
    GtkWidget *win = gtk_application_window_new(app);
    gtk_window_set_title(GTK_WINDOW(win), "Geosynchronous LGO Analysis");
    gtk_window_set_default_size(GTK_WINDOW(win), 660, 580);
    gtk_container_set_border_width(GTK_CONTAINER(win), 10);

    /* ── Root vertical box ────────────────────────────────────────────── */
    GtkWidget *root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
    gtk_container_add(GTK_CONTAINER(win), root);

    /* ══ PARAMETERS FRAME ════════════════════════════════════════════════ */
    GtkWidget *param_frame = gtk_frame_new("Parameters");
    gtk_box_pack_start(GTK_BOX(root), param_frame, FALSE, FALSE, 0);

    GtkWidget *grid = gtk_grid_new();
    gtk_grid_set_row_spacing(GTK_GRID(grid), 8);
    gtk_grid_set_column_spacing(GTK_GRID(grid), 10);
    gtk_widget_set_margin_start(grid, 10);
    gtk_widget_set_margin_end(grid, 10);
    gtk_widget_set_margin_top(grid, 8);
    gtk_widget_set_margin_bottom(grid, 10);
    gtk_container_add(GTK_CONTAINER(param_frame), grid);

    /* Row 0: Inclinations */
    GtkWidget *lbl_I = gtk_label_new("Inclinations I [deg]:");
    gtk_widget_set_halign(lbl_I, GTK_ALIGN_END);
    gtk_grid_attach(GTK_GRID(grid), lbl_I, 0, 0, 1, 1);

    w->entry_I = gtk_entry_new();
    gtk_entry_set_text(GTK_ENTRY(w->entry_I), "0, 63.4");
    gtk_entry_set_width_chars(GTK_ENTRY(w->entry_I), 30);
    gtk_widget_set_tooltip_text(w->entry_I, "Comma-separated inclination values in degrees");
    gtk_widget_set_hexpand(w->entry_I, TRUE);
    gtk_grid_attach(GTK_GRID(grid), w->entry_I, 1, 0, 1, 1);

    /* Row 1: m/n targets */
    GtkWidget *lbl_mn = gtk_label_new("m/n targets [-]:");
    gtk_widget_set_halign(lbl_mn, GTK_ALIGN_END);
    gtk_grid_attach(GTK_GRID(grid), lbl_mn, 0, 1, 1, 1);

    w->entry_mn0 = gtk_entry_new();
    gtk_entry_set_text(GTK_ENTRY(w->entry_mn0), "1.0, 1.5, 2.0");
    gtk_entry_set_width_chars(GTK_ENTRY(w->entry_mn0), 30);
    gtk_widget_set_tooltip_text(w->entry_mn0, "Comma-separated target m/n repetition factors");
    gtk_widget_set_hexpand(w->entry_mn0, TRUE);
    gtk_grid_attach(GTK_GRID(grid), w->entry_mn0, 1, 1, 1, 1);

    /* Row 2: separator */
    gtk_grid_attach(GTK_GRID(grid), gtk_separator_new(GTK_ORIENTATION_HORIZONTAL), 0, 2, 2, 1);

    /* Row 3: Case selector */
    GtkWidget *lbl_case = gtk_label_new("Case:");
    gtk_widget_set_halign(lbl_case, GTK_ALIGN_END);
    gtk_grid_attach(GTK_GRID(grid), lbl_case, 0, 3, 1, 1);

    GtkWidget *case_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 16);
    w->radio_special = gtk_radio_button_new_with_label(NULL, "Special  (φ = i)");
    w->radio_general = gtk_radio_button_new_with_label_from_widget(
                           GTK_RADIO_BUTTON(w->radio_special), "General");
    gtk_box_pack_start(GTK_BOX(case_box), w->radio_special, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(case_box), w->radio_general, FALSE, FALSE, 0);
    gtk_grid_attach(GTK_GRID(grid), case_box, 1, 3, 1, 1);

    /* Row 4: General case sub-frame */
    w->general_frame = gtk_frame_new("General case input");
    gtk_widget_set_sensitive(w->general_frame, FALSE);   /* greyed until selected */
    gtk_grid_attach(GTK_GRID(grid), w->general_frame, 0, 4, 2, 1);

    GtkWidget *gen_grid = gtk_grid_new();
    gtk_grid_set_row_spacing(GTK_GRID(gen_grid), 6);
    gtk_grid_set_column_spacing(GTK_GRID(gen_grid), 10);
    gtk_widget_set_margin_start(gen_grid, 10);
    gtk_widget_set_margin_end(gen_grid, 10);
    gtk_widget_set_margin_top(gen_grid, 6);
    gtk_widget_set_margin_bottom(gen_grid, 8);
    gtk_container_add(GTK_CONTAINER(w->general_frame), gen_grid);

    /* Input mode: ω or φ */
    GtkWidget *lbl_mode = gtk_label_new("Input via:");
    gtk_widget_set_halign(lbl_mode, GTK_ALIGN_END);
    gtk_grid_attach(GTK_GRID(gen_grid), lbl_mode, 0, 0, 1, 1);

    GtkWidget *mode_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 16);
    w->radio_om  = gtk_radio_button_new_with_label(NULL, "ω  arg. of perigee [deg]");
    w->radio_phi = gtk_radio_button_new_with_label_from_widget(
                       GTK_RADIO_BUTTON(w->radio_om), "φ  sub-satellite latitude [deg]");
    gtk_box_pack_start(GTK_BOX(mode_box), w->radio_om,  FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(mode_box), w->radio_phi, FALSE, FALSE, 0);
    gtk_grid_attach(GTK_GRID(gen_grid), mode_box, 1, 0, 1, 1);

    /* Angle values entry */
    w->lbl_angle_unit = gtk_label_new("ω values [deg]:");
    gtk_widget_set_halign(w->lbl_angle_unit, GTK_ALIGN_END);
    gtk_grid_attach(GTK_GRID(gen_grid), w->lbl_angle_unit, 0, 1, 1, 1);

    w->entry_angles = gtk_entry_new();
    gtk_entry_set_text(GTK_ENTRY(w->entry_angles), "135, 90");
    gtk_entry_set_placeholder_text(GTK_ENTRY(w->entry_angles), "e.g.  135, 90");
    gtk_widget_set_tooltip_text(w->entry_angles,
        "Comma-separated angle values (ω or φ depending on selection above)");
    gtk_widget_set_hexpand(w->entry_angles, TRUE);
    gtk_grid_attach(GTK_GRID(gen_grid), w->entry_angles, 1, 1, 1, 1);

    /* ══ RUN BUTTON + STATUS ═════════════════════════════════════════════ */
    GtkWidget *run_row = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 10);
    gtk_widget_set_margin_top(run_row, 2);
    gtk_box_pack_start(GTK_BOX(root), run_row, FALSE, FALSE, 0);

    w->btn_run = gtk_button_new_with_label("▶  Run");
    gtk_style_context_add_class(gtk_widget_get_style_context(w->btn_run), "run");
    gtk_box_pack_start(GTK_BOX(run_row), w->btn_run, FALSE, FALSE, 0);

    w->lbl_status = gtk_label_new("Ready");
    gtk_label_set_xalign(GTK_LABEL(w->lbl_status), 0.0);
    gtk_widget_set_hexpand(w->lbl_status, TRUE);
    gtk_box_pack_start(GTK_BOX(run_row), w->lbl_status, TRUE, TRUE, 0);

    /* ══ RESULTS FRAME ═══════════════════════════════════════════════════ */
    GtkWidget *res_frame = gtk_frame_new("Results");
    gtk_box_pack_start(GTK_BOX(root), res_frame, TRUE, TRUE, 0);

    GtkWidget *scroll = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroll),
                                   GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    gtk_widget_set_margin_start(scroll, 4);
    gtk_widget_set_margin_end(scroll, 4);
    gtk_widget_set_margin_top(scroll, 4);
    gtk_widget_set_margin_bottom(scroll, 4);
    gtk_container_add(GTK_CONTAINER(res_frame), scroll);

    w->textview = gtk_text_view_new();
    gtk_text_view_set_editable(GTK_TEXT_VIEW(w->textview), FALSE);
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(w->textview), FALSE);
    gtk_text_view_set_left_margin(GTK_TEXT_VIEW(w->textview), 6);
    gtk_text_view_set_top_margin(GTK_TEXT_VIEW(w->textview), 4);
    
    /* ── Force monospace so column padding aligns correctly ── */
    PangoFontDescription *mono = pango_font_description_from_string("Monospace 9");
    gtk_widget_override_font(w->textview, mono);
    pango_font_description_free(mono);

    gtk_container_add(GTK_CONTAINER(scroll), w->textview);

    /* Initial hint text */
    gtk_text_buffer_set_text(
        gtk_text_view_get_buffer(GTK_TEXT_VIEW(w->textview)),
        "Set parameters above and press ▶ Run.\n"
        "Results will appear here; the gnuplot window will open separately.", -1);

    /* ── Signal connections ───────────────────────────────────────────── */
    g_signal_connect(w->radio_general, "toggled", G_CALLBACK(on_case_toggled), w);
    g_signal_connect(w->radio_special, "toggled", G_CALLBACK(on_case_toggled), w);
    g_signal_connect(w->radio_om,      "toggled", G_CALLBACK(on_mode_toggled), w);
    g_signal_connect(w->radio_phi,     "toggled", G_CALLBACK(on_mode_toggled), w);
    g_signal_connect(w->btn_run,       "clicked", G_CALLBACK(on_run_clicked),  w);

    gtk_widget_show_all(win);
}

/* ══════════════════════════════════════════════════════════════════════════
 *  main
 * ══════════════════════════════════════════════════════════════════════════ */
int main(int argc, char **argv)
{
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    MKDIR("data");

    AppWidgets w = {0};

    GtkApplication *app = gtk_application_new("org.lgo.analysis",
                                               G_APPLICATION_DEFAULT_FLAGS);
    g_signal_connect(app, "activate", G_CALLBACK(activate), &w);

    int status = g_application_run(G_APPLICATION(app), argc, argv);
    g_object_unref(app);
    return status;
}
