// Originally written from Ben Jeurissen (ben.jeurissen@uantwerpen.be)
// under the supervision of Alexander Leemans (a.leemans@umcutrecht.nl)

#include <mex.h>
#include <string>
#include <cmath>
#include <limits>
// TEST
//#include <omp.h>

#define MAX_DIR_CHANGE 0.2
#define ANGLE_TOLERANCE 1e-4
#define M_PI 3.14159265358979323846264338327950288

static mxArray *precomp_legendre_persistent = NULL;
static float *precomp_legendre = NULL;
static int lmax_legendre = 0;
static int num_legendre_coefs = 0;
static int num_legendre_dirs = 0;
static float inc_legendre = 0.0;

inline int NforL (int lmax) { return ((lmax+1)*(lmax+2)/2); }
inline int index (int l, int m) { return (l*(l+1)/2 + m); }
inline int index_mpos (int l, int m) { return (l*l/4 + m); }
inline int LforN (int N) { return (2*(((int) (sqrt((float) (1+8*N)))-3)/4)); }
inline int L_mposforN(int N) { return (int)(2*sqrt((float)N)-2); }

class PrecomputedFraction {
public:
    PrecomputedFraction () : f1 (0.0), f2 (0.0), p1 (NULL), p2 (NULL) { }
    float f1, f2;
    float *p1, *p2;
};

inline void calc_index_fractions(PrecomputedFraction& f, float elevation) {
	f.f2 = elevation / inc_legendre;
    int index = (int) f.f2;
    if (index < 0) { index = 0; f.f1 = 1.0; f.f2 = 0.0; }
    else if (index >= num_legendre_dirs-1) { index = num_legendre_dirs-1; f.f1 = 1.0; f.f2 = 0.0; }
    else { f.f2 -= index; f.f1 = 1.0 - f.f2; }
    
    f.p1 = precomp_legendre + index*num_legendre_coefs;
    f.p2 = f.p1 + num_legendre_coefs;
}

inline float legendre_precomputed(const PrecomputedFraction& f, int l, int m) {
    int i (index_mpos (l,m));
	float retval = f.f1*f.p1[i];
	if (f.f2) retval += f.f2*f.p2[i];
    return (retval);
}

inline float eval(const float *SH, const int lmax, const float* dir) {
    float azimuth = atan2 (dir[1], dir[0]);
    float elevation = acos (dir[2]);
    
    float sel = sin(elevation);
    bool atpole = sel < 1e-4;
    float* legendre = new float[NforL(lmax)];
    
    float amplitude = 0.0;
    
    PrecomputedFraction f;
    calc_index_fractions (f, elevation);
    elevation = cos (elevation);
    
    for (int l = 0; l <= (int) lmax; l+=2) {
        for (int m = 0; m <= l; m++)
            legendre[index_mpos(l,m)] = legendre_precomputed (f, l, m);
        amplitude += SH[index(l,0)] * legendre[index_mpos(l,0)];
    }
    
    for (int m = 1; m <= lmax; m++) {
        float caz = cos (m*azimuth);
        float saz = sin (m*azimuth);
        for (int l = 2*((m+1)/2); l <= lmax; l+=2) {
            amplitude += SH[index(l,m)] * legendre[index_mpos(l,m)] * caz;
            amplitude += SH[index(l,-m)] * legendre[index_mpos(l,m)] * saz;
        }
    }
    delete[] legendre;
    return amplitude;
}

inline void rand_dir(const float* dir_in, const float dist_spread, float* dir_out) {
    float v[3];
    do { 
        v[0] = 2.0*((float)rand()/RAND_MAX) - 1.0; 
        v[1] = 2.0*((float)rand()/RAND_MAX) - 1.0; 
    } while (v[0]*v[0] + v[1]*v[1] > 1.0); 

	v[0] *= dist_spread;
	v[1] *= dist_spread;
	v[2] = 1.0 - (v[0]*v[0] + v[1]*v[1]);
	v[2] = v[2] < 0.0 ? 0.0 : sqrt (v[2]);

    if (dir_in[0]*dir_in[0] + dir_in[1]*dir_in[1] < 1e-4) {
        dir_out[0] = v[0];
        dir_out[1] = v[1];
        dir_out[2] = dir_in[2] > 0.0 ? v[2] : -v[2];
        return;
    }
    
    float y[] = { dir_in[0], dir_in[1], 0.0 };
    float a = sqrt((y[0]*y[0])+(y[1]*y[1])+(y[2]*y[2]));
    y[0] /= a; y[1] /= a; y[2] /= a;
	float x[] =  { -y[1], y[0], 0.0 };
	float y2[] = { -x[1]*dir_in[2], x[0]*dir_in[2], x[1]*dir_in[0] - x[0]*dir_in[1] };
	a = sqrt((y2[0]*y2[0])+(y2[1]*y2[1])+(y2[2]*y2[2]));
    y2[0] /= a; y2[1] /= a; y2[2] /= a;

    float cx = v[0]*x[0] + v[1]*x[1];
    float cy = v[0]*y[0] + v[1]*y[1];
    
    dir_out[0] = cx*x[0] + cy*y2[0] + v[2]*dir_in[0];
    dir_out[1] = cx*x[1] + cy*y2[1] + v[2]*dir_in[1];
    dir_out[2] = cy*y2[2] + v[2]*dir_in[2];
    return;
}

inline float sample_sh(const float* sh, const int lmax, const float* dir_in, float dist_spread, float min_val, float max_trials, float* dir_out) {
    float max_val = 0.0;
    float val = 0.0;
    for (int n = 0; n < 12; n++) {
        rand_dir(dir_in, dist_spread, dir_out);
        val = eval(sh, lmax, dir_out);
        if (val > max_val) max_val = val;
    }
    max_val *= 1.5;
    
    for (int n = 0; n < max_trials; n++) {
        rand_dir(dir_in, dist_spread, dir_out);
        val = eval(sh, lmax, dir_out);
        if (val > min_val){
            if (((float)rand()/RAND_MAX)*max_val < val) {
                return val;
            }
        }
    }
    dir_out[0] = std::numeric_limits<float>::quiet_NaN();
    dir_out[1] = std::numeric_limits<float>::quiet_NaN();
    dir_out[2] = std::numeric_limits<float>::quiet_NaN();
    return (std::numeric_limits<float>::quiet_NaN());
}

inline void derivatives (const float *SH, int lmax, float elevation, float azimuth, float &amplitude, float &dSH_del, float &dSH_daz, float &d2SH_del2, float &d2SH_deldaz, float &d2SH_daz2) {
    float sel = sin(elevation);
    bool atpole = sel < 1e-4;
    float* legendre = new float[NforL(lmax)];
    
    amplitude = dSH_del = dSH_daz = d2SH_del2 = d2SH_deldaz = d2SH_daz2 = 0.0;
    
    PrecomputedFraction f;
    calc_index_fractions (f, elevation);
    elevation = cos (elevation);
    
    for (int l = 0; l <= (int) lmax; l+=2) {
        for (int m = 0; m <= l; m++)
            legendre[index_mpos(l,m)] = legendre_precomputed (f, l, m);
        
        amplitude += SH[index(l,0)] * legendre[index_mpos(l,0)];
            
        if (l) {
            dSH_del += SH[index(l,0)] * sqrt((float) l*(l+1)/2.0) * legendre[index_mpos(l,1)];
            d2SH_del2 += SH[index(l,0)] * (
                    sqrt((float) l*(l+1)*(l-1)*(l+2)/2.0) * legendre[index_mpos(l,2)]
                    - l*(l+1) * legendre[index_mpos(l,0)] )/2.0;
        }
    }
    
    for (int m = 1; m <= lmax; m++) {
        float caz = cos (m*azimuth);
        float saz = sin (m*azimuth);
        for (int l = 2*((m+1)/2); l <= lmax; l+=2) {
            amplitude += SH[index(l,m)] * legendre[index_mpos(l,m)] * caz;
            amplitude += SH[index(l,-m)] * legendre[index_mpos(l,m)] * saz;
            
            float tmp = sqrt((float) (l+m)*(l-m+1))*legendre[index_mpos(l,m-1)];
            if (m == 1) tmp *= sqrt(2.0);
            if (l > m) tmp -= sqrt((float) (l-m)*(l+m+1))*legendre[index_mpos(l,m+1)];
            tmp /= -2.0;
            dSH_del += SH[index(l,m)] * tmp * caz;
            dSH_del += SH[index(l,-m)] * tmp * saz;
            
            float tmp2 = - ( (l+m)*(l-m+1) + (l-m)*(l+m+1) ) * legendre[index_mpos(l,m)];
            if (m == 1) {
                tmp2 -= (l+1)*l * legendre[index_mpos(l,1)];
            } else {
                if (m == 2) {
                    tmp2 += sqrt((float) (l+m)*(l-m+1)*(l+m-1)*(l-m+2)*2) * legendre[index_mpos(l,m-2)];
                } else {
                    tmp2 += sqrt((float) (l+m)*(l-m+1)*(l+m-1)*(l-m+2)) * legendre[index_mpos(l,m-2)];
                }
            }
            if (l > m+1) tmp2 += sqrt((float) (l-m)*(l+m+1)*(l-m-1)*(l+m+2)) * legendre[index_mpos(l,m+2)];
            tmp2 /= 4.0;
            d2SH_del2 += SH[index(l,m)] * tmp2 * caz;
            d2SH_del2 += SH[index(l,-m)] * tmp2 * saz;
            
            if (atpole) {
                dSH_daz -= SH[index(l,m)] * tmp * saz;
                dSH_daz += SH[index(l,-m)] * tmp * caz;
            } else {
                d2SH_deldaz -= m * SH[index(l,m)] * tmp * saz;
                d2SH_deldaz += m * SH[index(l,-m)] * tmp * caz;
                
                dSH_daz -= m * SH[index(l,m)] * legendre[index_mpos(l,m)] * saz;
                dSH_daz += m * SH[index(l,-m)] * legendre[index_mpos(l,m)] * caz;
                
                tmp =  m*m * legendre[index_mpos(l,m)];
                d2SH_daz2 -= SH[index(l,m)] * tmp * caz;
                d2SH_daz2 -= SH[index(l,-m)] * tmp * saz;
            }
        }
    }
    delete[] legendre;
    
    if (!atpole) {
        dSH_daz /= sel;
        d2SH_deldaz /= sel;
        d2SH_daz2 /= sel*sel;
    }
}

inline float get_peak (const float* SH, int lmax, float* unit_init_dir) {
    float amplitude, dSH_del, dSH_daz, d2SH_del2, d2SH_deldaz, d2SH_daz2;
    float az = atan2 (unit_init_dir[1], unit_init_dir[0]);
    float el = acos (unit_init_dir[2]);
    float del, daz, dSH_dt, d2SH_dt2, dt;
    
    for (int i = 0; i < 50; i++) {
        az = atan2 (unit_init_dir[1], unit_init_dir[0]);
        el = acos (unit_init_dir[2]);
        derivatives (SH, lmax, el, az, amplitude, dSH_del, dSH_daz, d2SH_del2, d2SH_deldaz, d2SH_daz2);
        
        del = sqrt (dSH_del*dSH_del + dSH_daz*dSH_daz);
        daz = dSH_daz/del;
        del = dSH_del/del;
        
        dSH_dt = daz*dSH_daz + del*dSH_del;
        d2SH_dt2 = daz*daz*d2SH_daz2 + 2.0*daz*del*d2SH_deldaz + del*del*d2SH_del2;
        dt = - dSH_dt / d2SH_dt2;
        
        if (dt < 0.0 || dt > MAX_DIR_CHANGE) dt = MAX_DIR_CHANGE;
        
        del *= dt;
        daz *= dt;
        
        unit_init_dir[0] += del*cos(az)*cos(el) - daz*sin(az);
        unit_init_dir[1] += del*sin(az)*cos(el) + daz*cos(az);
        unit_init_dir[2] += -del*sin(el);
        
        float a = sqrt((unit_init_dir[0]*unit_init_dir[0])+(unit_init_dir[1]*unit_init_dir[1])+(unit_init_dir[2]*unit_init_dir[2]));
        unit_init_dir[0] /= a; unit_init_dir[1] /= a; unit_init_dir[2] /= a;
        
        if (dt < ANGLE_TOLERANCE) return (amplitude);
    }
    
    unit_init_dir[0] = std::numeric_limits<float>::quiet_NaN();
    unit_init_dir[1] = std::numeric_limits<float>::quiet_NaN();
    unit_init_dir[2] = std::numeric_limits<float>::quiet_NaN();
    return (std::numeric_limits<float>::quiet_NaN());
}

inline void peaks(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 2) {
        mexErrMsgTxt("Need 2 input parameters.");
    }
    if (nlhs != 2) {
        mexErrMsgTxt("Need 2 output parameters.");
    }
    const mwSize* sh_size = mxGetDimensions(prhs[0]);
    const mwSize* dir_size = mxGetDimensions(prhs[1]);
    if (sh_size[1] != dir_size[1]) {
        mexErrMsgTxt("Dimensions not equal.");
    }
    if (!mxIsSingle(prhs[0]) || !mxIsSingle(prhs[1])) {
        mexErrMsgTxt("Parameters must be of type single.");
    }
    plhs[0] = mxDuplicateArray(prhs[1]);
    plhs[1] = mxCreateNumericMatrix(1, dir_size[1], mxSINGLE_CLASS, mxREAL);
    
    const mwSize* dir2_size = mxGetDimensions(plhs[0]);
    const mwSize* val_size = mxGetDimensions(plhs[1]);
    
    float* sh = (float*) mxGetData(prhs[0]);
    float* dir = (float*) mxGetData(plhs[0]);
    float* val = (float*) mxGetData(plhs[1]);
    
    int lmax = LforN(sh_size[0]);
    
    if (lmax != lmax_legendre) {
        mexErrMsgTxt("Correctly initialize this funtion first!");
    }
    
    //#pragma omp parallel for
    for (int pos = 0; pos < sh_size[1]; pos++) {
        val[pos] = get_peak(&sh[pos*sh_size[0]], lmax, &dir[pos*dir_size[0]]);
    }
}

inline void sample(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 5) {
        mexErrMsgTxt("Need 5 input parameters.");
    }
    
    if (nlhs != 2) {
        mexErrMsgTxt("Need 2 output parameters.");
    }
    
    const mwSize* sh_size = mxGetDimensions(prhs[0]);
    const mwSize* dir_size = mxGetDimensions(prhs[1]);    
    float* sh = (float*) mxGetData(prhs[0]);
    float* dir_in = (float*) mxGetData(prhs[1]);
    float* max_angle = (float*) mxGetData(prhs[2]);
    float* min_val = (float*) mxGetData(prhs[3]);
    float* max_trials = (float*) mxGetData(prhs[4]);
    
    if (sh_size[1] != dir_size[1]) {
        mexErrMsgTxt("Dimensions not equal.");
    }
    
    if (!mxIsSingle(prhs[0]) || !mxIsSingle(prhs[1]) || !mxIsSingle(prhs[2]) || !mxIsSingle(prhs[3]) || !mxIsSingle(prhs[4])) {
        mexErrMsgTxt("Parameters must be of type single.");
    }
    plhs[0] = mxCreateNumericMatrix(3, dir_size[1], mxSINGLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericMatrix(1, dir_size[1], mxSINGLE_CLASS, mxREAL);
    float* dir_out = (float*) mxGetData(plhs[0]);
    float* val = (float*) mxGetData(plhs[1]);
    
    float max_angle_ = (M_PI / 180.00) * max_angle[0];
    float sin_max_angle = sin(max_angle_);
    
    int lmax = LforN(sh_size[0]);
          
    for (int pos = 0; pos < sh_size[1]; pos++) {
        val[pos] = sample_sh(&sh[pos*sh_size[0]], lmax, &dir_in[pos*dir_size[0]], sin_max_angle, min_val[0], max_trials[0], &dir_out[pos*dir_size[0]]);
    }
}

inline void eval(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 2) {
        mexErrMsgTxt("Need 2 input parameters.");
    }
    
    if (nlhs != 1) {
        mexErrMsgTxt("Need 1 output parameters.");
    }
    
    const mwSize* F_sh_size = mxGetDimensions(prhs[0]);
    const mwSize* dir_size = mxGetDimensions(prhs[1]);
    
    if (F_sh_size[1] != dir_size[1]) {
        mexErrMsgTxt("Dimensions not equal.");
    }
    
    if (!mxIsSingle(prhs[0]) || !mxIsSingle(prhs[1])) {
        mexErrMsgTxt("Parameters must be of type single.");
    }
    plhs[0] = mxCreateNumericMatrix(1, dir_size[1], mxSINGLE_CLASS, mxREAL);
    
    float* F_sh = (float*) mxGetData(prhs[0]);
    float* dir = (float*) mxGetData(prhs[1]);
    float* val = (float*) mxGetData(plhs[0]);
    
    int lmax = LforN(F_sh_size[0]);
    
    for (int pos = 0; pos < F_sh_size[1]; pos++) {
        val[pos] = eval(&F_sh[pos*F_sh_size[0]], lmax, &dir[pos*dir_size[0]]);
    }
}

inline void cleanup(void) {
    mxDestroyArray(precomp_legendre_persistent);
}

inline void init(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 1) {
        mexErrMsgTxt("Need 1 input parameter.");
    }
    
    if (!mxIsSingle(prhs[0])) {
        mexErrMsgTxt("Parameter must be of type single.");
    }
    if (precomp_legendre) {
        mxDestroyArray(precomp_legendre_persistent);
    }
    precomp_legendre_persistent = mxDuplicateArray(prhs[0]);
    precomp_legendre = (float*) mxGetData(precomp_legendre_persistent);
    mexMakeArrayPersistent(precomp_legendre_persistent);
    mexAtExit(cleanup);
    lmax_legendre = L_mposforN(mxGetM(precomp_legendre_persistent));
    num_legendre_coefs = mxGetM(precomp_legendre_persistent);
    num_legendre_dirs = mxGetN(precomp_legendre_persistent);
    inc_legendre = M_PI/(num_legendre_dirs-1);
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray *prhs[]) {
    std::string mode;
    if (nrhs >= 1) {
        if (!mxIsChar(prhs[0])) {
            mexErrMsgTxt("First parameter must be a string!");
        }
        mode = std::string(mxArrayToString(prhs[0]));
    } else {
        mexErrMsgTxt("Please specify a valid mode of operation.\n");
    }
    if (mode == "init") {
        init(nlhs, plhs, nrhs-1, prhs+1);
    } else {
        if (!precomp_legendre) {
            mexErrMsgTxt("Initialize this funtion first!");
        }
        if (mode == "peaks") {
            peaks(nlhs, plhs, nrhs-1, prhs+1);
        } else if (mode == "sample") {
            sample(nlhs, plhs, nrhs-1, prhs+1);
        } else if (mode == "eval") {
            eval(nlhs, plhs, nrhs-1, prhs+1);
        } else {
            mexErrMsgTxt("Please specify a valid mode of operation.\n");
        }
    }
}