#ifndef MB_H
#define MB_H
#include <vector>

typedef int N;
typedef long double R_lo;
//typedef double R_lolo;
#define GD_NONE 0
#define GD_PDB 1
#define GD_KN 2
#define GD_KNPDB 3
#define GD_GERRIT 4
#define GD_GERRIT_B 5

namespace mpfr{
class mpreal;
}
typedef mpfr::mpreal R_hi;

template <typename T> class fcomplex;
typedef fcomplex<R_lo> C_lo;

class series_approximation;
class reference;
class perturbation;
class MBImage;
class MB{
private:
    series_approximation *mp_s;
    reference *mp_r;
    MBImage *mp_img;
    N m_skip;
    N m_precision;
    N m_m;
    N m_maxiters;
    R_lo m_sradius;
    mpfr::mpreal *mp_locRe, *mp_locIm;

public:
    bool m_solveGlitches;
    R_lo m_allowedErrorSA;
    R_lo m_allowedErrorGC;

    typedef enum {
        S_idle = 0,
        S_busy,
        S_paused,
        S_aborted,
        S_finished
    } status;
    typedef enum {
        R_nothing = 0,
        R_doWork,
        R_pause,
        R_abort,
        R_shutdown,
    } requests;

    status m_status;
    requests m_request;
    bool m_workRequested;
    N m_seriesNbr;
    N m_referenceNbr;
    N m_imageNbr;
    N m_glitchNbr;
    N m_glitchPasses;
    N m_maxGlitchPasses;
    unsigned int m_time;
    N m_minIterNbr;
    N m_maxIterNbr;
    N m_GDmethod;

public:
    MB(){
        mp_s=nullptr; mp_r=nullptr; mp_img=nullptr; mp_locRe=nullptr; mp_locIm=nullptr;
        m_maxiters = 60000; m_m = 32; m_precision = 1000; m_skip = -1;
        m_status = S_idle; m_request = R_nothing; m_workRequested = false;
        m_seriesNbr = 0; m_referenceNbr = 0; m_imageNbr = 0; m_glitchNbr = 0;
        m_time = 0; m_minIterNbr = 0; m_maxIterNbr = 0;
        m_solveGlitches = false;
        m_allowedErrorSA = 0.1;
        m_allowedErrorGC = 0.0001;
        m_GDmethod = GD_NONE;
        m_glitchPasses = 0;
        m_maxGlitchPasses = 64;
    }
    ~MB();
    //void ComputeMB(const char *sre, const char *sim, R_lo sradius,
    //              MBImage &img, N maxiters=60000, N m=64 /*series approx order*/, N precision=1000, N skip=-1);
    void setLocation(const char *sre, const char *sim);
    void setLocRel(R_lo dx, R_lo dy, R_lo zom);//(dx,dy) is assumed to be in pixel coordinates
    void win2tex(R_lo dx, R_lo dy, R_lo &tx, R_lo &ty);// (dx,dy) in C-plane relative to center of screen. (tx,ty) in image (pixel) coordinates
    void setSRadius(R_lo sradius){m_sradius = sradius;}
    R_lo getSRadius(){return m_sradius;}
    void setSRadius(const char *sradius);
    void setImage(MBImage &img){mp_img = &img;}
    void setSeriesApproxOrder(N m){m_m = m;}
    void setMaxIter(N maxiters){m_maxiters = maxiters;}
    N    getMaxIter(){return m_maxiters;}
    void setSkipIter(N skip){m_skip = skip;}
    bool isOK();
    int worker();
    void dowork(){
        m_status=S_busy;
        ComputeMB2(*mp_locRe, *mp_locIm, m_sradius, *mp_img, m_maxiters, m_m, m_precision, m_skip);
        m_status=S_finished;
    }
    void setPrecision(N precision);
    bool load(const char *fileName);
    bool save(const char *filename);
    unsigned int getNbrRoots();
    bool getRoot(unsigned int n, R_lo &x, R_lo &y, N &period, N &quality);
    bool getRootF(unsigned int n, float &px, float &py);
    void drawRoots();
    C_lo Scr2C(int x, int y);
    void setRequest(requests req){m_request = req;}
    status getStatus(){ return m_status;}

private:
    void newMP_s(R_lo allowedErrorSA, R_lo tmax);
    void newMP_r();
    void ComputeMB2(mpfr::mpreal &sre, mpfr::mpreal &sim, R_lo sradius,
                  MBImage &img, N maxiters=60000, N m=64 /*series approx order*/, N precision=1000, N skip=-1);
    solve_glitches(R_lo dt, std::vector<perturbation> &glitches);
};

#endif // MB_H
