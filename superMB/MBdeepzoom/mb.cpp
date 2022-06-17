// m.cpp 2016 Claude Heiland-Allen
// LICENSE: public domain, cc0, whatever
// this is demonstration code based on algorithms found in fractalforums.com
// glitch detection method: Pauldelbrot + knighty + Gerrit
// series stopping condition: knighty (modified)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Some modifications by knighty (just because I use msys32 under windows and don't have boost library):
//  - uses mpreal.h (depends on mpfr which depends on gmp: That is why I added gmp to the linked libraries) 
//  - Use scaling in series approximation. Now can use m=64 or more and zooms deeper :) . Issues with denormals and underflow solved. Now much faster. Update: Now much simpler formula :))
//  - Added different methods for glitch detection.
//  - Glitch correction changed to a very stupid method that works better than smart ones. ;o)
//  - A probe based SA stop criterion. It uses the periodic points that are found while computing the SA as new probes (actually nearby points are used as probes). Works well!
//  - Use custom lightweight complex<> template class to gain more speed? The standard complex<> template class seem to prevent some optimizations (unless using double floats)? todo: verify
//  - Added a GUI using ImGUI.
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// TODO:
//  - Fixme in perturbation class : check escape: What happens if the current pixel needs more iterations than the reference
//  - More refined tests. By divide and conquere? By solving for roots of the derivative of the SA?
//  - solve_glitches() : Do a non vector<> based one?
//  - Fix the case when the reference diverges prematurally in series approx and in perturbation ?
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//Windows specific?--------------------------------------------
#undef RC_INVOKED
#undef __STRICT_ANSI__
#include <float.h>
//--------------------------------------------

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <cmath>
#include <limits>

#include <iostream>
//#include <fstream>
//#include <array>

#if 0
#include <complex>
#define COMPLEX complex
#else
//A little problem with mpreal which includes <complex> so I had to change the name of the custom complex class to fcomplex
#include "complex.h"
#define COMPLEX fcomplex
#endif

#include <vector>

#include "mpfr/mpreal.h"

#include "mbimage.h"

#include "mb.h"

using namespace std;
using namespace mpfr;

// type aliases
//typedef mpreal R_hi;
typedef COMPLEX<R_lo> C_lo;
//typedef COMPLEX<R_lolo> C_lolo;
typedef COMPLEX<R_hi> C_hi;

//escape radius squared
#define ESCAPE_R2 256
#define DE_COLORING
#define SOLVE_GLITCHES
//#define RESCALE_ONCE

// Conversion routines
inline R_lo r_lo(const R_hi &z)
{
  return z.toLDouble();//MPFR_RNDN);//mpfr_get_ld(z.backend().data(), MPFR_RNDN);
}

inline R_lo r_lo(const char *s)
{
  return strtold(s, NULL);
}

inline C_lo c_lo(const C_hi &z)
{
  return C_lo(r_lo(real(z)), r_lo(imag(z)));
}

//-------------------------------------------------
// FPU control (experimental)
// ToDo: It seems those functions are spesific to Windows. How to do it on Unix systems?
//#ifdef windows
unsigned int GetFPUCW(){
    return _control87(0, 0);
}
void FPUPlay(){
    _control87(_PC_64, _MCW_PC);//Set precision to 64bits
    //_control87(0, _MCW_EM);//no exceptions... causes the program to halt prematurely
    //_control87(~_EM_DENORMAL & ~_EM_UNDERFLOW & ~_EM_INEXACT, _MCW_EM);
}
void FPURoundUp(){
    _control87(_RC_UP, _MCW_RC);//set denormal to save
}
void FPURoundNear(){
    _control87(_RC_NEAR, _MCW_RC);//set denormal to save
}
void RestoreFPUCW( unsigned int cw){
    _control87(cw, 0xFFFFFFFF);
}
//#endif
//-------------------------------------------------

//--------------------------------------------------------------------
// reference orbit for perturbation rendering
//--------------------------------------------------------------------
class reference
{
public:
  N n0;
  N n;
  C_hi z;
  C_hi c;
  vector<C_lo> z_ref;
  vector<R_lo> z_size;  // scaled absolute value for glitch detection

  // constructor
  reference(N n, C_hi z, C_hi c)
  : n0(n)
  , n(n)
  , z(z)
  , c(c)
  {
    C_lo z_lo(c_lo(z));
    z_ref.push_back(z_lo);
    z_size.push_back(1.0e-6 * norm(z_lo));
  }

  // reset
  // ToDo: remove?
  void reset(N _n, C_hi _z, C_hi _c)
  {
    n0 = _n; n = _n;
    z  = _z; c = _c;
    z_ref.clear(); z_size.clear();
    C_lo z_lo(c_lo(z));
    z_ref.push_back(z_lo);
    z_size.push_back(1.0e-6 * norm(z_lo));
  }

  // step a reference one iteration
  // returns false if reference escaped
  bool step()
  {
    n = n + 1;
    z = z * z + c;
    C_lo z_lo(c_lo(z));
    z_ref.push_back(z_lo);
    z_size.push_back(1.0e-6 * norm(z_lo));
    return (norm(z_lo) <= R_lo(ESCAPE_R2));
  }

  // step a reference until it escapes or maxiters reached
  // return false if escaped
  N run(N maxiters, N &n1, const MB::requests & request)
  {
    n1=0;
    while (n < maxiters && step() && request != MB::R_abort){
      n1++;
	}
    return n1;// n== maxiters;
  }
};
//--------------------------------------------------------------------
// perturbation technique per-pixel iterations with
//--------------------------------------------------------------------
class perturbation
{
    C_lo z;
    R_lo nz;
    C_lo dzdc1;
    C_lo dz1;
public:
  N x;
  N y;
  N n;
  C_lo dc;
  C_lo dz;
  C_lo dzdc;
  R_lo dt;//allowed error in Knighty's glitch detection
  R_lo err;
  bool glitched;
  bool escaped;
  N m_sn;
  R_lo m_sa;
  C_lo m_sz;
  C_lo m_sdz;

  // constructor initializes
  perturbation(N _x, N _y, N _n, C_lo _dc, C_lo _dz, C_lo _dzdc, R_lo _dt)
  : x(_x)
  , y(_y)
  , n(_n)
  , dc(_dc)
  , dz(_dz)
  , dzdc(_dzdc)
  , dt(_dt)
  , glitched(false)
  , escaped(false)
  , m_sn(0)
  , m_sa(1e100)
  { }

#define TOW_POW_PMAX64SQR 2.9387358770557187699218413430556e-39
#define TOW_POW_PMAX64 5.4210108624275221700372640043497e-20

  bool preblc(const reference &r, N maxiters){
      if(n - r.n0 >= N(r.z_ref.size()) ) {
          err = R_lo(maxiters - n);
          glitched = true;
          return true;
      }
      z = dz + r.z_ref[n - r.n0];
      nz = norm(z);
      dzdc1 = R_lo(2) * z * dzdc + C_lo(1);
      dz1 = (R_lo(2) * r.z_ref[n - r.n0] + dz) * dz + dc;
      return false;
  }

  bool haveEscaped(){
      if (nz > R_lo(ESCAPE_R2))
      {
        escaped = true;
        dz = z;
        return true;
      }
      return false;
  }

  void postblc(){
      if(m_sa > nz){
          m_sn = n;
          m_sa = nz;
          m_sz = dz;
          m_sdz = dzdc;
      }
      dzdc = dzdc1;
      dz   = dz1;
  }

  //Different glitch detection methods
  void runNoGD(const reference &r, N maxiters)
  {
    while (n < maxiters)
    {
        if (preblc(r,maxiters))
            break;

        if(haveEscaped())
            break;

        postblc();
        n = n + 1;
    }
  }

  void runPdbGD(const reference &r, N maxiters)
  {
    while (n < maxiters)
    {
        if (preblc(r,maxiters))
            break;

// Pauldelbrot's glitch detection heuristic
      if (nz < r.z_size[n - r.n0]){
        glitched = true;
        err = r.z_size[n - r.n0] / nz;
        dz = z;
        break;
      }

      if(haveEscaped())
          break;

      postblc();
      n = n + 1;
    }
  }

  void runKnGD(const reference &r, N maxiters)
  {
    while (n < maxiters)
    {
        if (preblc(r,maxiters))
            break;

      {
// Knighty's glitch detection heuristic... more accurate but much more expensive
// There is some room for optimizations
#if 0
          R_lo a2zn = R_lo(4.)*norm(r.z_ref[n - r.n0]);//norm(r.z_ref[n - r.n0]) is r.z_size[n - r.n0]
          R_lo adzn = norm(dz);
          R_lo adzn1= norm(dz1);
          R_lo a2znPdzn_adzdcn = norm((R_lo(2) * r.z_ref[n - r.n0] + dz) * dzdc1);
          if (TOW_POW_PMAX64SQR * (a2zn>adzn ? a2zn:adzn) * adzn1 > a2znPdzn_adzdcn * dt * dt)
          {
              glitched = true;
              err = r.z_size[n - r.n0] / nz;//TOW_POW_PMAX64SQR * (a2zn>adzn ? a2zn:adzn) * adzn1 / (a2znPdzn_adzdcn * dt * dt);
              dz = z;
              break;
          }
#else
          R_lo a2zn = R_lo(2.)*norm0(r.z_ref[n - r.n0]);//norm(r.z_ref[n - r.n0]) is r.z_size[n - r.n0]
          R_lo adzn = norm0(dz);
          R_lo adzn1= norm0(dz1);
          R_lo a2znPdzn_adzdcn = norm0((R_lo(2) * r.z_ref[n - r.n0] + dz) * dzdc1);
          if (TOW_POW_PMAX64 * (a2zn>adzn ? a2zn:adzn) * adzn1 > a2znPdzn_adzdcn * dt)
          {
              glitched = true;
              err = r.z_size[n - r.n0] / nz;
              dz = z;
              break;
          }
#endif
      }

      if(haveEscaped())
          break;

      postblc();
      n = n + 1;
    }
  }

  void runGerritGD(const reference &r, N maxiters)
  {
    while (n < maxiters)
    {
        if (preblc(r,maxiters))
            break;

      {
// Gerrit's glitch detection heuristic.
        R_lo azn = norm0(r.z_ref[n - r.n0]);
        R_lo adzn = norm0(dz);
        R_lo adzdc1= norm0(dzdc1);
        if (TOW_POW_PMAX64 * (adzn + R_lo(2) * azn) * adzn > adzdc1 * dt)
        {
            glitched = true;
            err = r.z_size[n - r.n0] / nz;//TOW_POW_PMAX64SQR * (a2zn>adzn ? a2zn:adzn) * adzn1 / (a2znPdzn_adzdcn * dt * dt);
            dz = z;
            break;
        }
      }

        if(haveEscaped())
            break;

        postblc();
      n = n + 1;
    }
  }

  void runGerritGDb(const reference &r, N maxiters)
  {
    R_lo S(0);
    while (n < maxiters)
    {
        if (preblc(r,maxiters))
            break;

      {
// Gerrit's glitch detection heuristic bu summing errors at each iteration:
        R_lo azn = norm0(r.z_ref[n - r.n0]);
        R_lo adzn = norm0(dz);
        R_lo numerator = adzn * (norm0(dz + R_lo(2) * r.z_ref[n - r.n0]) + adzn + R_lo(2)* azn) + norm0(dc);
        R_lo denominator = norm0(dzdc1);//norm0(R_lo(2) * dzdc * (dz + r.z_ref[n - r.n0]) + R_lo(1));
        S += numerator/denominator;
        if (S * TOW_POW_PMAX64 > dt)
        {
            glitched = true;
            err = r.z_size[n - r.n0] / nz;
            dz = z;
            break;
        }
      }

        if(haveEscaped())
            break;

        postblc();
      n = n + 1;
    }
  }

  void runKnPdbGD(const reference &r, N maxiters)
  {
    while (n < maxiters)
    {
        if (preblc(r,maxiters))
            break;

// Pauldelbrot's glitch detection heuristic as pre-condition
      if (nz < r.z_size[n - r.n0]){
// Knighty's glitch detection heuristic... more accurate but much more expensive
// Thereis some room for optimizations
        R_lo a2zn = R_lo(4.)*norm(r.z_ref[n - r.n0]);//norm(r.z_ref[n - r.n0]) is r.z_size[n - r.n0]
        R_lo adzn = norm(dz);
        R_lo adzn1= norm(dz1);
        R_lo a2znPdzn_adzdcn = norm((R_lo(2) * r.z_ref[n - r.n0] + dz) * dzdc1);
        if (TOW_POW_PMAX64SQR * (a2zn>adzn ? a2zn:adzn) * adzn1 > a2znPdzn_adzdcn * dt * dt)
        {
            glitched = true;
            err = r.z_size[n - r.n0] / nz;//TOW_POW_PMAX64SQR * (a2zn>adzn ? a2zn:adzn) * adzn1 / (a2znPdzn_adzdcn * dt * dt);
            dz = z;
            break;
        }
      }

      if(haveEscaped())
          break;

      postblc();
      n = n + 1;
    }
  }
// calculate pixel
// Fixme: what happens if the current pixel needs more iterations than the reference?
  void run(const reference &r, N maxiters)
  {
    while (n < maxiters)
    {
        if (preblc(r,maxiters))
            break;

// Pauldelbrot's glitch detection heuristic
#if 0
      if (nz < r.z_size[n - r.n0])
#endif
      {
#if 0
        glitched = true;
        err = r.z_size[n - r.n0] / nz;
        dz = z;
        break;
#else
// Knighty's glitch detection heuristic... more accurate but much more expensive
// There is some room for optimizations
        //#define TOW_POW_PMAX64SQR 2.9387358770557187699218413430556e-39
        R_lo a2zn = R_lo(4.)*norm(r.z_ref[n - r.n0]);//norm(r.z_ref[n - r.n0]) is r.z_size[n - r.n0]
        R_lo adzn = norm(dz);
        R_lo adzn1= norm(dz1);
        R_lo a2znPdzn_adzdcn = norm((R_lo(2) * r.z_ref[n - r.n0] + dz) * dzdc1);
        if (TOW_POW_PMAX64SQR * (a2zn>adzn ? a2zn:adzn) * adzn1 > a2znPdzn_adzdcn * dt * dt)
		{
			glitched = true;
            err = r.z_size[n - r.n0] / nz;//TOW_POW_PMAX64SQR * (a2zn>adzn ? a2zn:adzn) * adzn1 / (a2znPdzn_adzdcn * dt * dt);
			dz = z;
			break;
		}
#endif
      }

      if(haveEscaped())
          break;

      postblc();
      n = n + 1;
    }
	
  }
#undef TOW_POW_PMAX64SQR
#undef TOW_POW_PMAX64
};

//--------------------------------------------------------------------
// Class for probes used to determine when to stop SA computations
//--------------------------------------------------------------------
class probeClass
{
    C_lo m_dzz;
public:
    C_lo m_dc;
    C_lo m_dz;
    probeClass(){}
    probeClass(C_lo dc, C_lo dz):m_dc(dc),m_dz(dz){ m_dzz = C_lo(0);}
    probeClass(const probeClass & pC){ m_dc = pC.m_dc; m_dz = pC.m_dz; m_dzz = pC.m_dzz;}
    void advance(const C_lo &rz){
        m_dzz = m_dz * (R_lo(2) * rz + m_dz) + m_dc;
    }
    void shift(const C_lo &s, const C_lo &sz){
        m_dc -= s;
        m_dz -= sz;
        m_dzz = m_dz;//necessary to avoid the "floral fantasy" glitch. Also because areProbesOk uses dif2ex() below which compares against d_zz.
    }
    void apply(){m_dz = m_dzz;}
    C_lo dif2ex(const C_lo &ex){ return m_dzz - ex;}
};

//--------------------------------------------------------------------
// to store informations about the roots that are found during SA
//--------------------------------------------------------------------
class rootInfoClass
{
public:
    N n; //At which iteration this root was found
    N mult; //Multiplicity of the root... for later use
    C_lo RPos; //position of the root w.r.t. current reference point
    rootInfoClass(){}
    rootInfoClass(N _n, N _mult, const C_lo &_RPos):
        n(_n), mult(_mult), RPos(_RPos) {}
};

//--------------------------------------------------------------------
// series approximation with modification of
// obsolete: knighty's truncation error stopping condition
// Edit: scaled version. Edit02: the scaling is done every time R is near overflow. Now there are few if any denormals --> much faster than previous scaling.
// Edit: New SA stopping based on (modified) good old probes method.
//--------------------------------------------------------------------
class series_approximation
{
public:
  N n;             // iteration count
  R_lo dt;         // pixel size
  R_lo dy;         // tolerance on series approx relative error when evaluated
  R_lo tmax;
  R_lo Scl,iScl;   //Scale factors. iScl == 1/Scl
  R_lo RR;
  C_hi z;          // reference z at high precision
  C_hi c;          // reference c at high precision
  C_lo shift0;
  N m;             // series order
  vector<C_lo> a;  // coefficients (including z at start and R at end)
  vector<R_lo> abs_a; //For findRootsRel
  vector<R_lo> err_a;
  bool m_autoshift;
  vector<rootInfoClass> rts;
  vector<probeClass> m_probes;
  R_lo eps1;
  R_lo m_dt;//original (non scaled) dt

  // constructor initializes coefficients
  series_approximation(C_hi _c, R_lo _dt, R_lo _tmax, N _m, C_lo Shift0=C_lo(0), bool autoshift=true)
  : n(1)           // simplifies some things to start at iteration 1 with z = c
  , dt(_dt)
  , dy(1e-10)
  , tmax(_tmax)
  , Scl(1)
  , iScl(1)
  , RR(0)
  , z(_c)
  , c(_c)
  , shift0(Shift0)
  , m(_m)
  , a(_m + 2)
  , abs_a(_m + 2)
  , err_a(_m + 2)
  , m_autoshift(autoshift)
  , rts(0)
  , m_probes(0)
  {
    assert(_dt > 0);
    assert(_tmax > 0);
    assert(_m >= 1);
    c += shift0;
    z += shift0;
    a[0] = c_lo(c);
    a[1] = C_lo(1);
	eps1 = R_lo(1)+R_lo(16)*numeric_limits<R_lo>::epsilon();
    for (N k = 2; k <= m + 1; ++k)
      a[k] = C_lo(0); // includes truncation error R at end
    for (N k = 0; k <= m + 1; ++k)
      err_a[k] = R_lo(0);
    m_dt = dt;
#ifdef RESCALE_ONCE
    R_lo pScl(tmax);
    Scl *= pScl;
    iScl/= pScl;
    dt  /= pScl;
    tmax/= pScl;
    a[1]*= pScl;
#endif
  }

  series_approximation(const series_approximation & sa){
      n = sa.n;
      dt = sa.dt;
      dy = sa.dy;
      tmax = sa.tmax;
      Scl = sa.Scl;
      iScl = sa.iScl;
      RR = sa.RR;
      z = sa.z;
      c = sa.c;
      shift0 = sa.shift0;
      m = sa.m;
      a = sa.a;
      abs_a = sa.abs_a;
      err_a = sa.err_a;
      m_autoshift = sa.m_autoshift;
      rts = sa.rts;
      m_probes = sa.m_probes;
      m_dt = sa.m_dt;
      eps1 = sa.eps1;
  }

  void reset(C_hi _c, R_lo _dt, R_lo _tmax, N _m, C_lo Shift0=C_lo(0), bool autoshift=true){
      assert(_dt > 0);
      assert(_tmax > 0);
      assert(_m >= 1);

      dt = _dt;
      tmax = _tmax;
      m = _m;
      shift0 = Shift0;
      n = N(1);
      Scl = R_lo(1);
      iScl = R_lo(1);
      c=_c;
      c += shift0;
      z = c;
      rts.clear();
      m_probes.clear();

      a.resize(_m+2); abs_a.resize(_m+2); err_a.resize(_m+2);
      a[0] = c_lo(c);
      a[1] = C_lo(1);
      eps1 = R_lo(1)+R_lo(16)*numeric_limits<R_lo>::epsilon();
      for (N k = 2; k <= m + 1; ++k)
        a[k] = C_lo(0); // includes truncation error R at end
      for (N k = 0; k <= m + 1; ++k)
        err_a[k] = R_lo(0);
      m_autoshift = autoshift;
      m_dt = dt;
#ifdef RESCALE_ONCE
      R_lo pScl(tmax);
      Scl *= pScl;
      iScl/= pScl;
      dt  /= pScl;
      tmax/= pScl;
      a[1]*= pScl;
#endif
  }

  //Probes based Stop criterion functions---------------------------------------------------

  //Insert a probe
  void insertProbe(C_lo pdc, C_lo pdz){
      m_probes.push_back(probeClass(pdc,pdz));
  }

  //Do one iteration for all probes
  void advanceProbes(C_lo rzlo){
      for(unsigned int i = 0; i<m_probes.size(); i++){
          m_probes[i].advance(rzlo);
      }
  }

  //apply the computed iteration in advanceProbes()
  void applyProbes(){
      for(unsigned int i = 0; i<m_probes.size(); i++){
          m_probes[i].apply();
      }
  }

  //We found a root and want to shift the reference. the probes also have to be shifted.
  //this occures one time at most
  void shiftProbes(C_lo s, C_lo sz){
      for(unsigned int i = 0; i<m_probes.size(); i++){
          m_probes[i].shift(s,sz);
      }
  }

  //Check the difference between SA values and probes values
  bool areProbesOk(vector<C_lo> &an){
      for(unsigned int i = 0; i<m_probes.size(); i++){
          C_lo ex = series_dz(m_probes[i].m_dc, an);
          C_lo dez = series_dzdc(m_probes[i].m_dc, an);
          R_lo delta = norm0(m_probes[i].dif2ex(ex));
          if (delta > norm0(dez) * m_dt){
              //This is not necessary and maybe wrong. "Shaking" is sufficient.
              ////try to see if 2nd derivative is big enought
              //C_lo ddez = series_ddzddc(m_probes[i].m_dc, an);
              //if( delta > m_dt * norm0(dez + R_lo(.5) * m_dt * ddez) )
                return false;
          }
      }
      return true;
  }

  //findRootsRel()-----------------------------------------------------------------------------------------
  // Find the possible zeros at current iteration
  // the roots will be appended to rts[]
  // The precision is about numeric_limits<R_lo>::epsilon() relative to current zoom level

  //given a polynomial p(z) this function returns the polynomial q(z)==p(z+shift)
  void shiftPoly(vector<C_lo> & res, C_lo shift){
      vector<R_lo> bino(m+2);
      for(int i=0;i<=m;i++)
         bino[i]=1;
      for(int i=0; i<=m; i++){
          C_lo v(0);
          for(int j=m-i; j>=0; j--)
              v = a[j+i] * bino[j] + shift * v;
          for(int j=1; j<=m; j++)
              bino[j] += bino[j-1];
          res[i] = v;
      }
  }

  C_lo evalPoly(C_lo c){
      C_lo pc(0),ck(1);
      for(N k = 0; k<=m; k++){
          pc += a[k] * ck; ck *=c;
      }
      return pc;
  }

  C_lo evalDeriv(C_lo c){
      C_lo pc(0),ck(1);
      for(N k = 1; k<=m; k++){
          pc += R_lo(k) * a[k] * ck; ck *=c;
      }
      return pc;
  }

  bool solveNewton(C_lo &c, R_lo cr){//Use Newton method to find root. It is supposed that we are sure there is one and only one root near c. That is the outcome won't be far away from c.
      C_lo lc(c);
      for(N k = 0; k < 50; k++){
          C_lo nc = c - evalPoly(c) / evalDeriv(c);
          if ( max(std::abs(c.re - nc.re), std::abs(c.im - nc.im)) == R_lo(0)) break;// One would ask why not test for |p(c)/p'(c)|<eps ;o)
          c = nc;
      }
      //return true;
      //verify if it is inside the current square (there may be duplicates)
      cr+=cr*numeric_limits<R_lo>::epsilon()*1024.;
      if( c.re >= lc.re - cr && c.re < lc.re + cr &&
          c.im >= lc.im - cr && c.im < lc.im + cr) return true;
      return false;
  }

  void evalPolyIA(C_lo c, R_lo cr, R_lo &rmin, R_lo &rmax){
      C_lo pc(0),ck(c);
      R_lo r0(0), r1(0);
      R_lo ac(abs(c)), ack(ac);
      R_lo acpr(ac+cr), acprk(acpr);
      pc = a[0];
      for(N k = 1; k<=m; k++){
          pc += a[k] * ck; ck *=c;
          r0 += abs_a[k] * acprk; acprk *= acpr;
          r1 += abs_a[k] * ack; ack *= ac;
      }
      r0 -= r1; //r0 = max(R_lo(0), r0);
      R_lo apc(abs(pc));
      rmin = apc - r0; //max(R_lo(0), apc - r0);
      rmax = apc + r0;
  }

  void evalDerivIA(C_lo c, R_lo cr, R_lo &rmin, R_lo &rmax){
      vector<C_lo> na(m+2);
      shiftPoly(na, c);
      C_lo pc(0);
      R_lo r(0),crk(cr);
      for(N k = 2; k<=m; k++){
          r += R_lo(k) * abs(na[k]) * crk; crk *= cr;
      }
      pc = na[1];
      R_lo apc(abs(pc));
      //cout<<"cr= "<<cr<<"; r= "<<r<<"; apc= "<<apc;
      rmin = apc - r; //max(R_lo(0), apc - r0);
      rmax = apc + r;
      //cout<<"; rmin= "<<rmin<<"; rmax= "<<rmax<<endl;
  }

  void findRootsRecurse(C_lo c, R_lo cr, N depth=0){
      //First we check if current (disk) interval is farther than tmax
      // This is not necessary but may avoid problems... to verify!
      if (abs(c)-cr > tmax) return;

      R_lo updp,dndp;
      R_lo cr1(cr * R_lo(1.5));//std::sqrt(R_lo(2)));

      evalPolyIA(c, cr1, dndp, updp);
      //if p([z]) doesn't contain 0 then there are no roots: return.
      //This is not necessary per se but it may help skip unnecessary computations... to verify!
      if (dndp > R_lo(0)) return;

      //Check if cr < eps. If true exit: no roots in the interval [c]=c+[E]*cr
      //This is temporary code: we just ignore the roots that are very close to each other.
      if( depth > N(12) ){
          //This algorithm works pretty well most of the time but for some polynomials
          //It becomes too slow because the IA is too much over-conservative and it have to
          //subdivide a LOT. That's why we restrict the depth to 12.

          //try Newton
          if( solveNewton( c, cr ))
                rts.push_back(rootInfoClass(n,2,c * Scl));//the root is scaled back that's why the * Scl
          return;
      }

      //Compute p(c) and |p'([c])|
      R_lo apc(abs(evalPoly(c)));
      evalDerivIA(c, cr1, dndp, updp);

      //if |p(c)| > cr * up(|p'([c])|) then exit: no roots in the interval
      if (apc >= cr1 * updp) return;

      //else if |p(c)| > cr * down(|p'([c])|) then subdivide
      if (apc <= cr1 * dndp) {//use Newton method
          if( solveNewton( c, cr) )
                rts.push_back(rootInfoClass(n,1,c * Scl));//the root is scaled back that's why the * Scl
      } else {
          //subdivide
          cr *= R_lo(.5);
          findRootsRecurse(c + C_lo( cr, cr), cr, depth+1);
          findRootsRecurse(c + C_lo( cr,-cr), cr, depth+1);
          findRootsRecurse(c + C_lo(-cr, cr), cr, depth+1);
          findRootsRecurse(c + C_lo(-cr,-cr), cr, depth+1);
      }
  }

  void findRootsRel(){
      //We will look for zeroes in the square centered at a[0] and whose half-width is tmax
      //Don't forget that tmax is scaled by Scl. We need to unscale the roots at the end.

      //Computing the roots would be useful only when the minibrot is big enought w.r.t. screen
      //TODO: Maybe we need to compute some properties of the hyperbolic component like:
      //-type: disc or cardioid.
      //-size: one may be able to avoid to draw the inside.
      //-size of the atom domain. seems to be related to series approx and perturbation accuracies.
      //-parents periods. May give info about the Misiurevic points when at the scale of an imbedded julia and the like?
      //TODO: Don't compute roots if they are close to previously computed ones. This is also the
      //      case when a minibrot is visible. Close roots take longer to be calculated.

      //TOERASE: First, if the current n is multiple of the period of the first found root then skip (in order to avoid hangings)
      //if (rts.size()>0)
      //    if(n % rts[0].n == 0) return;

      for(N k=-1; k<=1; k++)
          for(N l=-1; l<=1; l++)
        findRootsRecurse(C_lo(k,l)*(tmax*R_lo(2)/R_lo(3)), tmax/R_lo(3));
      //TODO: We need to filter the duplicates
  }
//-----------------------------------------------------------------------------------------

  void updateScale(){
      N expo(0);
      R_lo rr = norm0(a[m]);
      frexpl(rr, &expo);
      if(expo >= 4000){
          expo /= m+1;
          R_lo pScl(1);
          pScl = ldexpl(pScl, -expo);
          Scl *= pScl;
          iScl/= pScl;
          dt  /= pScl;
          tmax/= pScl;
          R_lo pSclpi(pScl);
          for(N i = 1; i <= m+1; i++){
              a[i] *= pSclpi;
              pSclpi *= pScl;
          }

          RR = a[m+1].real();
          calcAbsA();
      }
  }

  //This will be used a lot of times.
  void calcAbsA(){
      for(N k=0; k<= m+1; k++) abs_a[k] = abs(a[k]);
  }

#if 1
#if 0
  //Attempt to parallelization with OpenMP: faster only for "super" deep zooms(up to 2x speedup). Maybe slower than non pallalelized otherwise.
  void calcCoefficients(vector<C_lo> &a_next){
        a_next[1] = R_lo(2) * a[0] * a[1] + Scl;//C_lo(1);//Adding Scl instead of 1 because of the scaling technique.
        for (N k = 2; k <= m; k++){
            C_lo sum(0);
            N l = (k-1)/2;
            for (N i = 0; i<= l; i++){
                C_lo ttt = a[i] * a[k-i];
                sum += ttt;
            }
            sum *= R_lo(2);
            if ((k & 1) == 0){
                C_lo ttt = a[k/2] * a[k/2];
                sum += ttt;
            }
            a_next[k] = sum;
        }
    }
  // advance if series approximation is valid.
  // Uses probes.
  bool step(bool manualSkip)
  {
      C_hi z_next;
      N lastRootsFound;
      //bool shiftFailed = false;
      //bool tooMuchRoots = false;
      calcAbsA();
      #pragma omp parallel sections
      {//To do in parallel
          #pragma omp section
          {//section 1 ... humm not possible when shift occures because c is updated while z_next is computed in parallel! any solution?

              //Find roots and recenter at the first root found if necessary
                    lastRootsFound = rts.size();
                    findRootsRel();
                    //if (N(rts.size()) > 2*m) tooMuchRoots = true;//return false;//safeguard!!! TODO: verify that the probes don't get glitched and/or escape or havin infinities?

                    //add the roots found as new probes
                    if(!(lastRootsFound == 0 && rts.size()==1))
                    for (unsigned int i = lastRootsFound; i < rts.size(); i++){
                        C_lo zz = rts[i].RPos + m_dt;//shake it a little by adding m_dt
                        C_lo sz = series_dz(zz);
                        insertProbe(zz, sz);
                    }
          }

      ///#pragma omp parallel sections
      ///{//To do in parallel
          /*{//sections 2&3. Split the computation of z_next into real and imaginary parts

            // Next iterate
            z_next = (z * z + c);
          }*/
          #pragma omp section
          {//section 2
              z_next.re = (z.re + z.im) * (z.re - z.im) + c.re;
          }
          #pragma omp section
          {//section 3
              z_next.im = R_lo(2.) * z.re * z.im + c.im;
          }

      }//synchronization happens here

    /*if(shiftFailed){
        cout << "shift failed\n";
        return false;
    }*/

#if 1
      if(m_autoshift){
          N RootsFound = rts.size();
          if(lastRootsFound==0 && RootsFound==1){
              C_lo shift = rts[0].RPos; cout << "shift0 = " << shift << endl;
              shiftProbes(shift, series_dz(shift));
              vector<C_lo> na(m+2);
              shiftPoly(na, shift * iScl);//the "*iScl" is because the coefficients of a and na are scaled. Never forget that!

              if(!areProbesOk(na)){
                  cout << "shift failed\n";
                  return false;
                  //shiftFailed = true;
              } else {
                  shift0 = shift;
                  c += shift0;
                  tmax += abs(shift0)*iScl;//The new disc have to cover the one that previousely was covering the rendered area.
                  z  = na[0]; cout << "z = " << na[0] << endl;
                  na[m+1] = a[m+1];
                  a = na;
                  rts[0].RPos=C_lo(0);
                  calcAbsA();
                  z_next = (z * z + c);//because if we do a shift the z_next computed above is no longuer valid

              }
          }
      }
#endif
    // update scale. Fixme: tmax^(m+1) should never become infinity.
    updateScale();


    // calculate coefficients
    vector<C_lo> a_next(m + 2);
    C_lo zlo = c_lo(z_next);
    if (norm(zlo) > R_lo(4))
        return false;
    advanceProbes(a[0]);//--------------------
    calcCoefficients(a_next);
    a_next[0] = zlo;

    //Temporary fix to avoid crashes when using user provided skip iteration number
    //---> tmax should not be too big: tmax^(m+1) < max long float ... or something like that.
    //if(std::isnan(sum) || std::isinf(sum)) return false;
    if(tmax > ldexpl(1,8000/(m+1))) return false;

    // check validity of next
    bool valid = areProbesOk(a_next);//isValid(a_next);

    // advance
    if (valid || manualSkip)
    {
      n = n + 1;
      z = z_next;
      a = a_next;
      applyProbes();
    }
    return valid || manualSkip;
  }
#else
  void calcCoefficients(vector<C_lo> &a_next, const C_lo &z_next){
      a_next[0] = z_next;
      a_next[1] = R_lo(2) * a[0] * a[1] + Scl;//C_lo(1);//Adding Scl instead of 1 because of the scaling technique.
      for (N k = 2; k <= m; k++){
          C_lo sum(0);
          N l = (k-1)/2;
          for (N i = 0; i<= l; i++){
              C_lo ttt = a[i] * a[k-i];
              sum += ttt;
          }
          sum *= R_lo(2);
          if ((k & 1) == 0){
              C_lo ttt = a[k/2] * a[k/2];
              sum += ttt;
          }
          a_next[k] = sum;
      }
  }
  // advance if series approximation is valid.
  // Uses probes.
  bool step(bool manualSkip)
  {
      calcAbsA();
      //Find roots and recenter at the first root found if necessary
            N lastRootsFound = rts.size();
            findRootsRel();
            if(N(rts.size()) > 2*m) return false;
            //add the roots found as new probes
            if(!(lastRootsFound == 0 && rts.size()==1))
            for (unsigned int i = lastRootsFound; i < rts.size(); i++){
                C_lo zz = rts[i].RPos + m_dt;//shake it a little by adding m_dt
                C_lo sz = series_dz(zz);
                insertProbe(zz, sz);
            }
      #if 1
            if(m_autoshift){
                N RootsFound = rts.size();
                if(lastRootsFound==0 && RootsFound==1){
                    C_lo shift = rts[0].RPos; cout << "shift0 = " << shift << endl;
                    shiftProbes(shift, series_dz(shift));
                    vector<C_lo> na(m+2);
                    shiftPoly(na, shift * iScl);//the "*iScl" is because the coefficients of a and na are scaled. Never forget that!

                    if(!areProbesOk(na)){
                        cout << "shift failed\n";
                        return false;
                    }
                    shift0 = shift;
                    c += shift0;
                    tmax += abs(shift0)*iScl;//The new disc have to cover the one that previousely was covering the rendered area.
                    z  = na[0]; cout << "z = " << na[0] << endl;
                    na[m+1] = a[m+1];
                    a = na;
                    rts[0].RPos=C_lo(0);
                    calcAbsA();
                }
            }
      #endif

    // Next iterate
	C_hi z_next(z * z + c);
	
    // update scale. Fixme: tmax^(m+1) should never become infinity.
#ifndef RESCALE_ONCE
    updateScale();
#endif


    // calculate coefficients
    vector<C_lo> a_next(m + 2);
    C_lo zlo = c_lo(z_next);
    if (norm(zlo) > R_lo(4))
        return false;
    advanceProbes(a[0]);//--------------------
    calcCoefficients(a_next, zlo);

    //Temporary fix to avoid crashes when using user provided skip iteration number
    //---> tmax should not be too big: tmax^(m+1) < max long float ... or something like that.
    //if(std::isnan(sum) || std::isinf(sum)) return false;
    if(tmax > ldexpl(1,8000/(m+1))) return false;

    // check validity of next
    bool valid = areProbesOk(a_next);//isValid(a_next);

    // advance
    if (valid || manualSkip)
    {
      n = n + 1;
      z = z_next;
      a = a_next;
      applyProbes();
    }
    return valid || manualSkip;
  }
#endif
#else
  R_lo calcTruncErr(){
      unsigned int FPU_CW = GetFPUCW();
      FPURoundUp();
      R_lo sum(0);
      for (N k = m+1; k>= 0 ; k--){//here doing the computation with complex numbers then taking the absolute value will give slightly(!) better estimate for R.
          sum *=tmax;//
          R_lo sum1(0);
          N i(0);
          for(; 2*i < m+1-k ; i++)
              sum1 += abs_a[i+k] * abs_a[m+1-i];
          sum1 *= R_lo(2);
          if (2*i == m+1-k)
              sum1 += abs_a[k+i] * abs_a[k+i];
          sum+=sum1;
      }
      RestoreFPUCW(FPU_CW);
      return sum;
  }

  //Check for validity of the given series approximation
  bool isValid(const vector<C_lo>& a_next){
    // if we suppose that the magnitude of the derivative is most of the time > 1... we also suppose that the areas where the derivative is small are much smaller than a pixel.
    // in this case we need that Deltay/y almost equal to Deltax/x
          // check validity of next
      vector<R_lo> a_abs(m + 2);
      for (N k = 0; k <= m + 1; ++k)
        a_abs[k] = abs(a_next[k]);
      bool valid;
      R_lo DeltaY(a_abs[m + 1]);
      R_lo SAmax(0);
      //R_lo err_a_sum(0);//------------
      R_lo tmaxn(tmax);
      for (N k = 1; k <= m; k++)
      {
          SAmax += tmaxn * a_abs[k];
          //err_a_sum += tmaxn * err_a[k];//------------
          tmaxn *= tmax;
      }
      DeltaY *= tmaxn;

      //Temporary fix: just trying to guess and quantify the effects of rounding errors
      //This works well but a little bit conservative:
      DeltaY *= expm1( numeric_limits<R_lo>::epsilon() * /*std::sqrt*/( R_lo(n * m * m) ))/numeric_limits<R_lo>::epsilon();

      valid = DeltaY < SAmax /* ldexpl(R_lo(1.),30)*/ * numeric_limits<R_lo>::epsilon() * dt;//R_lo(1024) means we can loose 10 bits of precision. In general we can give up to 20 bits (sometimes more) of precision. This is not really critical.
      return valid;
  }
  // advance if series approximation is valid
  // Truncation error based method using midpoint-radius interval (a sloppy one)
  bool step(bool manualSkip)
  {
      calcAbsA();
      //Find roots and recenter at the first root found if necessary
      //TODO: verify w.r.t. the other functions as: series-dz()...etc.
            N lastRootsFound = rts.size();
            findRootsRel();
      #if 1
            if(m_autoshift){
                N RootsFound = rts.size();
                if(lastRootsFound==0 && RootsFound==1){
                    shift0 = rts[0].RPos; cout << "shift0 = " << shift0 << endl;
                    vector<C_lo> na(m+2);
                    shiftPoly(na, shift0 * iScl);//the "*iScl" is because the coefficients of a and na are scaled. Never forget that!
                    c += shift0;
                    tmax += abs(shift0)*iScl;//The new disc have to cover the one that previousely was covering the rendered area.
                    z  = na[0]; cout << "z = " << na[0] << endl;
                    na[m+1] = a[m+1];
                    a = na;
                    rts[0].RPos=C_lo(0);
                    calcAbsA();
                }
            }
      #endif

    // Next iterate
    C_hi z_next(z * z + c);

    // update scale. Fixme: tmax^(m+1) should never become infinity.
#ifndef RESCALE_ONCE
    updateScale();
#endif


    // calculate coefficients
    vector<C_lo> a_next(m + 2);
    C_lo zlo = c_lo(z_next);
    if (norm(zlo) > R_lo(4))
        return false;
    calcCoefficients(a_next, zlo);

    // calculate truncation error
    // knighty's hyper complicated :o) formula.
    R_lo sum(0);
    sum = calcTruncErr();
    a_next[m + 1] = C_lo(sum);


    //Temporary fix to avoid crashes when using user provided skip iteration number
    //---> tmax should not be too big: tmax^(m+1) < max long float ... or something like that.
    //if(std::isnan(sum) || std::isinf(sum)) return false;
    if(tmax > ldexpl(1,8000/(m+1))) return false;

    // check validity of next
    bool valid = isValid(a_next);

    // advance
    if (valid || manualSkip)
    {
      n = n + 1;
      z = z_next;
      a = a_next;
      RR = sum;//
    }
    return valid || manualSkip;
  }
#endif

  // keep stepping until invalid
  // TODO FIXME check escape?
  // ToDo: move it to MB class in order to implement abort functionnality.
  void run(N skip, N &counter, const MB::requests & request)
  {
	bool manualSkip = skip >= 0;
    if(manualSkip)
        while (step( manualSkip) && n<skip && request != MB::R_abort)
            counter = n;
	else
        while (step( manualSkip) && request != MB::R_abort)
            counter = n;
    counter = n;
#if 1
    R_lo tmaxn(1);
        cout.precision(16);
        for (N k = 0; k <= m; k++)
        {
            cout << "a" << k << " = " << /*tmaxn */ a[k] << endl;
            tmaxn *= tmax;
        }
        cout << "R " << " = " << RR /* tmaxn*/ << endl;
        cout << "tmax = " << tmax << endl;
#endif
#if 0
        cout.precision(6);
        for (unsigned k = 0; k < rts.size(); k++)
        {
            cout << k << " : " << rts[k].n << " : " << rts[k].mult << " : " << rts[k].Scl << " : " << rts[k].RPos << endl;
        }
        cout << "Roots found = " << rts.size() << endl;
#endif
  }

  // initialize a pixel
  C_lo series_dz(const C_lo &dc1)
  {
    C_lo dc11((dc1) * iScl);
    C_lo sum(a[m]);
    for (N k = m-1; k >= 1; k--)
    {
      sum = sum * dc11 + a[k];
    }
    return sum*dc11;
  }

  // get SA using the newly computed series
  C_lo series_dz(const C_lo &dc1, const vector<C_lo> &an)
  {
    C_lo dc11((dc1) * iScl);
    C_lo sum(an[m]);
    for (N k = m-1; k >= 1; k--)
    {
      sum = sum * dc11 + an[k];
    }
    return sum*dc11;
  }

  // Get SA derivative
  C_lo series_dzdc(const C_lo &dc1)
  {
      C_lo dc11((dc1) * iScl);
      C_lo sum(R_lo(m)*a[m]);
      for (N k = m-1; k >= 1; k--)
      {
        sum = sum * dc11 + C_lo(k) * a[k];
      }
      return sum * iScl;
  }

  C_lo series_dzdc(const C_lo &dc1, const vector<C_lo> &an)
  {
    C_lo dc11((dc1) * iScl);
    //C_lo dcn(1);
    C_lo sum(R_lo(m)*an[m]);
    for (N k = m-1; k >= 1; k--)
    {
      sum = sum * dc11 + C_lo(k) * an[k];
    }
    return sum * iScl;
  }

  // Get SA 2nd derivative
  C_lo series_ddzddc(const C_lo &dc1, const vector<C_lo> &an)
  {
    C_lo dc11((dc1) * iScl);
    C_lo sum(R_lo(m)*R_lo(m-1)*an[m]);
    for (N k = m-1; k >= 2; k--)
    {
      sum = sum * dc11 + C_lo(k) * C_lo(k-1) * an[k];
    }
    return sum * iScl * iScl;
  }

  C_lo series_ddzddc(const C_lo &dc1)
  {
    C_lo dc11((dc1) * iScl);
    C_lo sum(R_lo(m)*R_lo(m-1)*a[m]);
    for (N k = m-1; k >= 2; k--)
    {
      sum = sum * dc11 + C_lo(k) * C_lo(k-1) * a[k];
    }
    return sum * iScl * iScl;
  }
  
  //SA truncation error at dc
  //Works with the "midpoint-radius" SA variation. see the deactivated step() memeber function.
  R_lo SAerr(const C_lo &_t)
  {
	C_lo t(_t * iScl);
	R_lo R(RR);
	R_lo abstn(1);
	R_lo abst = abs(t);
	C_lo tn(1);
	C_lo sum(0);
	//Compute derivative at t. Fixme: duplicate code. Use series_dzdc() instead.
	for(N i=1; i<=m; i++){
		sum = sum + R_lo(i) * a[i] * tn;
		tn *=t;
		abstn *=abst;
	}
#define TOW_POW_PMAX64 5.4210108624275221700372640043497e-20
#define TOW_POW_PMIN 0.000000059604644775390625
	R_lo rhs = (abs(sum) - R_lo(m+1) * R * abstn) * dt;//subtract the truncation part of the derivative: more conservative.
	if(rhs <= 0.) return 256.;
    R_lo lhs = R * abstn;// * abst;// Works better ???

    //Temporary fix: just trying to guess the effect of rounding errors
    //This works well but a little bit conservative:
    //lhs *= expm1( numeric_limits<R_lo>::epsilon() * /*std::sqrt*/( R_lo(n * m * m) ))/numeric_limits<R_lo>::epsilon();
    //lhs *= R_lo(n * m * m);//Works well for relatively small values of n. When n is big (>1000000) it is too weak

#if 0
    R_lo dy = abs(series_dzdc(_t));
    R_lo y  = abs(series_dz(_t));
    return sin(10. * -logl(TOW_POW_PMAX64 * y / (abs(_t) * dy)))*.5+.5;
#endif
    return lhs/rhs;
  }
  
  //Get a reference to carry on after series approximation initialisation
  reference get_reference()
  {
    return reference(n, z, c);
  }

  //Get a reference at delta. Used for glitch correction.
  //delta is relative to the current SA center.
  reference get_reference(C_lo delta)
  {
      C_hi rc(c), rz(z);
      C_lo zn = series_dz(delta);
      rc += delta; rz += zn;
      return reference(n, rz, rc);
  }

  //Get a new SA that is centered at delta.
  //delta is relative to the current SA center.
  series_approximation GetSAat(C_lo delta){//ToDo: finish it and test it.
      series_approximation SA(*this);
      vector<C_lo> na(m+2);
      SA.shiftPoly(na, delta * iScl);
      SA.shift0 += delta;
      SA.c += delta;
      SA.tmax += abs(delta)*iScl;//The new disc have to cover the one that previousely was covering the rendered area.
      SA.z  = na[0]; cout << "z = " << na[0] << endl;
      na[m+1] = a[m+1];
      SA.a = na;
      return SA;
  }

};

//--------------------------------------------------------------------
//Drawing functions
//--------------------------------------------------------------------
void drawRoots(MBImage &img, R_lo sradius, const series_approximation &p)
{
    N nbrRoots = p.rts.size();
    if(nbrRoots == 0) return;
    R_lo x,y;
    N px,py;
    for(N i=0; i<nbrRoots; i++){
        x = p.rts[i].RPos.re + p.shift0.re;//taking into account that the series approximation is not necessarily at the center
        y = p.rts[i].RPos.im + p.shift0.im;
        x /= R_lo(2) * sradius; y /= R_lo(2) * sradius;
        px = N(x * img.height + img.width * 0.5);
        py = N((y + 0.5) * img.height) ;
        //cout << i << " : " << px << " : " << py << endl;
        for(N k=-5; k<=5; k++){
            img.plot1(px+k,py+k,255,0,0);
            img.plot1(px+k,py-k,255,0,0);
        }
    }
}

// wrap around image plot to plot distance estimator colouring
void plot(MBImage &i, const perturbation &p, R_lo dt)
{
  N g(N(0.1*p.n) & 255);
  if (! p.escaped && ! p.glitched)
    g = 255;
  if (p.escaped)
  {
#ifdef DE_COLORING
    R_lo de(R_lo(2) * abs(p.dz) * log(abs(p.dz)) / abs(p.dzdc));
    g = 255 * tanh(.5 * de / dt);
#else
#endif
  }
  //i.plot(p.x, p.y, p.glitched ? 255 : g, g, g);
  if(p.glitched){
      R_lo err = 255. * 0.5 * (1. - std::sin(std::log(p.m_sa)));
      i.plot(p.x, p.y, err , p.m_sn % 256, 0);
  }
  else
    i.plot(p.x, p.y, g, g, g);
}

void plot(MBImage &i, const perturbation &p, R_lo dt, R_lo error)//Fixme: duplicate code!
{
  N g(N(0.1*p.n) & 255);
  if (! p.escaped && ! p.glitched)
    g = 255;
  if (p.escaped)
  {
#ifdef DE_COLORING
    R_lo de(R_lo(2) * abs(p.dz) * log(abs(p.dz)) / abs(p.dzdc));
    g = 255 * tanh(.5 * de / dt);
#else
#endif
  }
  if(error>1.)
      error=1.;
  error*=0.5;
  N r(error * 0 + (1.-error) * g);
  N gg(error * 255 + (1.-error) * g);
  N b(error * 0 + (1.-error) * g);
  
  if(p.glitched){
      R_lo err = 255. * 0.5 * (1. - std::sin(std::log(p.m_sa)));
      i.plot(p.x, p.y, err , p.m_sn % 256, 0);
  }
  else
	i.plot(p.x, p.y, r, gg, b);
}

//--------------------------------------------------------------------
// MB class member functions
//--------------------------------------------------------------------

//--------------------------------------------------------------------
//Glitch correction:
//Temporary solution that just works.
//ToDo: A glitch detection and correction that uses the screen, not the glitches vector
//--------------------------------------------------------------------

N MB::solve_glitches(R_lo dt, vector<perturbation> &glitches)
{
    C_lo curRefPosWrtSA(0,0);
    m_glitchPasses = 0;
    N prevGlchNbr = glitches.size();
    std::srand(prevGlchNbr);
    while(m_glitchPasses < m_maxGlitchPasses){
        if(m_request == MB::R_abort) return prevGlchNbr;

        m_glitchNbr = glitches.size();
        m_glitchPasses ++;
        //Find "best" reference... Just a random one! :)
        N bg(double(std::rand() * m_glitchNbr )/double(RAND_MAX+1));
        bg = std::min(m_glitchNbr-1, std::max(0,bg));

        //dc is the position of the next reference w.r.t. previous one
        C_lo dc = glitches[bg].dc;

        //curRefPosWrtSA is the position of the next reference w.r.t. the series approximation centre
        curRefPosWrtSA += dc;

        //We need this in order to rebase to the new reference
        C_lo dz = mp_s->series_dz(curRefPosWrtSA);
        //next reference
        reference r = mp_s->get_reference(curRefPosWrtSA);

        N bof=0;

        r.run(m_maxiters,bof,m_request);//todo: remove bof

        vector<perturbation> subglitches;
        N minIterNbr = m_minIterNbr;
        N maxIterNbr = m_maxIterNbr;
        #pragma omp parallel for schedule(static, 4) reduction(min : minIterNbr) reduction(max : maxIterNbr)
        for(unsigned k = 0; k< glitches.size(); k++)
        {
            m_glitchNbr--;

            // rebase to new reference
            const perturbation &p(glitches[k]);
            C_lo pc = Scr2C(p.x,p.y);//position of current perturbation / screen centre
            pc -= mp_s->shift0;//position of current perturbation / SA

            perturbation q(p.x, p.y, r.n0, pc - curRefPosWrtSA,
                           mp_s->series_dz(pc) - dz, mp_s->series_dzdc(pc), p.dt);

            // iterate
            switch(m_GDmethod){
              case GD_NONE:    q.runNoGD     (r, m_maxiters); break;
              case GD_PDB:     q.runPdbGD    (r, m_maxiters); break;
              case GD_KN:      q.runKnGD     (r, m_maxiters); break;
              case GD_KNPDB:   q.runKnPdbGD  (r, m_maxiters); break;
              case GD_GERRIT:  q.runGerritGD (r, m_maxiters); break;
              case GD_GERRIT_B:q.runGerritGDb(r, m_maxiters); break;
            }

            if(!q.glitched){
                minIterNbr = min(minIterNbr,q.n);
                maxIterNbr = max(maxIterNbr,q.n);
            }
            // output
            if (q.glitched){
                #pragma omp critical
                subglitches.push_back(q);
            }
            //else
                plot(*mp_img, q, dt);
        }
        m_minIterNbr = minIterNbr;
        m_maxIterNbr = maxIterNbr;
        if (subglitches.size() > 0)
            glitches = subglitches;
        else
            break;
    }
    return m_glitchNbr;
}

void MB::ComputeMB2(mpreal &sre, mpreal &sim, R_lo sradius,
              MBImage &img, N maxiters, N m, N precision, N skip)
{
    m_imageNbr = 0;
    m_minIterNbr = 0; m_maxIterNbr = 0;
    m_glitchNbr = 0; m_time = 0;
    m_glitchPasses = 0;
    // initial defaults
  N width = img.width;
  N height = img.height;

  R_lo tmax(sradius);
  precision = max(64, N(-log2l(tmax)) + 64);//set the precision == log2(radius) + Pmax (63)
  mpreal::set_default_prec(precision);//Was: mpfr_float::default_precision(precision);
  C_hi c((R_hi(sre)), (R_hi(sim)));

  R_lo dt( tmax / height);//dt == 1/2 pixel's width --> pixel's "radius"
  R_lo allowedErrorSA(R_lo(m_allowedErrorSA)*dt);
  R_lo allowedErrorGC(R_lo(m_allowedErrorGC)*dt);

  //Do series approximation
  if(mp_s==nullptr)
      mp_s = new series_approximation(c, allowedErrorSA, R_lo(1)*hypot(R_lo(1),R_lo(width)/R_lo(height))*tmax, m);
  else
      mp_s->reset(c,allowedErrorSA, R_lo(1)*hypot(R_lo(1),R_lo(width)/R_lo(height))*tmax, m);

  //set_default_prec() doesn't seem to affect the high precision numbers previousely declared/allocated
  mp_s->c.re.setPrecision(precision); mp_s->c.im.setPrecision(precision);

  //The 4 screen corner probes for SA stop criterion
  C_lo z = Scr2C(0,0); mp_s->insertProbe(z, z);
  z = Scr2C(width,0);  mp_s->insertProbe(z, z);
  z = Scr2C(width,height); mp_s->insertProbe(z, z);
  z = Scr2C(0,height); mp_s->insertProbe(z, z);

  //Run the SA
  cout << "Running the SA...\n";
  mp_s->run(skip, m_seriesNbr, m_request);

  if(m_request == MB::R_abort) return;

  //Compute remaining reference orbit
  cout << "Running the reference...\n";
  reference r(mp_s->get_reference());
  r.run(maxiters, m_referenceNbr,m_request);

  if(m_request == MB::R_abort) return;

  //this will hold the glitches
  vector<perturbation> glitches;

  //For statistics
  N minIterNbr = numeric_limits<N>::max();
  N maxIterNbr = numeric_limits<N>::min();

  //Render image
  cout << "Computing image...\n";
  #pragma omp parallel for schedule(static,2) reduction(min : minIterNbr) reduction(max : maxIterNbr)
  for (N y = 0; y < height; ++y)
  {
      if(m_request != MB::R_abort){
        for (N x = 0; x < width; ++x)
        {
          C_lo dc = Scr2C(x,y);
          dc -= mp_s->shift0;//taking into account that the series approximation is not necessarily at the center

          //R_lo error = mp_s->SAerr(dc);

          perturbation p(x, y, r.n0, dc, mp_s->series_dz(dc), mp_s->series_dzdc(dc), allowedErrorGC);

          switch(m_GDmethod){
            case GD_NONE:    p.runNoGD     (r, maxiters); break;
            case GD_PDB:     p.runPdbGD    (r, maxiters); break;
            case GD_KN:      p.runKnGD     (r, maxiters); break;
            case GD_KNPDB:   p.runKnPdbGD  (r, maxiters); break;
            case GD_GERRIT:  p.runGerritGD (r, maxiters); break;
            case GD_GERRIT_B:p.runGerritGDb(r, maxiters); break;
          }

          plot(img, p, dt,0.);//error);

          if(!p.glitched){
              minIterNbr = min(minIterNbr,p.n);
              maxIterNbr = max(maxIterNbr,p.n);
          }

    #ifdef SOLVE_GLITCHES
          if (m_solveGlitches)
          if (p.glitched)
          {
            #pragma omp critical
            glitches.push_back(p);
          }
    #endif
        }
        #pragma omp atomic
        ++m_imageNbr;
      }
  }

  m_minIterNbr = minIterNbr; m_maxIterNbr = maxIterNbr;

  if(m_request == MB::R_abort) return;

#ifdef SOLVE_GLITCHES
  cout << "Solving glitches...\n";
  if (m_solveGlitches)
  if (glitches.size() > 0)
    solve_glitches(dt, glitches);
#endif
  cout << "Finished.\n";
}

MB::~MB(){
    if(mp_locRe) delete mp_locRe;
    if(mp_locIm) delete mp_locIm;
    if(mp_s) delete mp_s;
    if(mp_r) delete mp_r;
}

void MB::setLocation(const char *sre, const char *sim){
    if(mp_locRe) delete mp_locRe;
    if(mp_locIm) delete mp_locIm;
    mp_locRe = new R_hi(sre);
    mp_locIm = new R_hi(sim);
}

void MB::setLocRel(R_lo dx, R_lo dy, R_lo zom){//(dx,dy) is assumed to be in pixel coordinates
    //make position relative to center and at correct scale
    dx = R_lo(2) * m_sradius * ((dx + R_lo(0.5)) - mp_img->width * R_lo(0.5)) / mp_img->height;
    dy = R_lo(2) * m_sradius * ((dy + R_lo(0.5)) / mp_img->height - R_lo(0.5));
    //update location
    *mp_locRe +=dx;
    *mp_locIm +=dy;
    //update zoom factor
    m_sradius *= zom;
    N precision = max(64, N(-log2l(m_sradius )) + 64);//set the precision == log2(radius) + Pmax (63)
    mp_locRe->setPrecision(precision); mp_locIm->setPrecision(precision);
}

void MB::win2tex(R_lo dx, R_lo dy, R_lo &tx, R_lo &ty){//from C plane to image coordinates
    //make position relative to center and at correct scale
    tx =  R_lo(0.5) * mp_img->height * dx / m_sradius + mp_img->width * R_lo(0.5) - R_lo(0.5) ;
    ty = (R_lo(0.5) * dy / m_sradius + R_lo(0.5)) * mp_img->height - R_lo(0.5) ;
}

void MB::setSRadius(const char *sradius){
    R_lo tmp = r_lo(sradius);
    if(tmp>4 || tmp<1e-4000L) return;
    m_sradius = tmp;
}

void MB::newMP_s(R_lo allowedErrorSA, R_lo tmax){
    if(mp_s) delete mp_s;
    mp_s = new series_approximation(C_hi(*mp_locRe,*mp_locIm), allowedErrorSA, hypot(R_lo(1),R_lo(mp_img->width)/R_lo(mp_img->height))*tmax, m_m);
}

void MB::newMP_r(){
    if(mp_r) delete mp_r;
    if(mp_s == nullptr) return;
    mp_r = new reference(mp_s->get_reference());
}

bool MB::isOK(){
    if(mp_locRe && mp_locIm && mp_img/*&& mp_s && mp_r*/) return true;
    return false;
}

int MB::worker(){
    for(;m_request != R_shutdown;){
        if(m_request==R_doWork && m_status==S_idle){
            m_status=S_busy;
            //call ComputeMB
            //if(isOK())
                ComputeMB2(*mp_locRe, *mp_locIm, m_sradius, *mp_img, m_maxiters, m_m, m_precision, m_skip);
            //
            m_request = R_nothing;
            m_status=S_finished;
        }
    }
    return 0;
}

void MB::setPrecision(N precision){
    m_precision = precision;
    mpreal::set_default_prec(precision);
}

bool MB::load(const char *fileName){
    FILE *in = fopen(fileName, "rb");
    if(!in) return false;
    char content[65536];
    N len=0;
    do {
        char b=fgetc(in);
        content[len]=b;
        len++;
        //cout << b;
    } while(! feof(in));
    fclose(in);
    //parse
    struct {N keyb, keye, valb, vale;} dict[128];
    N dlen = 0;
    dict[0].keyb = 0;
    for (N i=0; i<len; i++){
        if(content[i]==':') {content[i]=0; dict[dlen].keye = i-1; dict[dlen].valb=i+1;}
        if(content[i]=='\n') {content[i]=0; dict[dlen].vale = i-1; dlen++; dict[dlen].keyb = i+1;}
    }
    N ReInd=-1, ImInd=-1, ZoomInd=-1, IterationsInd=-1, FractalTypeInd=-1, PowerInd=-1;
    for(N i=0; i<dlen; i++){
        if (!strcmp(&content[dict[i].keyb],"Re"))
            ReInd = i;
        if (!strcmp(&content[dict[i].keyb],"Im"))
            ImInd = i;
        if (!strcmp(&content[dict[i].keyb],"Zoom"))
            ZoomInd = i;
        if (!strcmp(&content[dict[i].keyb],"Iterations"))
            IterationsInd = i;
        if (!strcmp(&content[dict[i].keyb],"FractalType"))
            FractalTypeInd = i;
        if (!strcmp(&content[dict[i].keyb],"Power"))
            PowerInd = i;
        //cout << "val = " << dict[i].keyb <<","<< dict[i].keye << ";" << dict[i].valb <<","<< dict[i].vale << endl;
        //cout << &content[dict[i].keyb] << ":: " << &content[dict[i].valb] << endl;
    }

    if(ReInd==-1 || ImInd==-1 || ZoomInd==-1 || IterationsInd==-1){
        //cout << "error: doesn't seem to be a kalles fraktaler file " << endl;
        return false;
    }
    if(FractalTypeInd != -1 && PowerInd !=-1){
        if(atoi(&content[dict[FractalTypeInd].valb])!=0){
            cout << "This is not mandelbrot set" << endl;
            return false;
        }
        if(atoi(&content[dict[PowerInd].valb])!=2){
            cout << "This is not a power 2 mandelbrot set" << endl;
            return false;
        }
    }
    //cout << "This is a power 2 mandelbrot set" << endl;

    R_lo Zoom = r_lo(&content[dict[ZoomInd].valb]);
    Zoom = R_lo(2) / Zoom; //cout << "Zoom is :" << Zoom << endl;
    setSRadius(Zoom);

    N Iterations = atoi(&content[dict[IterationsInd].valb]); //cout << "Iterations is :" << Iterations << endl;
    setMaxIter(Iterations);

    N ReLen = strlen(&content[dict[ReInd].valb]);
    N ImLen = strlen(&content[dict[ImInd].valb]);
    N precision = max(ReLen,ImLen);
    precision = N(3.322 * precision) + 64; //cout << "precision is :" << precision << " bits" << endl;
    setPrecision(precision);

    setLocation(&content[dict[ReInd].valb], &content[dict[ImInd].valb]);

    if (mp_s) {delete mp_s; mp_s = nullptr;}

    /*ifstream file(fileName);
    std::array<char, 5000> buf;
    for(; file.getline(&buf[0],5000,':'); ){
        string key(&buf[0]);
        file.getline(&buf[0],5000);
        string val(&buf[0]);
        cout << key << ":" << val << endl;
    }*/
    return true;
}

bool MB::save(const char *filename)
{
    char fName[8192];
    strcpy(fName,filename);
    int l=strlen(filename)-4;
    if(strcmp(&filename[l],".kfr"))
      strcat(fName,".kfr");
    printf("Saving to : %s\n",fName);

    FILE *out = fopen(fName, "wb");
    if(!out) return false;

    //I'm getting a problem when using string type: compiles fines but craches.
    /*char format[4] = "%Rg";

    char *sval = NULL; //mp_locRe->toString(format);
    mpfr_asprintf(&sval, format, mp_locRe->mpfr_srcptr());
    fprintf(out, "Re: %s \n", sval);

    mpfr_asprintf(&sval, format, mp_locIm->mpfr_srcptr());
    fprintf(out, "Im: %s \n", sval);*/

    mpfr_fprintf(out,"Re: %Re\r\n", mp_locRe->mpfr_srcptr());
    mpfr_fprintf(out,"Im: %Re\r\n", mp_locIm->mpfr_srcptr());

    R_lo zom = R_lo(2)/m_sradius;
    int  expo = int(std::log10(std::abs(zom)));
    fprintf(out, "Zoom: %gE%d\r\n", double(zom * std::pow(10.L, -expo)), expo);

    fprintf(out, "Iterations: %d\r\n", m_maxiters);

    fprintf(out, "IterDiv: 1.000000\r\n");
    fprintf(out, "SmoothMethod: 0\r\n");
    fprintf(out, "ColorMethod: 0\r\n");
    fprintf(out, "ColorOffset: 0\r\n");
    fprintf(out, "Rotate: 0.000000\r\n");
    fprintf(out, "Ratio: 360.000000\r\n");
    fprintf(out, "Colors: 0,0,0,41,35,190,132,225,108,214,174,82,144,73,241,241,187,233,235,179,166,219,60,135,12,62,153,36,94,13,28,6,183,71,222,179,18,77,200,67,187,139,166,31,3,90,125,9,56,37,31,93,212,203,252,150,245,69,59,19,13,137,10,28,219,174,50,32,154,80,238,64,120,54,253,18,73,50,246,158,125,73,220,173,79,20,242,68,64,102,208,107,196,48,183,50,59,161,34,246,34,145,157,225,139,31,218,176,202,153,2,185,114,157,73,44,128,126,197,153,213,233,128,178,234,201,204,83,191,103,214,191,20,214,126,45,220,142,102,131,239,87,73,97,255,105,143,97,205,209,30,157,156,22,114,114,230,29,240,132,79,74,119,2,215,232,57,44,83,203,201,18,30,51,116,158,12,244,213,212,159,212,164,89,126,53,207,50,34,244,204,207,211,144,45,72,211,143,117,230,217,29,42,229,192,247,43,120,129,135,68,14,95,80,0,212,97,141,190,123,5,21,7,59,51,130,31,24,112,146,218,100,84,206,177,133,62,105,21,248,70,106,4,150,115,14,217,22,47,103,104,212,247,74,74,208,87,104,118,250,22,187,17,173,174,36,136,121,254,82,219,37,67,229,60,244,69,211,216,40,206,11,245,197,96,89,61,151,39,138,89,118,45,208,194,201,205,104,212,73,106,121,37,8,97,64,20,177,59,106,165,17,40,193,140,214,169,11,135,151,140,47,241,21,29,154,149,193,155,225,192,126,233,168,154,167,134,194,181,84,191,154,231,217,35,209,85,144,56,40,209,217,108,161,102,94,78,225,48,156,254,217,113,159,226,165,226,12,155,180,71,101,56,42,70,137,169,130,121,122,118,120,194,99,\r\n");
    fprintf(out, "Smooth: 1\r\n");
    fprintf(out, "MultiColor: 0\r\n");
    fprintf(out, "BlendMC: 0\r\n");
    fprintf(out, "MultiColors: \r\n");
    fprintf(out, "Power: 2\r\n");
    fprintf(out, "FractalType: 0\r\n");
    fprintf(out, "Slopes: 0\r\n");
    fprintf(out, "SlopePower: 50\r\n");
    fprintf(out, "SlopeRatio: 50\r\n");
    fprintf(out, "SlopeAngle: 45\r\n");
    fprintf(out, "imag: 1\r\n");
    fprintf(out, "real: 1\r\n");
    fprintf(out, "SeedR: 0\r\n");
    fprintf(out, "SeedI: 0\r\n");
    fprintf(out, "FactorAR: 1\r\n");
    fprintf(out, "FactorAI: 0\r\n");

    fflush(out);
    fclose(out);
    return true;
}

unsigned int MB::getNbrRoots(){
        if (mp_s) return mp_s->rts.size();//mp_s->m_probes.size();//
        return 0;
}

bool MB::getRoot(unsigned int n, R_lo &x, R_lo &y, N &period, N &quality)
{
    if (n > getNbrRoots()-1) return false;
#if 1
    x = mp_s->rts[n].RPos.re + mp_s->shift0.re;//taking into account that the series approximation is not necessarily at the center
    y = mp_s->rts[n].RPos.im + mp_s->shift0.im;
#else
    x = mp_s->m_probes[n].m_dc.re + mp_s->shift0.re;//taking into account that the series approximation is not necessarily at the center
    y = mp_s->m_probes[n].m_dc.im + mp_s->shift0.im;
#endif
    period = mp_s->rts[n].n;
    quality = mp_s->rts[n].mult;
    return true;
    /*x = mp_s->rts[n].RPos.re; y = mp_s->rts[n].RPos.im;
    x /= R_lo(2) * m_sradius; y /= R_lo(2) * m_sradius;
    x = x * mp_img->height + mp_img->width * 0.5 - 0.5; //
    y = (y + 0.5) * mp_img->height - 0.5 ;              //
    return true;*/
}

bool MB::getRootF(unsigned int n, float &px, float &py)
{
    R_lo x,y;
    N p,q;
    if(getRoot( n, x, y, p, q)){
        px = (x * mp_img->height * R_lo(0.5) / (m_sradius * mp_img->width) + R_lo(0.5)) * mp_img->width - R_lo(0.5);
        py = (y * R_lo(0.5) / m_sradius + R_lo(0.5)) * mp_img->height - R_lo(0.5);
        return true;
    } else return false;
}

void MB::drawRoots()
{
    if(mp_s)
        ::drawRoots(*mp_img, m_sradius, *mp_s);
}

C_lo MB::Scr2C(int x, int y){
    return R_lo(2) * m_sradius *
           C_lo(R_lo(mp_img->width) * ((R_lo(x) + R_lo(0.5)) / R_lo(mp_img->width) - R_lo(0.5)) / R_lo(mp_img->height),
                                       (R_lo(y) + R_lo(0.5)) / R_lo(mp_img->height) - R_lo(0.5));

}
