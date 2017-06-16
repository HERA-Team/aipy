/*
 *  This file is part of Healpix_cxx.
 *
 *  Healpix_cxx is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Healpix_cxx is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Healpix_cxx; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix, see http://healpix.jpl.nasa.gov
 */

/*
 *  Healpix_cxx is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Copyright (C) 2003, 2004, 2005 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "alm_map_tools.h"
#include "alm.h"
#include "healpix_map.h"
#include "fftpack_support.h"
#include "ylmgen.h"
#include "xcomplex.h"

using namespace std;

namespace {

void init_lam_fact_1d (int m, arr<double> &lam_fact)
  {
  for (int l=m; l<lam_fact.size(); ++l)
    lam_fact[l] = (l<2) ? 0. : 2*sqrt((2*l+1.)/(2*l-1.) * (l*l-m*m));
  }

void init_lam_fact_deriv_1d (int m, arr<double> &lam_fact)
  {
  lam_fact[m]=0;
  for (int l=m+1; l<lam_fact.size(); ++l)
    lam_fact[l] = sqrt((2*l+1.)/(2*l-1.) * (l*l-m*m));
  }

void init_normal_l (arr<double> &normal_l)
  {
  for (int l=0; l<normal_l.size(); ++l)
    normal_l[l] = (l<2) ? 0. : sqrt(1./((l+2.)*(l+1.)*l*(l-1.)));
  }

void get_chunk_info (int nrings, int &nchunks, int &chunksize)
  {
  nchunks = nrings/max(100,nrings/10) + 1;
  chunksize = (nrings+nchunks-1)/nchunks;
  }

void fill_work (const xcomplex<double> *datain, int nph, int mmax,
  bool shifted, const arr<xcomplex<double> > &shiftarr,
  arr<xcomplex<double> > &work)
  {
  for (int m=1; m<nph; ++m) work[m]=0;
  work[0]=datain[0];

  int cnt1=0, cnt2=nph;
  for (int m=1; m<=mmax; ++m)
    {
    if (++cnt1==nph) cnt1=0;
    if (--cnt2==-1) cnt2=nph-1;
    xcomplex<double> tmp = shifted ? (datain[m]*shiftarr[m]) : datain[m];
    work[cnt1] += tmp;
    work[cnt2] += conj(tmp);
    }
  }

void read_work (const arr<xcomplex<double> >& work, int nph, int mmax,
  bool shifted, const arr<xcomplex<double> > &shiftarr,
  xcomplex<double> *dataout)
  {
  int cnt2=0;
  for (int m=0; m<=mmax; ++m)
    {
    dataout[m] = work[cnt2];
    if (++cnt2==nph) cnt2=0;
    }
  if (shifted)
    for (int m=0; m<=mmax; ++m) dataout[m] *= shiftarr[m];
  }

void recalc_map2alm (int nph, int mmax, rfft &plan,
  arr<xcomplex<double> > &shiftarr)
  {
  if (plan.size() == nph) return;
  plan.Set (nph);
  double f1 = pi/nph;
  for (int m=0; m<=mmax; ++m)
    {
    if (m<nph)
      shiftarr[m].Set (cos(m*f1),-sin(m*f1));
    else
      shiftarr[m]=-shiftarr[m-nph];
    }
  }

template<typename T> void fft_map2alm (int nph, int mmax, bool shifted,
  double weight, rfft &plan, T *mapN, T *mapS,
  xcomplex<double> *phas_n, xcomplex<double> *phas_s,
  const arr<xcomplex<double> > &shiftarr, arr<xcomplex<double> > &work)
  {
  for (int m=0; m<nph; ++m) work[m] = mapN[m]*weight;
  plan.forward_c(work);
  read_work (work, nph, mmax, shifted, shiftarr, phas_n);
  if (mapN!=mapS)
    {
    for (int m=0; m<nph; ++m) work[m] = mapS[m]*weight;
    plan.forward_c(work);
    read_work (work, nph, mmax, shifted, shiftarr, phas_s);
    }
  else
    for (int m=0; m<=mmax; ++m) phas_s[m]=0;
  }

void recalc_alm2map (int nph, int mmax, rfft &plan,
  arr<xcomplex<double> > &shiftarr)
  {
  if (plan.size() == nph) return;
  plan.Set (nph);
  double f1 = pi/nph;
  for (int m=0; m<=mmax; ++m)
    {
    if (m<nph)
      shiftarr[m].Set (cos(m*f1),sin(m*f1));
    else
      shiftarr[m]=-shiftarr[m-nph];
    }
  }

template<typename T> void fft_alm2map (int nph, int mmax, bool shifted,
  rfft &plan, T *mapN, T *mapS, xcomplex<double> *b_north,
  xcomplex<double> *b_south, const arr<xcomplex<double> > &shiftarr,
  arr<xcomplex<double> > &work)
  {
  fill_work (b_north, nph, mmax, shifted, shiftarr, work);
  plan.backward_c(work);
  for (int m=0; m<nph; ++m) mapN[m] = work[m].re;
  if (mapN==mapS) return;
  fill_work (b_south, nph, mmax, shifted, shiftarr, work);
  plan.backward_c(work);
  for (int m=0; m<nph; ++m) mapS[m] = work[m].re;
  }

} // namespace

template<typename T> void map2alm (const Healpix_Map<T> &map,
  Alm<xcomplex<T> > &alm, const arr<double> &weight, bool add_alm)
  {
  planck_assert (map.Scheme()==RING, "map2alm: map must be in RING scheme");
  planck_assert (weight.size()>=2*map.Nside(),
    "map2alm: weight array has too few entries");

  int lmax = alm.Lmax(), mmax = alm.Mmax(), nside = map.Nside();

  int nchunks, chunksize;
  get_chunk_info(2*nside,nchunks,chunksize);

  arr2<xcomplex<double> > phas_n(chunksize,mmax+1), phas_s(chunksize,mmax+1);
  arr<double> cth(chunksize), sth(chunksize);
  double normfact = pi/(3*nside*nside);

  if (!add_alm) alm.SetToZero();

  for (int chunk=0; chunk<nchunks; ++chunk)
    {
    int llim=chunk*chunksize, ulim=min(llim+chunksize,2*nside);

#pragma omp parallel
{
    arr<xcomplex<double> > shiftarr(mmax+1), work(4*nside);
    rfft plan;

    int ith;
#pragma omp for schedule(dynamic,1)
    for (ith=llim; ith<ulim; ++ith)
      {
      int istart_north, istart_south, nph;
      bool shifted;
      map.get_ring_info (ith+1,istart_north,nph,cth[ith-llim],sth[ith-llim],
                         shifted);
      istart_south = 12*nside*nside - istart_north - nph;

      recalc_map2alm (nph, mmax, plan, shiftarr);
      fft_map2alm (nph, mmax, shifted, weight[ith]*normfact, plan,
        &map[istart_north], &map[istart_south], phas_n[ith-llim],
        phas_s[ith-llim], shiftarr, work);
      }
} // end of parallel region

#pragma omp parallel
{
    Ylmgen generator(lmax,mmax,1e-30);
    arr<double> Ylm;
    arr<xcomplex<double> > alm_tmp(lmax+1);
    int m;
#pragma omp for schedule(dynamic,1)
    for (m=0; m<=mmax; ++m)
      {
      for (int l=m; l<=lmax; ++l) alm_tmp[l].Set(0.,0.);
      for (int ith=0; ith<ulim-llim; ++ith)
        {
        int l;
        generator.get_Ylm(cth[ith],sth[ith],m,Ylm,l);
        if (l<=lmax)
          {
          xcomplex<double> p1 = phas_n[ith][m]+phas_s[ith][m],
                           p2 = phas_n[ith][m]-phas_s[ith][m];

          if ((l-m)&1) goto middle;
start:    alm_tmp[l].re += p1.re*Ylm[l]; alm_tmp[l].im += p1.im*Ylm[l];
          if (++l>lmax) goto end;
middle:   alm_tmp[l].re += p2.re*Ylm[l]; alm_tmp[l].im += p2.im*Ylm[l];
          if (++l<=lmax) goto start;
end:      ;
          }
        }
      xcomplex<T> *palm = alm.mstart(m);
      for (int l=m; l<=lmax; ++l)
        { palm[l].re += alm_tmp[l].re; palm[l].im += alm_tmp[l].im; }
      }
} // end of parallel region
    }
  }

template void map2alm (const Healpix_Map<float> &map,
  Alm<xcomplex<float> > &alm, const arr<double> &weight,
  bool add_alm);
template void map2alm (const Healpix_Map<double> &map,
  Alm<xcomplex<double> > &alm, const arr<double> &weight,
  bool add_alm);

template<typename T> void map2alm_iter (const Healpix_Map<T> &map,
  Alm<xcomplex<T> > &alm, int num_iter, const arr<double> &weight)
  {
  map2alm(map,alm,weight);
  for (int iter=1; iter<=num_iter; ++iter)
    {
    Healpix_Map<T> map2(map.Nside(),map.Scheme(),SET_NSIDE);
    alm2map(alm,map2);
    for (int m=0; m<map.Npix(); ++m)
      map2[m] = map[m]-map2[m];
    map2alm(map2,alm,weight,true);
    }
  }

template void map2alm_iter (const Healpix_Map<float> &map,
  Alm<xcomplex<float> > &alm, int num_iter,
  const arr<double> &weight);
template void map2alm_iter (const Healpix_Map<double> &map,
  Alm<xcomplex<double> > &alm, int num_iter,
  const arr<double> &weight);

template<typename T> void map2alm_iter2 (const Healpix_Map<T> &map,
  Alm<xcomplex<T> > &alm, double err_abs, double err_rel)
  {
  double x_err_abs=1./err_abs, x_err_rel=1./err_rel;
  arr<double> wgt(2*map.Nside());
  wgt.fill(1);
  Healpix_Map<T> map2(map);
  alm.SetToZero();
  while(true)
    {
    map2alm(map2,alm,wgt,true);
    alm2map(alm,map2);
    double errmeasure=0;
    for (int m=0; m<map.Npix(); ++m)
      {
      double err = abs(map[m]-map2[m]);
      double rel = (map[m]!=0) ? abs(err/map[m]) : 1e300;
      errmeasure = max(errmeasure,min(err*x_err_abs,rel*x_err_rel));
      map2[m] = map[m]-map2[m];
      }
cout << "map error measure: " << errmeasure << endl;
    if (errmeasure<1) break;
    }
  }

template void map2alm_iter2 (const Healpix_Map<double> &map,
  Alm<xcomplex<double> > &alm, double err_abs, double err_rel);

#define SETUP_LAMBDA \
  const double t1  = lam_lm1m*lam_fact[l]; \
  const double a_w = (l-m2)*two_on_s2 + l*(l-1); \
  const double a_x = twocth*(l-1)*lam_lm; \
  const double lambda_w = a_w*lam_lm - t1*c_on_s2; \
  const double lambda_x = m_on_s2 * (a_x-t1);

#define MAP2ALM_POL_MACRO(T1,Q1,Q2,U1,U2) \
  { \
  double lam_lm1m=lam_lm; \
  lam_lm=Ylm[l]; \
  alm_tmp[l][0].re += T1.re*lam_lm; alm_tmp[l][0].im += T1.im*lam_lm; \
  SETUP_LAMBDA \
  alm_tmp[l][1].re += Q1.re*lambda_w - U2.im*lambda_x; \
  alm_tmp[l][1].im += Q1.im*lambda_w + U2.re*lambda_x; \
  alm_tmp[l][2].re += U1.re*lambda_w + Q2.im*lambda_x; \
  alm_tmp[l][2].im += U1.im*lambda_w - Q2.re*lambda_x; \
  }

template<typename T> void map2alm_pol
  (const Healpix_Map<T> &mapT,
   const Healpix_Map<T> &mapQ,
   const Healpix_Map<T> &mapU,
   Alm<xcomplex<T> > &almT,
   Alm<xcomplex<T> > &almG,
   Alm<xcomplex<T> > &almC,
   const arr<double> &weightT,
   const arr<double> &weightQ,
   const arr<double> &weightU,
   bool add_alm)
  {
  planck_assert (mapT.Scheme()==RING,
    "map2alm_pol: maps must be in RING scheme");
  planck_assert (mapT.conformable(mapQ) && mapT.conformable(mapU),
    "map2alm_pol: maps are not conformable");
  planck_assert (almT.conformable(almG) && almT.conformable(almC),
    "map2alm_pol: a_lms are not conformable");
  planck_assert ((weightT.size()>=2*mapT.Nside()) &&
    (weightQ.size()>=2*mapT.Nside()) && (weightU.size()>=2*mapT.Nside()),
    "map2alm_pol: at least one weight array has too few entries");

  int lmax = almT.Lmax(), mmax = almT.Mmax(), nside = mapT.Nside();

  arr<double> normal_l (lmax+1);
  init_normal_l (normal_l);

  int nchunks, chunksize;
  get_chunk_info(2*nside,nchunks,chunksize);

  arr2<xcomplex<double> > phas_nT(chunksize,mmax+1),phas_sT(chunksize,mmax+1),
                          phas_nQ(chunksize,mmax+1),phas_sQ(chunksize,mmax+1),
                          phas_nU(chunksize,mmax+1),phas_sU(chunksize,mmax+1);

  arr<double> cth(chunksize), sth(chunksize);
  double normfact = pi/(3*nside*nside);

  if (!add_alm)
    { almT.SetToZero(); almG.SetToZero(); almC.SetToZero(); }

  for (int chunk=0; chunk<nchunks; ++chunk)
    {
    int llim=chunk*chunksize, ulim=min(llim+chunksize,2*nside);

#pragma omp parallel
{
    arr<xcomplex<double> > shiftarr(mmax+1), work(4*nside);
    rfft plan;

    int ith;
#pragma omp for schedule(dynamic,1)
    for (ith=llim; ith<ulim; ++ith)
      {
      int istart_north, istart_south, nph;
      bool shifted;
      mapT.get_ring_info (ith+1,istart_north,nph,cth[ith-llim],sth[ith-llim],
                          shifted);
      istart_south = 12*nside*nside - istart_north - nph;

      recalc_map2alm (nph, mmax, plan, shiftarr);
      fft_map2alm (nph, mmax, shifted, weightT[ith]*normfact, plan,
        &mapT[istart_north], &mapT[istart_south], phas_nT[ith-llim],
        phas_sT[ith-llim], shiftarr, work);
      fft_map2alm (nph, mmax, shifted, weightQ[ith]*normfact, plan,
        &mapQ[istart_north], &mapQ[istart_south], phas_nQ[ith-llim],
        phas_sQ[ith-llim], shiftarr, work);
      fft_map2alm (nph, mmax, shifted, weightU[ith]*normfact, plan,
        &mapU[istart_north], &mapU[istart_south], phas_nU[ith-llim],
        phas_sU[ith-llim], shiftarr, work);
      }
} // end of parallel region

#pragma omp parallel
{
    Ylmgen generator(lmax,mmax,1e-30);
    arr<double> Ylm;
    arr<double> lam_fact(lmax+1);
    arr<xcomplex<double>[3] > alm_tmp(lmax+1);
    int m;
#pragma omp for schedule(dynamic,1)
    for (m=0; m<=mmax; ++m)
      {
      init_lam_fact_1d (m,lam_fact);

      for (int l=m; l<alm_tmp.size(); ++l)
        alm_tmp[l][0]=alm_tmp[l][1]=alm_tmp[l][2] = 0;

      for (int ith=0; ith<ulim-llim; ++ith)
        {
        int l;
        generator.get_Ylm(cth[ith],sth[ith],m,Ylm,l);
        if (l<=lmax)
          {
          double one_on_s2 = 1/(sth[ith]*sth[ith]);
          double c_on_s2 = cth[ith] * one_on_s2;
          double two_on_s2 = 2*one_on_s2;
          double twocth = 2*cth[ith];
          int m2 = m*m;
          double m_on_s2 = m*one_on_s2;

          xcomplex<double> T1 = phas_nT[ith][m]+phas_sT[ith][m],
                           T2 = phas_nT[ith][m]-phas_sT[ith][m],
                           Q1 = phas_nQ[ith][m]+phas_sQ[ith][m],
                           Q2 = phas_nQ[ith][m]-phas_sQ[ith][m],
                           U1 = phas_nU[ith][m]+phas_sU[ith][m],
                           U2 = phas_nU[ith][m]-phas_sU[ith][m];

          double lam_lm = 0;
          if ((l-m)&1) goto middle;
start:    MAP2ALM_POL_MACRO(T1,Q1,Q2,U1,U2)
          if (++l>lmax) goto end;
middle:   MAP2ALM_POL_MACRO(T2,Q2,Q1,U2,U1)
          if (++l<=lmax) goto start;
end:      ;
          }
        }
      xcomplex<T> *palmT=almT.mstart(m), *palmG=almG.mstart(m),
                  *palmC=almC.mstart(m);
      for (int l=m;l<=lmax;++l)
        {
        palmT[l].re += alm_tmp[l][0].re;
        palmT[l].im += alm_tmp[l][0].im;
        palmG[l].re += alm_tmp[l][1].re*normal_l[l];
        palmG[l].im += alm_tmp[l][1].im*normal_l[l];
        palmC[l].re += alm_tmp[l][2].re*normal_l[l];
        palmC[l].im += alm_tmp[l][2].im*normal_l[l];
        }
      }
} // end of parallel region
    }
  }

template void map2alm_pol
  (const Healpix_Map<float> &mapT,
   const Healpix_Map<float> &mapQ,
   const Healpix_Map<float> &mapU,
   Alm<xcomplex<float> > &almT,
   Alm<xcomplex<float> > &almG,
   Alm<xcomplex<float> > &almC,
   const arr<double> &weightT,
   const arr<double> &weightQ,
   const arr<double> &weightU,
   bool add_alm);
template void map2alm_pol
  (const Healpix_Map<double> &mapT,
   const Healpix_Map<double> &mapQ,
   const Healpix_Map<double> &mapU,
   Alm<xcomplex<double> > &almT,
   Alm<xcomplex<double> > &almG,
   Alm<xcomplex<double> > &almC,
   const arr<double> &weightT,
   const arr<double> &weightQ,
   const arr<double> &weightU,
   bool add_alm);

template<typename T> void map2alm_pol_iter
  (const Healpix_Map<T> &mapT,
   const Healpix_Map<T> &mapQ,
   const Healpix_Map<T> &mapU,
   Alm<xcomplex<T> > &almT,
   Alm<xcomplex<T> > &almG,
   Alm<xcomplex<T> > &almC,
   int num_iter,
   const arr<double> &weightT,
   const arr<double> &weightQ,
   const arr<double> &weightU)
  {
  map2alm_pol(mapT,mapQ,mapU,almT,almG,almC,weightT,weightQ,weightU);
  for (int iter=1; iter<=num_iter; ++iter)
    {
    Healpix_Map<T> mapT2(mapT.Nside(),mapT.Scheme(),SET_NSIDE),
                   mapQ2(mapT.Nside(),mapT.Scheme(),SET_NSIDE),
                   mapU2(mapT.Nside(),mapT.Scheme(),SET_NSIDE);

    alm2map_pol(almT,almG,almC,mapT2,mapQ2,mapU2);
    for (int m=0; m<mapT.Npix(); ++m)
      {
      mapT2[m] = mapT[m]-mapT2[m];
      mapQ2[m] = mapQ[m]-mapQ2[m];
      mapU2[m] = mapU[m]-mapU2[m];
      }
    map2alm_pol(mapT2,mapQ2,mapU2,almT,almG,almC,weightT,weightQ,weightU,true);
    }
  }

template void map2alm_pol_iter
  (const Healpix_Map<float> &mapT,
   const Healpix_Map<float> &mapQ,
   const Healpix_Map<float> &mapU,
   Alm<xcomplex<float> > &almT,
   Alm<xcomplex<float> > &almG,
   Alm<xcomplex<float> > &almC,
   int num_iter,
   const arr<double> &weightT,
   const arr<double> &weightQ,
   const arr<double> &weightU);
template void map2alm_pol_iter
  (const Healpix_Map<double> &mapT,
   const Healpix_Map<double> &mapQ,
   const Healpix_Map<double> &mapU,
   Alm<xcomplex<double> > &almT,
   Alm<xcomplex<double> > &almG,
   Alm<xcomplex<double> > &almC,
   int num_iter,
   const arr<double> &weightT,
   const arr<double> &weightQ,
   const arr<double> &weightU);

template<typename T> void map2alm_pol_iter2
  (const Healpix_Map<T> &mapT,
   const Healpix_Map<T> &mapQ,
   const Healpix_Map<T> &mapU,
   Alm<xcomplex<T> > &almT,
   Alm<xcomplex<T> > &almG,
   Alm<xcomplex<T> > &almC,
   double err_abs, double err_rel)
  {
  arr<double> wgt(2*mapT.Nside());
  wgt.fill(1);
  Healpix_Map<T> mapT2(mapT), mapQ2(mapQ), mapU2(mapU);
  almT.SetToZero(); almG.SetToZero(); almC.SetToZero();
  while(true)
    {
    map2alm_pol(mapT2,mapQ2,mapU2,almT,almG,almC,wgt,wgt,wgt,true);
    alm2map_pol(almT,almG,almC,mapT2,mapQ2,mapU2);
    double errmeasure=0;
    for (int m=0; m<mapT.Npix(); ++m)
      {
      double err = abs(mapT[m]-mapT2[m]);
      double rel = (mapT[m]!=0) ? abs(err/mapT[m]) : 1e300;
      errmeasure = max(errmeasure,min(err/err_abs,rel/err_rel));
      mapT2[m] = mapT[m]-mapT2[m];
      err = abs(mapQ[m]-mapQ2[m]);
      rel = (mapQ[m]!=0) ? abs(err/mapQ[m]) : 1e300;
      errmeasure = max(errmeasure,min(err/err_abs,rel/err_rel));
      mapQ2[m] = mapQ[m]-mapQ2[m];
      err = abs(mapU[m]-mapU2[m]);
      rel = (mapU[m]!=0) ? abs(err/mapU[m]) : 1e300;
      errmeasure = max(errmeasure,min(err/err_abs,rel/err_rel));
      mapU2[m] = mapU[m]-mapU2[m];
      }
cout << "map error measure: " << errmeasure << endl;
    if (errmeasure<1) break;
    }
  }

template void map2alm_pol_iter2
  (const Healpix_Map<double> &mapT,
   const Healpix_Map<double> &mapQ,
   const Healpix_Map<double> &mapU,
   Alm<xcomplex<double> > &almT,
   Alm<xcomplex<double> > &almG,
   Alm<xcomplex<double> > &almC,
   double err_abs, double err_rel);

template<typename T> void alm2map (const Alm<xcomplex<T> > &alm,
  Healpix_Map<T> &map)
  {
  int lmax = alm.Lmax(), mmax = alm.Mmax(), nside = map.Nside();

  int nchunks, chunksize;
  get_chunk_info(2*nside,nchunks,chunksize);

  arr2<xcomplex<double> > b_north(chunksize,mmax+1), b_south(chunksize,mmax+1);
  arr<double> cth(chunksize),sth(chunksize);

  for (int chunk=0; chunk<nchunks; ++chunk)
    {
    int llim=chunk*chunksize, ulim=min(llim+chunksize,2*nside);
    for (int ith=llim; ith<ulim; ++ith)
      {
      int nph, istart_north;
      bool shifted;
      map.get_ring_info (ith+1,istart_north,nph, cth[ith-llim],sth[ith-llim],
                         shifted);
      }

#pragma omp parallel
{
    Ylmgen generator(lmax,mmax,1e-30);
    arr<double> Ylm;
    arr<xcomplex<double> > alm_tmp(lmax+1);
    int m;
#pragma omp for schedule(dynamic,1)
    for (m=0; m<=mmax; ++m)
      {
      for (int l=m; l<=lmax; ++l)
        alm_tmp[l]=alm(l,m);

      for (int ith=0; ith<ulim-llim; ++ith)
        {
        int l;
        generator.get_Ylm(cth[ith],sth[ith],m,Ylm,l);
        if (l<=lmax)
          {
          xcomplex<double> p1=0, p2=0;

          if ((l-m)&1) goto middle;
start:    p1.re += alm_tmp[l].re*Ylm[l]; p1.im += alm_tmp[l].im*Ylm[l];
          if (++l>lmax) goto end;
middle:   p2.re += alm_tmp[l].re*Ylm[l]; p2.im += alm_tmp[l].im*Ylm[l];
          if (++l<=lmax) goto start;
end:      b_north[ith][m] = p1+p2; b_south[ith][m] = p1-p2;
          }
        else
          {
          b_north[ith][m] = b_south[ith][m] = 0;
          }
        }
      }
} // end of parallel region

#pragma omp parallel
{
    arr<xcomplex<double> > shiftarr(mmax+1), work(4*nside);
    rfft plan;
    int ith;
#pragma omp for schedule(dynamic,1)
    for (ith=llim; ith<ulim; ++ith)
      {
      int istart_north, istart_south, nph;
      double dum1, dum2;
      bool shifted;
      map.get_ring_info (ith+1,istart_north,nph,dum1,dum2,shifted);
      istart_south = map.Npix()-istart_north-nph;
      recalc_alm2map (nph, mmax, plan, shiftarr);
      fft_alm2map (nph, mmax, shifted, plan, &map[istart_north],
        &map[istart_south], b_north[ith-llim], b_south[ith-llim],
        shiftarr, work);
      }
} // end of parallel region
    }
  }

template void alm2map (const Alm<xcomplex<double> > &alm,
  Healpix_Map<double> &map);
template void alm2map (const Alm<xcomplex<float> > &alm,
  Healpix_Map<float> &map);

#define ALM2MAP_POL_MACRO(bT1,bQ1,bQ2,bU1,bU2) \
  { \
  double lam_lm1m = lam_lm; \
  lam_lm = Ylm[l]; \
  SETUP_LAMBDA \
  bT1.re+=alm_tmp[l][0].re*lam_lm;   bT1.im+=alm_tmp[l][0].im*lam_lm; \
  bQ1.re+=alm_tmp[l][1].re*lambda_w; bQ1.im+=alm_tmp[l][1].im*lambda_w; \
  bQ2.re-=alm_tmp[l][2].im*lambda_x; bQ2.im+=alm_tmp[l][2].re*lambda_x; \
  bU1.re-=alm_tmp[l][2].re*lambda_w; bU1.im-=alm_tmp[l][2].im*lambda_w; \
  bU2.re-=alm_tmp[l][1].im*lambda_x; bU2.im+=alm_tmp[l][1].re*lambda_x; \
  }

template<typename T> void alm2map_pol
  (const Alm<xcomplex<T> > &almT,
   const Alm<xcomplex<T> > &almG,
   const Alm<xcomplex<T> > &almC,
   Healpix_Map<T> &mapT,
   Healpix_Map<T> &mapQ,
   Healpix_Map<T> &mapU)
  {
  planck_assert (mapT.Scheme()==RING,
    "alm2map_pol: maps must be in RING scheme");
  planck_assert (mapT.conformable(mapQ) && mapT.conformable(mapU),
    "alm2map_pol: maps are not conformable");
  planck_assert (almT.conformable(almG) && almT.conformable(almC),
    "alm2map_pol: a_lms are not conformable");

  int lmax = almT.Lmax(), mmax = almT.Mmax(), nside = mapT.Nside();

  arr<double> normal_l (lmax+1);
  init_normal_l (normal_l);

  int nchunks, chunksize;
  get_chunk_info(2*nside,nchunks,chunksize);

  arr2<xcomplex<double> >
    b_north_T(chunksize,mmax+1), b_south_T(chunksize,mmax+1),
    b_north_Q(chunksize,mmax+1), b_south_Q(chunksize,mmax+1),
    b_north_U(chunksize,mmax+1), b_south_U(chunksize,mmax+1);

  arr<double> cth(chunksize),sth(chunksize);

  for (int chunk=0; chunk<nchunks; ++chunk)
    {
    int llim=chunk*chunksize, ulim=min(llim+chunksize,2*nside);
    for (int ith=llim; ith<ulim; ++ith)
      {
      int nph, istart_north;
      bool shifted;
      mapT.get_ring_info (ith+1,istart_north,nph, cth[ith-llim],sth[ith-llim],
                          shifted);
      }

#pragma omp parallel
{
    Ylmgen generator(lmax,mmax,1e-30);
    arr<double> Ylm;
    arr<double> lam_fact (lmax+1);
    arr<xcomplex<double>[3]> alm_tmp(lmax+1);
    int m;
#pragma omp for schedule(dynamic,1)
    for (m=0; m<=mmax; ++m)
      {
      int m2 = m*m;
      init_lam_fact_1d (m,lam_fact);
      for (int l=m; l<=lmax; ++l)
        {
        alm_tmp[l][0] = almT(l,m);
        alm_tmp[l][1] = almG(l,m)*(-normal_l[l]);
        alm_tmp[l][2] = almC(l,m)*(-normal_l[l]);
        }
      for (int ith=0; ith<ulim-llim; ++ith)
        {
        int l;
        generator.get_Ylm(cth[ith],sth[ith],m,Ylm,l);
        if (l<=lmax)
          {
          double one_on_s2 = 1/(sth[ith]*sth[ith]);
          double c_on_s2 = cth[ith] * one_on_s2;
          double two_on_s2 = 2*one_on_s2;
          double m_on_s2 = m*one_on_s2;
          double twocth = 2*cth[ith];

          xcomplex<double> bT1=0, bT2=0, bQ1=0, bQ2=0, bU1=0, bU2=0;

          double lam_lm = 0;
          if ((l-m)&1) goto middle;
start:    ALM2MAP_POL_MACRO(bT1,bQ1,bQ2,bU1,bU2)
          if (++l>lmax) goto end;
middle:   ALM2MAP_POL_MACRO(bT2,bQ2,bQ1,bU2,bU1)
          if (++l<=lmax) goto start;
end:      b_north_T[ith][m] = bT1+bT2; b_south_T[ith][m] = bT1-bT2;
          b_north_Q[ith][m] =-bQ1-bQ2; b_south_Q[ith][m] =-bQ1+bQ2;
          b_north_U[ith][m] = bU1+bU2; b_south_U[ith][m] = bU1-bU2;
          }
        else
          {
          b_north_T[ith][m] = b_south_T[ith][m] = 0;
          b_north_Q[ith][m] = b_south_Q[ith][m] = 0;
          b_north_U[ith][m] = b_south_U[ith][m] = 0;
          }
        }
      }
} // end of parallel region

#pragma omp parallel
{
    arr<xcomplex<double> > shiftarr(mmax+1), work(4*nside);
    rfft plan;
    int ith;
#pragma omp for schedule(dynamic,1)
    for (ith=llim; ith<ulim; ++ith)
      {
      int nph, istart_north, istart_south;
      double dum1, dum2;
      bool shifted;
      mapT.get_ring_info (ith+1,istart_north,nph,dum1,dum2,shifted);
      istart_south = mapT.Npix()-istart_north-nph;
      recalc_alm2map (nph, mmax, plan, shiftarr);
      fft_alm2map (nph, mmax, shifted, plan, &mapT[istart_north],
        &mapT[istart_south], b_north_T[ith-llim], b_south_T[ith-llim],
        shiftarr, work);
      fft_alm2map (nph, mmax, shifted, plan, &mapQ[istart_north],
        &mapQ[istart_south], b_north_Q[ith-llim], b_south_Q[ith-llim],
        shiftarr, work);
      fft_alm2map (nph, mmax, shifted, plan, &mapU[istart_north],
        &mapU[istart_south], b_north_U[ith-llim], b_south_U[ith-llim],
        shiftarr, work);
      }
} // end of parallel region
    }
  }

template void alm2map_pol (const Alm<xcomplex<double> > &almT,
                           const Alm<xcomplex<double> > &almG,
                           const Alm<xcomplex<double> > &almC,
                           Healpix_Map<double> &mapT,
                           Healpix_Map<double> &mapQ,
                           Healpix_Map<double> &mapU);

template void alm2map_pol (const Alm<xcomplex<float> > &almT,
                           const Alm<xcomplex<float> > &almG,
                           const Alm<xcomplex<float> > &almC,
                           Healpix_Map<float> &mapT,
                           Healpix_Map<float> &mapQ,
                           Healpix_Map<float> &mapU);

#define ALM2MAP_DER1_MACRO(b1,bdth1,bdph1) \
  { \
  double lam_lm1m = lam_lm; \
  lam_lm = Ylm[l]; \
  const double t1 = alm_tmp[l].re*lam_lm; \
  const double t2 = alm_tmp[l].im*lam_lm; \
  const double t3 = l*cotanth; \
  const double t4 = one_on_s*lam_lm1m*lam_fact[l]; \
  b1.re+=t1; b1.im+=t2; \
  bdth1.re+=t3*t1-t4*alm_tmp[l].re; bdth1.im+=t3*t2-t4*alm_tmp[l].im; \
  bdph1.re-=m*t2; bdph1.im+=m*t1; \
  }

template<typename T> void alm2map_der1
  (const Alm<xcomplex<T> > &alm,
   Healpix_Map<T> &map,
   Healpix_Map<T> &mapdth,
   Healpix_Map<T> &mapdph)
  {
  planck_assert (map.Scheme()==RING,
    "alm2map_der1: maps must be in RING scheme");
  planck_assert (map.conformable(mapdth) && map.conformable(mapdph),
    "alm2map_der1: maps are not conformable");

  int lmax = alm.Lmax(), mmax = alm.Mmax(), nside = map.Nside();

  int nchunks, chunksize;
  get_chunk_info(2*nside,nchunks,chunksize);

  arr2<xcomplex<double> >
    b_north(chunksize,mmax+1), b_south(chunksize,mmax+1),
    b_north_dth(chunksize,mmax+1), b_south_dth(chunksize,mmax+1),
    b_north_dph(chunksize,mmax+1), b_south_dph(chunksize,mmax+1);

  arr<double> cth(chunksize),sth(chunksize);

  for (int chunk=0; chunk<nchunks; ++chunk)
    {
    int llim=chunk*chunksize, ulim=min(llim+chunksize,2*nside);
    for (int ith=llim; ith<ulim; ++ith)
      {
      int nph, istart_north;
      bool shifted;
      map.get_ring_info (ith+1,istart_north,nph, cth[ith-llim],sth[ith-llim],
                         shifted);
      }

#pragma omp parallel
{
    Ylmgen generator(lmax,mmax,1e-30);
    arr<double> Ylm;
    arr<double> lam_fact (lmax+1);
    int m;
#pragma omp for schedule(dynamic,1)
    for (m=0; m<=mmax; ++m)
      {
      const xcomplex<T> *alm_tmp=alm.mstart(m);
      init_lam_fact_deriv_1d (m,lam_fact);
      for (int ith=0; ith<ulim-llim; ++ith)
        {
        int l;
        generator.get_Ylm(cth[ith],sth[ith],m,Ylm,l);
        if (l<=lmax)
          {
          double one_on_s = 1/sth[ith];
          double cotanth = cth[ith]*one_on_s;

          xcomplex<double> b1=0, b2=0, bdth1=0, bdth2=0, bdph1=0, bdph2=0;

          double lam_lm = 0;
          if ((l-m)&1) goto middle;
start:    ALM2MAP_DER1_MACRO(b1,bdth1,bdph1)
          if (++l>lmax) goto end;
middle:   ALM2MAP_DER1_MACRO(b2,bdth2,bdph2)
          if (++l<=lmax) goto start;
end:      b_north[ith][m] = b1+b2; b_south[ith][m] = b1-b2;
          b_north_dth[ith][m] = bdth1+bdth2; b_south_dth[ith][m] = bdth2-bdth1;
          b_north_dph[ith][m] = (bdph1+bdph2)*one_on_s;
          b_south_dph[ith][m] = (bdph1-bdph2)*one_on_s;
          }
        else
          {
          b_north[ith][m] = b_south[ith][m] = 0;
          b_north_dth[ith][m] = b_south_dth[ith][m] = 0;
          b_north_dph[ith][m] = b_south_dph[ith][m] = 0;
          }
        }
      }
} // end of parallel region

#pragma omp parallel
{
    arr<xcomplex<double> > shiftarr(mmax+1), work(4*nside);
    rfft plan;
    int ith;
#pragma omp for schedule(dynamic,1)
    for (ith=llim; ith<ulim; ++ith)
      {
      int nph, istart_north, istart_south;
      double dum1, dum2;
      bool shifted;
      map.get_ring_info (ith+1,istart_north,nph,dum1,dum2,shifted);
      istart_south = map.Npix()-istart_north-nph;
      recalc_alm2map (nph, mmax, plan, shiftarr);
      fft_alm2map (nph, mmax, shifted, plan, &map[istart_north],
        &map[istart_south], b_north[ith-llim], b_south[ith-llim],
        shiftarr, work);
      fft_alm2map (nph, mmax, shifted, plan, &mapdth[istart_north],
        &mapdth[istart_south], b_north_dth[ith-llim], b_south_dth[ith-llim],
        shiftarr, work);
      fft_alm2map (nph, mmax, shifted, plan, &mapdph[istart_north],
        &mapdph[istart_south], b_north_dph[ith-llim], b_south_dph[ith-llim],
        shiftarr, work);
      }
} // end of parallel region
    }
  }

template void alm2map_der1 (const Alm<xcomplex<double> > &alm,
                            Healpix_Map<double> &map,
                            Healpix_Map<double> &map_dth,
                            Healpix_Map<double> &map_dph);

template void alm2map_der1 (const Alm<xcomplex<float> > &alm,
                            Healpix_Map<float> &map,
                            Healpix_Map<float> &map_dth,
                            Healpix_Map<float> &map_dph);
