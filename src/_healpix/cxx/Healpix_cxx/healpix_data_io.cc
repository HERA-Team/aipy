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
 *  Copyright (C) 2003, 2005 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "healpix_data_io.h"
#include "arr.h"
#include "fitshandle.h"
#include "paramfile.h"
#include "simparams.h"

using namespace std;

void read_weight_ring (const string &dir, int nside, arr<double> &weight_T)
  {
  fitshandle inp;
  inp.open(dir+"/weight_ring_n"+intToString(nside,5)+".fits");
  inp.goto_hdu(2);
  weight_T.alloc (2*nside);
  inp.read_column(1,weight_T);
  }
void read_weight_ring (const string &dir, int nside, arr<double> &weight_T,
  arr<double> &weight_Q, arr<double> &weight_U)
  {
  fitshandle inp;
  inp.open(dir+"/weight_ring_n"+intToString(nside,5)+".fits");
  inp.goto_hdu(2);
  weight_T.alloc (2*nside); weight_Q.alloc (2*nside); weight_U.alloc (2*nside);
  inp.read_column(1,weight_T);
  weight_Q.fill(0);
  weight_U.fill(0);
  }

void get_ring_weights (paramfile &params, simparams &par, int nside,
  arr<double> &weight_T)
  {
  bool weighted = params.find<bool>("weighted",false);
  par.add ("weighted","WEIGHTED",weighted,"ring weights used?");
  weight_T.alloc (2*nside);
  if (weighted)
    {
    string datadir = params.find<string>("healpix_data");
    read_weight_ring (datadir, nside, weight_T);
    for (int m=0; m<weight_T.size(); ++m) weight_T[m]+=1;
    }
  else
    weight_T.fill(1);
  }
void get_ring_weights (paramfile &params, simparams &par, int nside,
  arr<double> &weight_T, arr<double> &weight_Q, arr<double> &weight_U)
  {
  bool weighted = params.find<bool>("weighted",false);
  par.add ("weighted","WEIGHTED",weighted,"ring weights used?");
  weight_T.alloc (2*nside); weight_Q.alloc (2*nside); weight_U.alloc (2*nside);
  if (weighted)
    {
    string datadir = params.find<string>("healpix_data");
    read_weight_ring (datadir, nside, weight_T, weight_Q, weight_U);
    for (int m=0; m<weight_T.size(); ++m) weight_T[m]+=1;
    for (int m=0; m<weight_Q.size(); ++m) weight_Q[m]+=1;
    for (int m=0; m<weight_U.size(); ++m) weight_U[m]+=1;
    }
  else
    { weight_T.fill(1); weight_Q.fill(1); weight_U.fill(1); }
  }

void read_pixwin (const string &dir, int nside, arr<double> &temp)
  {
  fitshandle inp;
  inp.open(dir+"/pixel_window_n"+intToString(nside,4)+".fits");
  inp.goto_hdu(2);
  if (temp.size()==0)
    inp.read_entire_column(1,temp);
  else
    inp.read_column(1,temp);
  }
void read_pixwin (const string &dir, int nside, arr<double> &temp,
  arr<double> &pol)
  {
  fitshandle inp;
  inp.open(dir+"/pixel_window_n"+intToString(nside,4)+".fits");
  inp.goto_hdu(2);
  if (temp.size()==0)
    inp.read_entire_column(1,temp);
  else
    inp.read_column(1,temp);
  if (pol.size()==0)
    inp.read_entire_column(2,pol);
  else
    inp.read_column(2,pol);
  }

void get_pixwin (paramfile &params, simparams &par, int lmax,
  int nside, arr<double> &pixwin)
  {
  bool do_pixwin = params.find<bool>("pixel_window",false);
  par.add("pixel_window","PIXWIN",do_pixwin,"pixel window used?");
  pixwin.alloc(lmax+1);
  pixwin.fill(1);
  if (do_pixwin)
    {
    string datadir = params.find<string>("healpix_data");
    read_pixwin (datadir,nside,pixwin);
    }
  }
void get_pixwin (paramfile &params, simparams &par, int lmax,
  int nside, arr<double> &pixwin, arr<double> &pixwin_pol)
  {
  bool do_pixwin = params.find<bool>("pixel_window",false);
  par.add("pixel_window","PIXWIN",do_pixwin,"pixel window used?");
  pixwin.alloc(lmax+1);
  pixwin.fill(1);
  pixwin_pol.alloc(lmax+1);
  pixwin_pol.fill(1);
  if (do_pixwin)
    {
    string datadir = params.find<string>("healpix_data");
    read_pixwin (datadir,nside,pixwin,pixwin_pol);
    }
  }
