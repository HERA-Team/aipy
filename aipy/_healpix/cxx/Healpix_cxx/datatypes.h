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
 *  This file defines various platform-independent data types.
 *  If any of the requested types is not available, compilation aborts
 *  with an error (unfortunately a rather obscure one).
 *
 *  Copyright (C) 2004 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef PLANCK_DATATYPES_H
#define PLANCK_DATATYPES_H

#include <string>
#include "message_error.h"

// Template magic to select the proper data types. These templates
// should not be used outside this file.

template <typename T, bool equalSize> struct sizeChooserHelper
  { typedef void TYPE; };

template <typename T> struct sizeChooserHelper<T,true>
  { typedef T TYPE; };

template <typename T1, typename T2, typename T3> struct sizeChooserHelper2
  { typedef T1 TYPE; };

template <typename T2, typename T3> struct sizeChooserHelper2 <void, T2, T3>
  { typedef T2 TYPE; };

template <typename T3> struct sizeChooserHelper2 <void, void, T3>
  { typedef T3 TYPE; };

template <> struct sizeChooserHelper2 <void, void, void>
  { };

template <int sz, typename T1, typename T2=char, typename T3=char>
  struct sizeChooser
  {
  typedef typename sizeChooserHelper2
    <typename sizeChooserHelper<T1,sizeof(T1)==sz>::TYPE,
     typename sizeChooserHelper<T2,sizeof(T2)==sz>::TYPE,
     typename sizeChooserHelper<T3,sizeof(T3)==sz>::TYPE >::TYPE TYPE;
  };

typedef sizeChooser<4, int, long, short>::TYPE
  int32;
typedef sizeChooser<4, unsigned int, unsigned long, unsigned short>::TYPE
  uint32;

typedef sizeChooser<8, long, long long>::TYPE
  int64;
typedef sizeChooser<8, unsigned long, unsigned long long>::TYPE
  uint64;

typedef sizeChooser<4, float, double>::TYPE
  float32;
typedef sizeChooser<8, double, long double>::TYPE
  float64;

// mapping of types to integer constants
enum { PLANCK_INT32   = 0,
       PLANCK_UINT32  = 1,
       PLANCK_INT64   = 2,
       PLANCK_UINT64  = 3,
       PLANCK_FLOAT32 = 4,
       PLANCK_FLOAT64 = 5,
       PLANCK_BOOL    = 6,
       PLANCK_STRING  = 7 };

template<typename T> struct typehelper {};

template<> struct typehelper<int32>
  { enum { id=PLANCK_INT32 }; };
template<> struct typehelper<uint32>
  { enum { id=PLANCK_UINT32 }; };
template<> struct typehelper<int64>
  { enum { id=PLANCK_INT64 }; };
template<> struct typehelper<uint64>
  { enum { id=PLANCK_UINT64 }; };
template<> struct typehelper<float32>
  { enum { id=PLANCK_FLOAT32 }; };
template<> struct typehelper<float64>
  { enum { id=PLANCK_FLOAT64 }; };
template<> struct typehelper<bool>
  { enum { id=PLANCK_BOOL }; };
template<> struct typehelper<std::string>
  { enum { id=PLANCK_STRING }; };

inline int type2size (int type)
  {
  switch (type)
    {
    case PLANCK_INT32  : return 4;
    case PLANCK_UINT32 : return 4;
    case PLANCK_INT64  : return 8;
    case PLANCK_UINT64 : return 8;
    case PLANCK_FLOAT32: return 4;
    case PLANCK_FLOAT64: return 8;
    case PLANCK_BOOL   : return 1;
    case PLANCK_STRING : return 1;
    default: throw Message_error ("unsupported data type");
    }
  }

inline int string2type(const std::string &type)
  {
  if (type=="FLOAT64") return PLANCK_FLOAT64;
  if (type=="FLOAT32") return PLANCK_FLOAT32;
  if (type=="INT32")   return PLANCK_INT32;
  if (type=="UINT32")  return PLANCK_UINT32;
  if (type=="INT64")   return PLANCK_INT64;
  if (type=="UINT64")  return PLANCK_UINT64;
  if (type=="BOOL")    return PLANCK_BOOL;
  if (type=="STRING")  return PLANCK_STRING;
  throw Message_error ("unknown data type "+type);
  }

#endif
