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
 *  This file contains the implementation of various convenience functions
 *  used by the Planck LevelS package.
 *
 *  Copyright (C) 2002, 2003, 2004, 2005 Max-Planck-Society
 *  Authors: Martin Reinecke, Reinhard Hell
 */

// if we are using g++, check for version 3.0 or higher
#ifdef __GNUC__
#if (__GNUC__<3)
#error your C++ compiler is too old. g++ version 3.0 or higher is required.
#endif
#endif

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdio>
#include "cxxutils.h"

using namespace std;

bool file_present (const string &filename)
  {
  ifstream dummy(filename.c_str());
  return dummy;
  }

void assert_present (const string &filename)
  {
  if (file_present(filename)) return;
  throw Message_error ("Error: file " + filename + " does not exist!");
  }

void assert_not_present (const string &filename)
  {
  if (!file_present(filename)) return;
  throw Message_error ("Error: file " + filename + " already exists!");
  }

void remove_file (const string &filename)
  {
  remove (filename.c_str());
  }

string trim (const string &orig)
  {
  string::size_type p1=orig.find_first_not_of(" \t");
  if (p1==string::npos) return "";
  string::size_type p2=orig.find_last_not_of(" \t");
  return orig.substr(p1,p2-p1+1);
  }

template<typename T> string dataToString (const T &x)
  {
  ostringstream strstrm;
  strstrm << x;
  return trim(strstrm.str());
  }

template<> string dataToString (const bool &x)
  { return x ? "T" : "F"; }
template<> string dataToString (const string &x)
  { return trim(x); }
template<> string dataToString (const float &x)
  {
  ostringstream strstrm;
  strstrm << setprecision(8) << x;
  return trim(strstrm.str());
  }
template<> string dataToString (const double &x)
  {
  ostringstream strstrm;
  strstrm << setprecision(16) << x;
  return trim(strstrm.str());
  }

template string dataToString (const int &x);
template string dataToString (const unsigned int &x);
template string dataToString (const long &x);
template string dataToString (const unsigned long long &x);
template string dataToString (const long long &x);
template string dataToString (const unsigned long &x);

string intToString(int x, int width)
  {
  ostringstream strstrm;
  strstrm << setw(width) << setfill('0') << x;
  return trim(strstrm.str());
  }

template<typename T> void stringToData (const string &x, T &value)
  {
  istringstream strstrm(x);
  strstrm >> value;
  planck_assert(strstrm, "conversion error in stringToData()");
// FIXME: disable the test below for now, since some compilers choke on it
//  planck_assert(strstrm.tellg()==streampos(x.length()),
//    "parse error in stringToData()");
  }

template<> void stringToData (const string &x, string &value)
  { value = trim(x); }

template<> void stringToData (const string &x, bool &value)
  {
  if ( x=="F" || x=="f" || x=="n" || x=="N" || x=="false" || x==".false."
       || x=="FALSE" || x==".FALSE.")
    value=false;
  else if (x=="T" || x=="t" || x=="y" || x=="Y" || x=="true" || x==".true."
       || x=="TRUE" || x==".TRUE.")
    value=true;
  else throw Message_error ("stringToData<bool>: error parsing argument");
  }

template void stringToData (const string &x, int &value);
template void stringToData (const string &x, long &value);
template void stringToData (const string &x, float &value);
template void stringToData (const string &x, double &value);
template void stringToData (const string &x, unsigned long long &value);
template void stringToData (const string &x, long long &value);
template void stringToData (const string &x, unsigned long &value);
template void stringToData (const string &x, unsigned int &value);

void announce_progress (int now, int total)
  {
  if ((now%(max(total/100,1)))==0)
    cout << "\r " << setw(3) << planck_nint ((now*100.)/total)
         << "% done\r" << flush;
  }

void announce_progress (double now, double last, double total)
  {
  int lastpercent = int((last/total)*100),
      nowpercent  = int(( now/total)*100);
  if (nowpercent>lastpercent)
    cout << "\r " << setw(3) << nowpercent << "% done\r" << flush;
  }

void announce (const string &name)
  {
  cout << endl << "+-";
  for (unsigned int m=0; m<name.length(); ++m) cout << "-";
  cout << "-+" << endl;
  cout << "| " << name << " |" << endl;
  cout << "+-";
  for (unsigned int m=0; m<name.length(); ++m) cout << "-";
  cout << "-+" << endl;
  }

void module_startup (const std::string &name, int argc, const char **,
  int argc_expected, const std::string &argv_expected)
  {
  announce (name);
  if (argc==argc_expected) return;
  cerr << "Usage: " << name << " " << argv_expected << endl;
  throw Message_error();
  }

void parse_file (const string &filename, map<string,string> &dict)
  {
  string line;
  int lineno=0;
  dict.clear();
  ifstream inp(filename.c_str());
  planck_assert (inp,"Could not open parameter file "+filename);
  while (inp)
    {
    getline(inp, line);
    ++lineno;
    line=line.substr(0,line.find_first_of("#"));
    line=trim(line);
    if (line.size()>0)
      {
      string::size_type eqpos=line.find("=");
      if (eqpos!=string::npos)
        {
        string key=trim(line.substr(0,eqpos)),
               value=trim(line.substr(eqpos+1,string::npos));
        if (key=="")
          cerr << "Warning: empty key in " << filename << ", line "
               << lineno << endl;
        else
          {
          if (dict.find(key)!=dict.end())
            cerr << "Warning: key " << key << " multiply defined in "
                 << filename << ", line " << lineno << endl;
          dict[key]=value;
          }
        }
      else
        cerr << "Warning: unrecognized format in " << filename << ", line "
             << lineno << ":\n" << line << endl;
      }
    }
  }
