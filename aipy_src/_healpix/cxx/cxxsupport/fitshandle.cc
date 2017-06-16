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
 *  This file contains the implementation of the FITS I/O helper class
 *  used by the Planck LevelS package.
 *
 *  Copyright (C) 2002, 2003, 2004 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include "fitsio.h"
#include "fitshandle.h"
#include "cxxutils.h"

using namespace std;

string type2char (int type)
  {
  if (type==TFLOAT) return "E";
  if (type==TDOUBLE) return "D";
  if (type==TINT) return "J";
  if (type==TINT32BIT) return "J";
  if (type==TSTRING) return "A";
  throw Message_error("wrong datatype in type2char()");
  }

string type2asciiform (int type)
  {
  if (type==TFLOAT) return "E14.7";
  if (type==TDOUBLE) return "E23.15";
  if (type==TINT) return "I11";
  if (type==TINT32BIT) return "I11";
  throw Message_error("wrong datatype in type2asciiform()");
  }

void fitshandle::check_errors()
  {
  if (status==0) return;
  char msg[81];
  fits_get_errstatus (status, msg);
  cerr << msg << endl;
  while (fits_read_errmsg(msg)) cerr << msg << endl;
  throw Message_error("FITS error");
  }

void fitshandle::clean_data()
  {
  if (!fptr) return;
  axes_.clear();
  columns_.clear();
  hdutype_=INVALID;
  bitpix_=INVALID;
  nrows_=0;
  }

void fitshandle::clean_all()
  {
  if (!fptr) return;
  clean_data();
  fits_close_file (fptr, &status);
  check_errors();
  fptr=0;
  }

void fitshandle::init_image()
  {
  int naxis;
  fits_get_img_type (fptr, &bitpix_, &status);
  fits_get_img_dim (fptr, &naxis, &status);
  check_errors();
  arr<long> naxes(naxis);
  fits_get_img_size (fptr, naxis, naxes, &status);
  for (long m=0; m<naxis; ++m) axes_.push_back(naxes[naxis-m-1]);
  check_errors();
  }

void fitshandle::init_asciitab()
  {
  char ttype[81], tunit[81], tform[81];
  int ncol, typecode;
  fits_get_num_cols (fptr, &ncol, &status);
  fits_get_num_rows (fptr, &nrows_, &status);
  check_errors();
  for (int m=1; m<=ncol; ++m)
    {
    fits_get_acolparms (fptr, m, ttype, 0, tunit, tform,
                        0, 0, 0, 0, &status);
    fits_ascii_tform (tform, &typecode, 0,0, &status);
    check_errors();
    columns_.push_back (fitscolumn (ttype,tunit,1,typecode));
    }
  }

void fitshandle::init_bintab()
  {
  char ttype[81], tunit[81], tform[81];
  long repc;
  int ncol, typecode;
  fits_get_num_cols (fptr, &ncol, &status);
  fits_get_num_rows (fptr, &nrows_, &status);
  check_errors();
  for (int m=1; m<=ncol; ++m)
    {
    fits_get_bcolparms (fptr, m, ttype, tunit, tform, &repc,
                        0, 0, 0, 0, &status);
    fits_binary_tform (tform, &typecode, 0,0, &status);
    check_errors();
    columns_.push_back (fitscolumn (ttype,tunit,repc,typecode));
    }
  }

void fitshandle::init_data()
  {
  switch(hdutype_)
    {
    case IMAGE_HDU:
      init_image();
      break;
    case ASCII_TBL:
      init_asciitab();
      break;
    case BINARY_TBL:
      init_bintab();
      break;
    default:
      throw Message_error("init_data(): wrong HDU type");
      break;
    }
  }

void fitshandle::read_col (int colnum, void *data, long ndata, int dtype,
  long offset)
  {
  assert_table_hdu("fitshandle::read_column()",colnum);
  long repc = columns_[colnum-1].repcount();
  planck_assert (ndata<=(repc*nrows_-offset),
                 "read_column(): array too large");
  long frow = offset/repc+1;
  long felem = offset%repc+1;
  fits_read_col (fptr, dtype, colnum, frow, felem, ndata, 0, data, 0, &status);
  check_errors();
  }
void fitshandle::write_col (int colnum, const void *data, long ndata,
  int dtype, long offset)
  {
  assert_table_hdu("fitshandle::write_column()",colnum);
  long repc = columns_[colnum-1].repcount();
  long frow = offset/repc+1;
  long felem = offset%repc+1;
  fits_write_col (fptr, dtype, colnum, frow, felem, ndata,
    const_cast<void *>(data), &status);
  fits_get_num_rows (fptr, &nrows_, &status);
  check_errors();
  }

void fitshandle::open (const string &fname, int rwmode)
  {
  clean_all();
  fits_open_file(&fptr, fname.c_str(), rwmode, &status);
  check_errors();
  }

void fitshandle::create (const string &fname)
  {
  clean_all();
  fits_create_file(&fptr, fname.c_str(), &status);
  fits_write_imghdr(fptr, 8, 0, 0, &status); // insert empty PHDU
  fits_write_date(fptr, &status);
  check_errors();
  }

void fitshandle::goto_hdu (int hdu)
  {
  clean_data();
  fits_movabs_hdu(fptr, hdu, &hdutype_, &status);
  check_errors();
  init_data();
  }

void fitshandle::insert_bintab (const vector<fitscolumn> &cols)
  {
  clean_data();
  int ncol=cols.size();
  arr2b<char> ttype(ncol,81), tform(ncol,81), tunit(ncol,81);

  for (long m=0; m<ncol; ++m)
    {
    strcpy (ttype[m], cols[m].name().c_str());
    strcpy (tunit[m], cols[m].unit().c_str());
    ostringstream x;
    x << cols[m].repcount() << type2char (cols[m].type());
    strcpy (tform[m], x.str().c_str());
    }
  fits_insert_btbl (fptr, nrows_, ncol, ttype, tform, tunit,
    const_cast<char *>("xtension"), 0, &status);
  check_errors();
  hdutype_=BINARY_TBL;
  columns_=cols;
  }

void fitshandle::insert_asctab (const vector<fitscolumn> &cols)
  {
  clean_data();
  int ncol=cols.size();
  arr2b<char> ttype(ncol,81), tform(ncol,81), tunit(ncol,81);

  for (long m=0; m<ncol; ++m)
    {
    strcpy (ttype[m], cols[m].name().c_str());
    strcpy (tunit[m], cols[m].unit().c_str());
    ostringstream x;
    if (cols[m].type()!=TSTRING)
      {
      planck_assert (cols[m].repcount()==1,"bad repcount for ASCII table");
      x << type2asciiform (cols[m].type());
      }
    else
      {
      x << "A" << dataToString(cols[m].repcount());
      }
    strcpy (tform[m], x.str().c_str());
    }
  fits_insert_atbl (fptr, 0, nrows_, ncol, ttype, 0, tform, tunit,
    const_cast<char *>("xtension"), &status);
  check_errors();
  hdutype_=ASCII_TBL;
  columns_=cols;
  }

void fitshandle::insert_image (int btpx, const vector<long> &Axes)
  {
  clean_data();
  hdutype_=IMAGE_HDU;
  bitpix_=btpx;
  axes_=Axes;
  arr<long> tmpax(Axes.size());
  for (long m=0; m<long(Axes.size()); m++) tmpax[m]=Axes[Axes.size()-1-m];
  fits_insert_img (fptr, btpx, Axes.size(), tmpax, &status);
  check_errors();
  }

template<typename T>
  void fitshandle::insert_image (int btpx, const arr2<T> &data)
  {
  clean_data();
  hdutype_=IMAGE_HDU;
  bitpix_=btpx;
  axes_.push_back(data.size1());
  axes_.push_back(data.size2());
  arr<long> tmpax(axes_.size());
  for (long m=0; m<long(axes_.size()); m++) tmpax[m]=axes_[axes_.size()-1-m];
  fits_insert_img (fptr, btpx, axes_.size(), tmpax, &status);
  arr2<T> &tmparr = const_cast<arr2<T> &> (data);
  fits_write_img (fptr, FITSUTIL<T>::DTYPE, 1, axes_[0]*axes_[1],
    &tmparr[0][0], &status);
  check_errors();
  }

template void fitshandle::insert_image (int btpx, const arr2<double> &data);
template void fitshandle::insert_image (int btpx, const arr2<float> &data);
template void fitshandle::insert_image (int btpx, const arr2<int> &data);

void fitshandle::copy_historified_header (const fitshandle &orig)
  {
  const char *inclist[] = { "*" };
  const char *exclist[] = {
      "SIMPLE","BITPIX","NAXIS","NAXIS#","PCOUNT","GCOUNT",
      "EXTEND","ORIGIN","DATE*","TFIELDS","TTYPE#","TFORM#",
      "TUNIT#","EXTNAME","CTYPE#","CRVAL#","CRPIX#","CDELT#",
      "XTENSION","INSTRUME","TELESCOP","PDMTYPE","TBCOL#" };
  char card[81];
  string card2;
  orig.assert_connected("fitshandle::copy_historified_header()");
  assert_connected("fitshandle::copy_historified_header()");
  fits_read_record (orig.fptr, 0, card, &status);
  check_errors();
  while (true)
    {
    fits_find_nextkey (orig.fptr, const_cast<char **>(inclist), 1,
      const_cast<char **>(exclist), 23, card, &status);
    if (status!=0) break;
    card2=trim(card);
    if (card2!="END" && card2!="COMMENT" && card2!="HISTORY")
      {
      if (card2.find("COMMENT")==0)
        card2.replace(0,7,"HISTORY");
      if (card2.find("HISTORY")!=0)
        card2.insert(0,"HISTORY ");
      if (card2.length()<=80)
        fits_write_record (fptr, card2.c_str(), &status);
      else
        {
        fits_write_record (fptr, card2.substr(0,80).c_str(), &status);
        card2=card2.substr(80,string::npos);
        card2.insert(0,"HISTORY ");
        fits_write_record (fptr, card2.c_str(),&status);
        }
      }
    check_errors();
    }
  if (status==KEY_NO_EXIST) { fits_clear_errmsg(); status=0; }
  check_errors();
  }

void fitshandle::copy_header (const fitshandle &orig)
  {
  const char *inclist[] = { "*" };
  const char *exclist[] = {
      "SIMPLE","BITPIX","NAXIS","NAXIS#","PCOUNT","GCOUNT",
      "EXTEND","ORIGIN","DATE*","TFIELDS","TTYPE#","TFORM#",
      "TUNIT#","EXTNAME","CTYPE#","CRVAL#","CRPIX#","CDELT#",
      "XTENSION","INSTRUME","TELESCOP","PDMTYPE","TBCOL#" };
  char card[81];
  string card2;
  orig.assert_connected("fitshandle::copy_header()");
  assert_connected("fitshandle::copy_header()");
  fits_read_record (orig.fptr, 0, card, &status);
  check_errors();
  while (true)
    {
    fits_find_nextkey (orig.fptr, const_cast<char **>(inclist), 1,
      const_cast<char **>(exclist), 23, card, &status);
    if (status!=0) break;
    card2=trim(card);
    if (card2!="END" && card2!="COMMENT" && card2!="HISTORY")
      fits_write_record (fptr, card, &status);
    check_errors();
    }
  if (status==KEY_NO_EXIST) { fits_clear_errmsg(); status=0; }
  check_errors();
  }

void fitshandle::check_key_present(const string &name)
  {
  char card[81];
  fits_read_card(fptr, const_cast<char *>(name.c_str()), card, &status);
  if (status==KEY_NO_EXIST)
    { fits_clear_errmsg(); status=0; return; }
  check_errors();
// FIXME: disabled for now; but this issue needs to be resolved!
//  cerr << "Warning: key " << name << " set more than once!" << endl;
  }

template<typename T> void fitshandle::add_key (const string &name,
  const T &value, const string &comment)
  {
  assert_connected("fitshandle::add_key()");
  check_key_present (name);
  fits_write_key (fptr, FITSUTIL<T>::DTYPE, const_cast<char *>(name.c_str()),
    const_cast<T *>(&value), const_cast<char *>(comment.c_str()), &status);
  check_errors();
  }

template void fitshandle::add_key(const string &name,
  const int &value, const string &comment);
template void fitshandle::add_key(const string &name,
  const long &value, const string &comment);
template void fitshandle::add_key(const string &name,
  const float &value, const string &comment);
template void fitshandle::add_key(const string &name,
  const double &value, const string &comment);
template<> void fitshandle::add_key(const string &name,
  const bool &value, const string &comment)
  {
  assert_connected("fitshandle::add_key()");
  check_key_present (name);
  int val=value;
  fits_write_key (fptr, TLOGICAL, const_cast<char *>(name.c_str()),
    &val, const_cast<char *>(comment.c_str()),
    &status);
  check_errors();
  }
template<> void fitshandle::add_key (const string &name,
  const string &value, const string &comment)
  {
  assert_connected("fitshandle::add_key()");
  check_key_present (name);
  fits_write_key (fptr, TSTRING, const_cast<char *>(name.c_str()),
    const_cast<char *>(value.c_str()), const_cast<char *>(comment.c_str()),
    &status);
  check_errors();
  }

template<typename T> void fitshandle::update_key (const string &name,
  const T &value, const string &comment)
  {
  assert_connected("fitshandle::update_key()");
  fits_update_key (fptr, FITSUTIL<T>::DTYPE, const_cast<char *>(name.c_str()),
    const_cast<T *>(&value), const_cast<char *>(comment.c_str()), &status);
  check_errors();
  }

template void fitshandle::update_key(const string &name,
  const int &value, const string &comment);
template void fitshandle::update_key(const string &name,
  const long &value, const string &comment);
template void fitshandle::update_key(const string &name,
  const float &value, const string &comment);
template void fitshandle::update_key(const string &name,
  const double &value, const string &comment);
template<> void fitshandle::update_key(const string &name,
  const bool &value, const string &comment)
  {
  assert_connected("fitshandle::update_key()");
  int val=value;
  fits_update_key (fptr, TLOGICAL, const_cast<char *>(name.c_str()),
    &val, const_cast<char *>(comment.c_str()),
    &status);
  check_errors();
  }
template<> void fitshandle::update_key (const string &name,
  const string &value, const string &comment)
  {
  assert_connected("fitshandle::update_key()");
  fits_update_key (fptr, TSTRING, const_cast<char *>(name.c_str()),
    const_cast<char *>(value.c_str()), const_cast<char *>(comment.c_str()),
    &status);
  check_errors();
  }

void fitshandle::delete_key (const string &name)
  {
  assert_connected("fitshandle::delete_key()");
  fits_delete_key (fptr, const_cast<char *>(name.c_str()), &status);
  check_errors();
  }

void fitshandle::add_comment(const string &comment)
  {
  assert_connected("fitshandle::add_comment()");
  fits_write_comment(fptr,const_cast<char *>(comment.c_str()),&status);
  check_errors();
  }

template<typename T> void fitshandle::get_key (const string &name, T &value)
  {
  assert_connected("fitshandle::get_key()");
  fits_read_key (fptr, FITSUTIL<T>::DTYPE, const_cast<char *>(name.c_str()),
    &value, 0, &status);
  if (status==KEY_NO_EXIST) throw Message_error
    ("Fitshandle::get_key(): key "+name+" not found");
  check_errors();
  }
template void fitshandle::get_key(const string &name,int &value);
template void fitshandle::get_key(const string &name,long &value);
template void fitshandle::get_key(const string &name,float &value);
template void fitshandle::get_key(const string &name,double &value);
template<> void fitshandle::get_key(const string &name,bool &value)
  {
  assert_connected("fitshandle::get_key()");
  int val;
  fits_read_key (fptr, TLOGICAL, const_cast<char *>(name.c_str()), &val, 0,
    &status);
  if (status==KEY_NO_EXIST) throw Message_error
    ("Fitshandle::get_key(): key "+name+" not found");
  check_errors();
  value=val;
  }
template<> void fitshandle::get_key (const string &name,string &value)
  {
  char tmp[2049];
  assert_connected("fitshandle::get_key()");
  fits_read_key (fptr, TSTRING, const_cast<char *>(name.c_str()), tmp, 0,
    &status);
  if (status==KEY_NO_EXIST) throw Message_error
    ("Fitshandle::get_key(): key "+name+" not found");
  check_errors();
  value=tmp;
  }

bool fitshandle::key_present(const string &name)
  {
  char card[81];
  fits_read_card(fptr, const_cast<char *>(name.c_str()), card, &status);
  if (status==KEY_NO_EXIST)
    { fits_clear_errmsg(); status=0; return false; }
  check_errors();
  return true;
  }

void fitshandle::assert_pdmtype (const string &pdmtype)
  {
  string type;
  get_key("PDMTYPE",type);
  if (pdmtype==type) return;
  cerr << "PDMTYPE " << pdmtype << " expected, but found " << type << endl;
  }

void fitshandle::read_column_raw_void
  (int colnum, void *data, int type, long num, long offset)
  {
  switch (type)
    {
    case PLANCK_INT32:
      read_col (colnum, data, num, TINT, offset); break;
    case PLANCK_FLOAT32:
      read_col (colnum, data, num, TFLOAT, offset); break;
    case PLANCK_FLOAT64:
      read_col (colnum, data, num, TDOUBLE, offset); break;
    case PLANCK_STRING:
      {
      string *data2 = static_cast<string *> (data);
      assert_table_hdu("fitshandle::read_column()",colnum);
      planck_assert (num<=(nrows_-offset),
        "read_column(): array too large");
      arr2b<char> tdata(num, columns_[colnum-1].repcount()+1);
      fits_read_col (fptr, TSTRING, colnum, offset+1, 1, num,
        0, tdata, 0, &status);
      check_errors();
      for (long m=0;m<num;++m) data2[m]=tdata[m];
      break;
      }
    default:
      throw Message_error ("unsupported data type in read_column_raw_void()");
    }
  }

void fitshandle::write_column_raw_void
  (int colnum, const void *data, int type, long num, long offset)
  {
  switch (type)
    {
    case PLANCK_INT32:
      write_col (colnum, data, num, TINT, offset); break;
    case PLANCK_FLOAT32:
      write_col (colnum, data, num, TFLOAT, offset); break;
    case PLANCK_FLOAT64:
      write_col (colnum, data, num, TDOUBLE, offset); break;
    case PLANCK_STRING:
      {
      const string *data2 = static_cast<const string *> (data);
      assert_table_hdu("fitshandle::write_column()",colnum);
      int stringlen = columns_[colnum-1].repcount()+1;
      arr2b<char> tdata(num, stringlen);
      for (long m=0;m<num;++m)
        {
        strncpy(tdata[m],data2[m].c_str(),stringlen-1);
        tdata[m][stringlen-1] = '\0';
        }
      fits_write_col (fptr, TSTRING, colnum, offset+1, 1, num,
        tdata, &status);
      fits_get_num_rows (fptr, &nrows_, &status);
      check_errors();
      break;
      }
    default:
      throw Message_error ("unsupported data type in write_column_raw_void()");
    }
  }

template<typename T> void fitshandle::write_image (const arr2<T> &data)
  {
  assert_image_hdu("fitshandle::write_image()");
  planck_assert (axes_.size()==2, "wrong number of dimensions");
  planck_assert (axes_[0]==data.size1(), "wrong size of dimension 1");
  planck_assert (axes_[1]==data.size2(), "wrong size of dimension 2");

  fits_write_img (fptr, FITSUTIL<T>::DTYPE, 1, axes_[0]*axes_[1],
    const_cast<T *>(&data[0][0]), &status);
  check_errors();
  }

template void fitshandle::write_image (const arr2<float> &data);
template void fitshandle::write_image (const arr2<double> &data);
template void fitshandle::write_image (const arr2<int> &data);

template<typename T> void fitshandle::write_subimage
  (const arr<T> &data, long offset)
  {
  assert_image_hdu("fitshandle::write_subimage()");
  fits_write_img (fptr, FITSUTIL<T>::DTYPE, 1+offset,
      data.size(), const_cast<T *>(&data[0]), &status);
  check_errors();
  }

template void fitshandle::write_subimage (const arr<float> &data, long offset);
template void fitshandle::write_subimage
  (const arr<double> &data, long offset);
template void fitshandle::write_subimage (const arr<int> &data, long offset);

template<typename T> void fitshandle::read_image (arr2<T> &data)
  {
  assert_image_hdu("fitshandle::read_image()");
  planck_assert (axes_.size()==2, "wrong number of dimensions");
  data.alloc(axes_[0], axes_[1]);
  fits_read_img (fptr, FITSUTIL<T>::DTYPE, 1, axes_[0]*axes_[1], 0,
    &data[0][0], 0, &status);
  check_errors();
  }

template void fitshandle::read_image (arr2<float> &data);
template void fitshandle::read_image (arr2<double> &data);
template void fitshandle::read_image (arr2<int> &data);

template<typename T> void fitshandle::read_subimage
  (arr2<T> &data, int xl, int yl)
  {
  assert_image_hdu("fitshandle::read_subimage()");
  planck_assert (axes_.size()==2, "wrong number of dimensions");
  for (int m=0; m<data.size1(); ++m)
    fits_read_img (fptr, FITSUTIL<T>::DTYPE, (xl+m)*axes_[1]+yl+1,
      data.size2(), 0, &data[m][0], 0, &status);
  check_errors();
  }

template void fitshandle::read_subimage (arr2<float> &data, int xl, int yl);
template void fitshandle::read_subimage (arr2<double> &data, int xl, int yl);
template void fitshandle::read_subimage (arr2<int> &data, int xl, int yl);

template<typename T> void fitshandle::read_subimage (arr<T> &data, long offset)
  {
  assert_image_hdu("fitshandle::read_subimage()");
  fits_read_img (fptr, FITSUTIL<T>::DTYPE, 1+offset,
      data.size(), 0, &data[0], 0, &status);
  check_errors();
  }

template void fitshandle::read_subimage (arr<float> &data, long offset);
template void fitshandle::read_subimage (arr<double> &data, long offset);
template void fitshandle::read_subimage (arr<int> &data, long offset);

void fitshandle::add_healpix_keys (int datasize)
  {
  int nside = isqrt(datasize/12);
  planck_assert (12*nside*nside==datasize, "Wrong Healpix map size");

  add_key ("PIXTYPE",string("HEALPIX"),"HEALPIX pixelisation");
  add_key ("ORDERING",string("RING"),
           "Pixel ordering scheme, either RING or NESTED");
  add_key ("NSIDE",nside,"Resolution parameter for HEALPIX");
  add_key ("FIRSTPIX",0,"First pixel # (0 based)");
  add_key ("LASTPIX",datasize-1,"Last pixel # (0 based)");
  add_key ("INDXSCHM",string("IMPLICIT"), "Indexing : IMPLICIT or EXPLICIT");
  add_key ("GRAIN",0,"Grain of pixel indexing");
  }
