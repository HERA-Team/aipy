/************************************************************************/
/*									*/
/*  Subroutines to aid interfacing between C and FORTRAN.		*/
/*									*/
/*  History:								*/
/*    rjs  ??????? Original version.					*/
/*    rjs  23dec92 Broke out into separate file.			*/
/*    rjs  05dec95 Comment out a dirty trick in zterm that was screwing */
/*		   up with some compilers!! 				*/
/*    pjt  17jun02 MIR4 prototypes                                      */
/************************************************************************/

#include <string.h>

void pad(char *string,int length)
/*
  This takes a zero-terminated string, and pads it with blanks up a certain
  length.

  Input:
    length	Length to pad to.
  Input/Output:
    string	Output is the blank padded version of the input.
------------------------------------------------------------------------*/
{
  int len0,i;
  char *s;

  len0 = strlen(string);
  s = string + len0;
  for(i=len0; i < length; i++) *s++ = ' ';
}
/************************************************************************/
char *zterm(char *string,int length)
/*
    This returns a pointer to a nul terminated string. This is usually
    called to convert strings from a FORTRAN to C manner. Its algorithm
    consists of firstly checking if the string is already null-terminated
    (without trailing blanks). If so return the string address. Otherwise
    trim back trailing blanks and copy the string (null terminating) to
    a buffer. The buffer is a circular one, with no checks for "overflows".
    Overflows are unlikely because the FORTRAN to C boundary is only
    spanned at most once at any given moment.

  Input:
    string	Points to the string of interest.
    length	The FORTRAN length of the string (may include blank padding).

  Output:
    zterm	Pointer to null terminated string.
------------------------------------------------------------------------*/
{
#define CIRBUFSIZE 2048
  static char buffer[CIRBUFSIZE];
  static int offset=0;

  char *s;

/* Trim back over blanks, and check if its already null terminated. */
/* If its already null terminated, there is nothing to do. */

  s = string + length;
  while(*--s == ' ' && length)length--;
/*  if(*(string+length) == 0)return(string); */

/* We have to put it in our circular buffer. Determine where to put it. */

  if(offset + length + 1 > CIRBUFSIZE) offset = 0;
  s = buffer + offset;
  memcpy(s,string,length);
  *(s+length) = 0;
  offset += length + 1;
  return(s);
}
