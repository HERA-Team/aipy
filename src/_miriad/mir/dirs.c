/************************************************************************/
/*									*/
/*	Various directory related routines - VMS only.			*/
/*									*/
/*  History:								*/
/*    rjs  Dark-ages Original version.					*/
/*    rjs  15nov89   Added the dexpand_c routine.			*/
/*    rjs  27apr90   Added the ddelete_c routine.			*/
/*    pjt  16dec90   Added comment that this is VMS only...             */
/*    rjs  25jun91   Improved the dexpand_c routine substantially.	*/
/*    rjs  26aug93   Add drmdir.					*/
/************************************************************************/

#define private static
#define MAXPATH 128
#include <descrip.h>

#define NULL 0
#define STRING struct dsc$descriptor_s
#define STATIC_STRING(name,string) struct dsc$descriptor_s \
	name = { sizeof(string)-1, DSC$K_DTYPE_T, DSC$K_CLASS_S, string}
#define DYNAMIC_STRING(name,string) \
	(name).dsc$w_length = strlen(string);\
	(name).dsc$b_dtype  = DSC$K_DTYPE_T;\
	(name).dsc$b_class  = DSC$K_CLASS_S;\
	(name).dsc$a_pointer = string

struct context { STRING dir;
		 char *contxt,name[MAXPATH];};

static STATIC_STRING(asterisk,"*.*;0");
static STATIC_STRING(dirspec,"*.DIR;0");
static STATIC_STRING(delete_default,".;*");

#define Strcpy(a,b) (void)strcpy(a,b)
#define Malloc(x) malloc((unsigned)x)

char *strcpy(),*malloc();
void dopendir_c(),dreaddir_c(),dclosedir_c();
private void dheadtail();

/************************************************************************/
void drmdir_c(path,iostat)
char *path;
int *iostat;
/*
  Remove this directory.

  Input:
    path	A directory name. It will look something like
		device:[alpha.beta]. This routine converts it to
		device:[alpha]beta, and then invokes rmdir.
------------------------------------------------------------------------*/
{
  char name[MAXPATH],*s;

  strcpy(name,"rmdir ");
  strcat(name,path);
  s = name + strlen(name) - 1;
  *s = 0;
  while(*s != '.')s--;
  *s = ']';
  *iostat = system(name);
}
/************************************************************************/
void ddelete_c(path,iostat)
char *path;
int *iostat;
/*
  This deletes a file.
------------------------------------------------------------------------*/
{
  char name[MAXPATH];
  STRING Path;

  DYNAMIC_STRING(Path,path);
  *iostat = lib$delete_file(&Path,&delete_default);
  if(*iostat == 1) *iostat = 0;
}
/************************************************************************/
int dexpand_c(template,output,length)
char *template,*output;
int length;
/*
  This expands wildcards, matching them with files.

  Input:
    template	The input character string, containing the wildcards.
    length	The length of the output buffer.
  Output:
    output	All the files matching "template". Filenames are separated
		by commas.
------------------------------------------------------------------------*/
{
  int l,iostat;
  char *contxt,name[MAXPATH],head[MAXPATH],tail[MAXPATH],*p,*s;
  char temp[MAXPATH],head2[MAXPATH],tail2[MAXPATH];
  STRING Template;
  STATIC_STRING(Name,name);

/* Break the name into a head and tail and translate the head. */

  dheadtail(template,temp,tail);
  if(*tail == 0) return(-1);
  iostat = 0;
  if(*temp != 0)dtrans_c(temp,head,&iostat);
  else		*head = 0;
  if(iostat)return(-1);

/* Concatenate the head and tail, and form a descriptor. */

  strcpy(temp,head); strcat(temp,tail);
  DYNAMIC_STRING(Template,temp);

/* Loop the loop. */

  contxt = NULL;
  s = output;
  *s = 0;
  while(lib$find_file(&Template,&Name,&contxt,&dirspec) % 2){

/* Strip trailing blanks. */

    p = name + Name.dsc$w_length;
    while(*(p-1) == ' ')p--;
    *p = 0;

/* Form the head and tail again, and append the old head to the new tail. */

    dheadtail(name,head2,tail2);
    if(*tail2 == 0){lib$find_file_end(&contxt); return(-1);}
    strcpy(name,head); strcat(name,tail2);
    l = strlen(name);

/* Make sure we have enough space, and then copy it. */

    if(length < l){lib$find_file_end(&contxt); return(-1);}
    Strcpy(s,name);
    s += l;
    length -= l;
    *s++ = ',';
    length--;
  }
  if(s != output) *--s = 0;
  lib$find_file_end(&contxt);
  return(s-output);
}
/************************************************************************/
private void dheadtail(in,head,tail)
char *in,*head,*tail;
/*
  This performs converts the file spec to lower case, strips off a
  trailing .dir;1 specification, and then breaks the name into a head
  and tail.
------------------------------------------------------------------------*/
{
#define U2L ('a' - 'A')
  char *p,*q;
  int length;

/* Copy the input to the head array, converting to lower case as we
   go. */

  p = head;
  q = in;
  while(*q && *q != ';'){
    if(*q >= 'A' && *q <= 'Z') *p++ = *q + U2L;
    else		       *p++ = *q;
    q++;
  }
  length = q - in;
  *p = 0;

/* Discard a trailing .dir portion. */

  if(length > 4) if(!strcmp(".dir",head+length-4)){
    *(head+length-4) = 0;
    length -= 4;
  }

/* Trim back from the end of the string, looking for the first "/", "."
   or "]". */

  p = head + length;
  while(p > head && *(p-1) != '.' && *(p-1) != '/' && *(p-1) != ']')p--;
  strcpy(tail,p);
  if(p != head) p--;
  *p = 0;
}  
/************************************************************************/
void dopendir_c(contxt,directory)
char **contxt,*directory;
/*
  Prepare to search a directory.
------------------------------------------------------------------------*/
{
  struct context *s;
  s = (struct context *)Malloc(sizeof(struct context));
  Strcpy(s->name,directory);
  DYNAMIC_STRING(s->dir,s->name);
  s->contxt = NULL;
  *contxt = (char *)s;
}
/************************************************************************/
void dreaddir_c(contxt,line,length)
char *contxt,*line;
int length;
{
  struct context *s;
  char path[MAXPATH],*p;
  STATIC_STRING(Path,path);

  s = (struct context *)contxt;
  if(lib$find_file(&s->dir,&Path,&s->contxt,&asterisk)%2){
    p = path + Path.dsc$w_length;

/* Trim back across blanks and ther version. */

    while(*--p != ';');

/* Zero terminate it. Directories are made to end in /. */

    if(*(p-1) == '.') *(p-1) = 0;
    else if(!strncmp(p-4,".DIR",4)){
      *(p-4) = '/';
      *(p-3) = 0;
    }else *p = 0;

/* Skip past the device/directory spec, and copy to the output,
   converting to lower case. */

    p = path;
    while(*p++ != ']');
    while(*p){
      if(*p >= 'A' && *p <= 'Z')*line++ = *p++ - 'A' + 'a';
      else			*line++ = *p++;
    }
  }
  *line = 0;
}
/************************************************************************/
void dclosedir_c(contxt)
char *contxt;
/*
  Delete memory associated with this directopry search.
------------------------------------------------------------------------*/
{
  struct context *s;
  s = (struct context *)contxt;
  lib$find_file_end(&s->contxt);
  free(contxt);
}
