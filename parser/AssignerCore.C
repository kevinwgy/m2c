#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "Assigner.h"
#include "Dictionary.h"
void
Assigner::assignInt(int)
{
  fprintf(stderr, " *** ERROR: %s cannot be assigned an integer\n", name.c_str());
  exit(-1);
}

void
Assigner::assignDouble(double)
{
  fprintf(stderr, " *** ERROR: %s cannot be assigned a double\n", name.c_str());
  exit(-1);
}

void
Assigner::assignToken(int)
{
  fprintf(stderr, " *** ERROR: %s cannot be assigned a token\n", name.c_str());
  exit(-1);
}

Assigner *
Assigner::findSubToken(int)
{
  fprintf(stderr, " *** ERROR: %s does not have this subtoken\n", name.c_str());
  exit(-1);
  return 0;
}

void
Assigner::assignString(const char *)
{
  fprintf(stderr, " *** ERROR: %s cannot be assigned a string\n", name.c_str());
  exit(-1);
}

Assigner * Assigner::findIndexObject(int)  {

  fprintf(stderr, " *** ERROR: %s cannot be assigned a string\n", name.c_str());
  exit(-1);

}

SysIntObj::SysIntObj(const char *n, int *p)
: Assigner(n)
{
  val = p;
  addSysSymbol(n, this);
}


SysDoubleObj::SysDoubleObj(const char *n, double *p)
: Assigner(n)
{
  val = p;
  addSysSymbol(n, this);
}

SysStrObj::SysStrObj(const char *n, const char **p)
: Assigner(n)
{
 val = p;
 addSysSymbol(n, this);
}

SysTokenObj::SysTokenObj(const char *n, int *p, int nt, ...) : tk(nt), val(nt),
 Assigner(n)
{
  ptr = p;
  addSysSymbol(n, this);
  va_list vl;
  va_start(vl, nt);
  
  for(int i = 0; i < nt; ++i) {
    const char *tokenString = va_arg(vl, const char *);
    int v = va_arg(vl, int);
    tk[i] = findSysToken(tokenString);
    val[i] = v;
    
  }
  va_end(vl);
}

void
SysTokenObj::assignToken(int t)
{
  for(int i = 0; i < tk.size(); ++i)
    if(tk[i] == t) {
      *ptr = val[i];
      return;
    }
  fprintf(stderr, "ERROR: Token not understood: %s\n", 
		  dictionary->word(t).c_str());
  exit(1);
}
/*
ClassAssigner::ClassAssigner(const char *n, int ns)
		: subAsgn(ns), subToken(ns)
{
  addSysSymbol(n, this);
  cn = 0;
}
*/
ClassAssigner::ClassAssigner(const char *n, int ns, ClassAssigner *p)
		: subAsgn(ns), subToken(ns), Assigner(n)
{
  cn = 0;
  if (p) p->addSmb(n, this);
  else addSysSymbol(n, this);
}

Assigner *
ClassAssigner::findSubToken(int t)
{
  for(int i = 0; i < cn; ++i)
    if(subToken[i] == t)
      return subAsgn[i];
/*
  fprintf(stderr, "ERROR: Structure element not found: %s\n", 
		  dictionary->word(t).c_str());
*/
  return 0;
}

void
ClassAssigner::addSmb(const char *smb, Assigner *asgn)
{
  subToken[cn] = findSysToken(smb);
  subAsgn[cn] = asgn;
  ++cn;
}
