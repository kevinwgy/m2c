#ifndef _ASSIGNER_H_
#define _ASSIGNER_H_

#include <stdio.h>

#ifdef OLD_STL
#include <vector.h>
#include <string.h>
#include <map.h>
#else
#include <vector>
#include <string>
#include <map>
using std::string;
using std::vector;
using std::map;
#endif

class Assigner {
 public:
   string name;
   Assigner(const char *n) { name = n; }
   virtual void assignInt(int);
   virtual void assignDouble(double);
   virtual void assignToken(int);
   virtual void assignString(const char *);
   virtual Assigner *findSubToken(int);

   virtual Assigner *findIndexObject(int);
};

class SysIntObj : public Assigner {
   int *val;
  public:
   SysIntObj(const char *n, int *p);
   void assignInt(int v) { *val = v; }
   void assignDouble(double v) { *val = int(v); }
};

class SysDoubleObj : public Assigner {
   double *val;
  public:
   SysDoubleObj(const char *n, double *p);
   void assignInt(int v) { *val = v; }
   void assignDouble(double v) { *val = v; }
};

class SysTokenObj : public Assigner {
   int *ptr;
   vector<int> tk;
   vector<int> val;
 public:
   SysTokenObj(const char *n, int *ptr, int nt, ...);
   void assignToken(int);
};

class SysStrObj : public Assigner {
  const char **val;
public:
  SysStrObj(const char *n, const char **p);
  virtual void assignString(const char *p) { *val = p; }
};

template<class Target>
class SysMapObj : public Assigner  {

    map<int, Target *> *mapObj;

  public:
    SysMapObj(const char *n, map<int, Target *> *); 
    Assigner *findIndexObject(int);
};

// template <class T>
class ClassAssigner : public Assigner {
    // T *ptr;
    vector<Assigner *> subAsgn;
    vector<int> subToken;
    int cn;

    vector<map<int, Assigner *> > subAsgnMap;
    vector<int> subTokenMap;
    int cnMap;
  public:
  //ClassAssigner(char *n, int ns);             
    ClassAssigner(const char *n, int ns, ClassAssigner * = 0);
    virtual void addSmb(const char *, Assigner *);
    Assigner *findSubToken(int);
};

class RootClassAssigner: public ClassAssigner  {

  public:
    RootClassAssigner() : ClassAssigner("", 0)  {};
    void addSmb(const char *, Assigner *) {};
};


template <class T>
class ClassInt : public Assigner {
    T *ptr;
    int T::*sp;
  public:
    ClassInt(ClassAssigner *, const char *n, T *ptr, int T::*sp);
    void assignInt(int v) { ptr->*sp = v; }
};

template <class T>
class ClassDouble : public Assigner {
    T *ptr;
    double T::*sp;
  public:
    ClassDouble(ClassAssigner *, const char *n, T *ptr, double T::*sp);
    void assignInt(int v) { ptr->*sp = v; }
    void assignDouble(double v) { ptr->*sp = v; }
};


template <class T>
class ClassToken : public Assigner {
    T *ptr;
    int T::*token;
    vector<int> tk;
    vector<int> val;
  public:
    ClassToken(ClassAssigner *, const char *n, T *ptr, int T::*sp, int nt, ...);
    void assignToken(int);
};

template <class T>
class ClassStr  : public Assigner {
  T *ptr;
  const char *T::*str;
public:
    ClassStr(ClassAssigner *, const char *n, T *ptr, const char *T::*sp);
    void assignString(const char *str);
};

#include "Assigner.C"

#endif
