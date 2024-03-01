#include <stdarg.h>
#include "Assigner.h"
#include "Dictionary.h"

template <class T>
ClassInt<T>::ClassInt(ClassAssigner *ca, const char *n, T *_ptr, int T::*_sp)
: Assigner(n)
{
  ptr = _ptr;
  sp = _sp;
  ca->addSmb(n, this);
}

template <class T>
ClassDouble<T>::ClassDouble(ClassAssigner *ca, const char *n, T *_ptr, double T::*_sp)
: Assigner(n)
{
  ptr = _ptr;
  sp = _sp;
  ca->addSmb(n, this);
}

template <class T>
ClassToken<T>::
ClassToken(ClassAssigner *ca, const char *n, T *_ptr, int T::*_sp, int nt, ...)
		: tk(nt), val(nt), Assigner(n)
{
  ptr = _ptr;
  token = _sp;
  ca->addSmb(n, this);
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

template <class T>
void
ClassToken<T>::assignToken(int t)
{
  for(int i = 0; i < (int)tk.size(); ++i)
    if(tk[i] == t) {
      ptr->*token = val[i];
      return;
    }
  fprintf(stderr, "ERROR: Token not understood: %s\n", 
		  dictionary->word(t).c_str());
  exit(-1);
}

template <class T>
ClassStr<T>::ClassStr(ClassAssigner *ca, const char *n, T *_ptr, const char *T::*_sp)
: Assigner(n)
{
 ptr = _ptr;
 str = _sp;
 ca->addSmb(n, this);
}

template <class T>
void
ClassStr<T>::assignString(const char *string)
{
 ptr->*str = string;
}

template <class Target>
SysMapObj<Target>::SysMapObj(const char *n, map<int, Target *> *_mapObj) : Assigner(n)  {

  mapObj = _mapObj;

}

template <class Target>
Assigner * SysMapObj<Target>::findIndexObject(int index)  {

  typename map<int, Target *>::iterator it = mapObj->find(index);
  if (it == mapObj->end())  {
    (*mapObj)[index] = new Target;
    it = mapObj->find(index);
  }
 
  return it->second->getAssigner();

}  
