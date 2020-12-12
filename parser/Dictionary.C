#include <stdio.h>
#include "Dictionary.h"

SysSmbMap *sysSmb = 0;
Dictionary *dictionary = 0;

int addSysSymbol(const char *name, Assigner *a)
{
 if(sysSmb == 0)
    sysSmb = new SysSmbMap;
 (*sysSmb)[findSysToken(name) ] = a;
 return 1;
}

int findSysToken(const char *str)
{
  if(dictionary == 0)
    dictionary = new Dictionary;
  return dictionary->token(str);
}

Assigner *
findSysObj(int tk)
{
  SysSmbMap::iterator it;
  it = sysSmb->find(tk);
  if(it == sysSmb->end()) {
    fprintf(stderr, "Error: Symbol not found: %s\n", dictionary->word(tk).c_str());
    return 0;
  }
  return it->second;
}
