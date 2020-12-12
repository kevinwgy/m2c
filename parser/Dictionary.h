#ifndef _DICTIONARY
#define _DICTIONARY

#ifdef OLD_STL
#include <map.h>
#include <string.h>
#else
#include <map>
#include <string>
using std::map;
using std::string;
#endif

class Assigner;

typedef map<string, int> DictMap;
typedef map<int, string> TokenMap;
typedef map<int, Assigner *> SysSmbMap;

class Dictionary {
   int nToken;
   DictMap map;
   TokenMap tkm;
 public:
   Dictionary() { nToken = 0; }
   int token(const char *text) { 
      DictMap::iterator it = map.find(text);
      if(it != map.end() ) return it->second;
      map[text] = nToken;
      tkm[nToken] = text;
      return nToken++;
   }
   string &word(int tk) { return tkm[tk]; }
};

extern Dictionary *dictionary;
extern SysSmbMap *sysSmb;

int addSysSymbol(const char *, Assigner *);
int findSysToken(const char *);
Assigner * findSysObj(int tk);

#endif
