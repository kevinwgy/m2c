%{
#include <math.h>
#include <stdio.h>
#include <string.h>
extern int  yyCmdflex(void);
extern void yyCmdferror(const char  *);
#include "Assigner.h"
#include "Dictionary.h"
#include "ResizeArray.h"

ResizeArray<Assigner *> underVar(0);
int underLevel = 0;

%}

%union
{
 int        ival;
 double     dval;
 int	    token;
 char *sval;
 Assigner *asgn;
}

%token UNDER EoF IntConstant DblConstant Symbol String END BOGUS

%left '+' '-'
%left '*' '/'
%left '^'

%type <ival>	IntConstant IntExpr
%type <dval>    DblConstant DblExpr
%type <token>   Symbol
%type <asgn>	Assignable
%type <sval>	String

%%

File: ValidInput EoF
          { return 0; }
	| ValidInput END
	{ return 0; }

ValidInput: 
	| ValidInput Assignment
	| ValidInput GroupInput

Assignment: Assignable '=' Symbol ';'
	{ $1->assignToken($3); }
	| Assignable '=' IntExpr ';'
	{ $1->assignInt($3); }
	| Assignable '=' DblExpr ';'
	{ $1->assignDouble($3); }
	| Assignable '=' String ';'
	{ $1->assignString(strdup($3)); }

GroupInput: UNDER Assignable '{'
	{
          underVar[underLevel] = $2; underLevel++;
       }
	 ValidInput '}'
	{ underLevel--; }


Assignable:
	Symbol
	{ 
           $$ = 0;
           for(int i = underLevel; i-- ; ) {
             $$ = underVar[i]->findSubToken($1);
             if($$ != 0)  
                break; 
           } 
           if($$ == 0)
             $$ = findSysObj($1); 
        }
	| Assignable '.' Symbol
	{ 
          $$ = $1->findSubToken($3); 
          if($$ == 0) {
	     fprintf(stderr, "ERROR: Structure element not found: %s\n",
                  dictionary->word($3).c_str());
             exit(-1);
          }
        }
        | Assignable '[' IntExpr ']'
        { $$ = $1->findIndexObject($3); 
          if($$ == 0) {
            fprintf(stderr, "ERROR: Object is not an array\n");
            exit(-1);
          }
        }
        


IntExpr:
	IntConstant
	| '(' IntExpr ')'
	{ $$ = $2; }
	| IntExpr '+' IntExpr
	{ $$ = $1 + $3; }
	| IntExpr '-' IntExpr
	{ $$ = $1 - $3; }
	| IntExpr '*' IntExpr
	{ $$ = $1 * $3; }
	| IntExpr '/' IntExpr
	{ $$ = $1 / $3; }

DblExpr:
	DblConstant
	| '(' DblExpr ')'
        { $$ = $2; }
        | DblExpr '+' DblExpr
        { $$ = $1 + $3; }
        | DblExpr '+' IntExpr
	{ $$ = $1 + $3; }
        | IntExpr '+' DblExpr
	{ $$ = $1 + $3; }
        | DblExpr '-' DblExpr
        { $$ = $1 - $3; }
        | DblExpr '-' IntExpr
        { $$ = $1 - $3; }
        | IntExpr '-' DblExpr
        { $$ = $1 - $3; }
        | DblExpr '*' DblExpr
        { $$ = $1 * $3; }
        | IntExpr '*' DblExpr
        { $$ = $1 * $3; }
        | DblExpr '*' IntExpr
        { $$ = $1 * $3; }
        | DblExpr '/' DblExpr
        { $$ = $1 / $3; }
        | IntExpr '/' DblExpr
        { $$ = $1 / $3; }
        | DblExpr '/' IntExpr
        { $$ = $1 / $3; }
%%
