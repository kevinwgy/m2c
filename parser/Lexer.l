%{
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Dictionary.h"
#include "Parser.tab.h"
int yyCmdfline_num = 1;
void C_comments(void);
void CPP_comments(void);
%}
DIGIT  [0-9]
SPACES [ \t]+

%%

"\n"                 { yyCmdfline_num++; }

"\/\/" while(yyinput() != '\n') ; yyCmdfline_num++;

"\/\*"        {
                     for(;;)
                     {
			    int c;
                             while((c = yyinput())!='*')
                                     if(c == '\n') yyCmdfline_num++;
                             if(yyinput()=='/')
                                     break;
                             else

                                     unput(yytext[yyleng-1]);
                     }
             }


"\""[A-z0-9._\/\-]*"\""	{
     char *string = new char[ strlen(yyCmdftext)-1]; 
     string[0] = 0;
     strncat(string, yyCmdftext+1, strlen(yyCmdftext)-2);
     yyCmdflval.sval = string; return String; }

"-"|"+"|"*"|"["|"]"|"^"|"{"|"}"|"="|","|"("|")"|":"|";"|"."|"/"   {
    return *yyCmdftext; }

{SPACES}                                 /* Just discard spaces */

"end"	return END;

"under" return UNDER;

[_A-Za-z][A-Za-z0-9_]*  {
     yyCmdflval.token = dictionary->token(yytext); return Symbol; 
    }

[+-]?{DIGIT}+ {
    yyCmdflval.ival = atoi(yyCmdftext) ; return IntConstant ; 
  }

[+-]?{DIGIT}+"."{DIGIT}*                  |
[+-]?"."{DIGIT}+                          |
[+-]?{DIGIT}+"."{DIGIT}*[eE][+-]?{DIGIT}+ |
[+-]?{DIGIT}+[eE][+-]?{DIGIT}+            |
[+-]?"."{DIGIT}+[eE][+-]?{DIGIT}+        {
                     yyCmdflval.dval = atof(yyCmdftext) ; return DblConstant ;
		      }

[+-]?{DIGIT}+"."{DIGIT}*[dD][+-]?{DIGIT}+ |
[+-]?{DIGIT}+[dD][+-]?{DIGIT}+            |
[+-]?"."{DIGIT}+[dD][+-]?{DIGIT}+  {
	int il;
	for(il = 0; yyCmdftext[il] != 'd' && yyCmdftext[il] != 'D';) ++il;
	yyCmdftext[il] = 'e'; yyCmdflval.dval = atof(yyCmdftext) ; return DblConstant; 
	}

<<EOF>> { return EoF; }
%%

void
yyCmdferror(const char  *)
{
 fprintf(stderr,"Syntax Error near line %d : %s\n", yyCmdfline_num,yyCmdftext);
}

int yyCmdfwrap()
{
 // End of file is terminal
 return 1;
}

