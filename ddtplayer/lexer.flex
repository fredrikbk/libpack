%{
#include "parser.tab.hpp"
#define YY_NO_UNPUT

%}

digit         [0-9]
letter        [a-zA-Z]

%%
"byte"               { return BYTE_; }
"char"               { return CHAR_; }
"int"                { return INT_; }
"double"             { return DOUBLE_; }
"float"              { return FLOAT_; }

"ctg"                { return CONTIGUOUS; }
"vec"                { return VECTOR; }
"hvec"               { return HVECTOR; }
"hidx"               { return HINDEXED; }
"idxb"               { return INDEXEDBLOCK; }

[-]                  { return ELEM; }

{digit}+|-{digit}+   { yylval.val = atoi(yytext); return NUM; }


[\[\]\(\):,;]        { return yytext[0]; }

[ \t\n\r]            /* skip whitespace */

"//".*               { /*Ignore comment*/ }

.                    { printf("Unknown character [%c]\n", yytext[0]); return UNKNOWN; }
%%

int yywrap(void){return 1;}

