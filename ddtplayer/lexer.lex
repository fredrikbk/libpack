%{
#include "parser.hpp"
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

[-]                  { return ELEM; }

{digit}+|-{digit}+   { yylval.val = atoi(yytext); return NUM; }


[\[\]\(\):,]         { return yytext[0]; }

[ \t\n\r]            /* skip whitespace */
.                    { printf("Unknown character [%c]\n", yytext[0]); return UNKNOWN; }
%%

int yywrap(void){return 1;}



/* int main(int argc, char **argv) { */
/* 	int token; */

/* 	if (argc < 2) { */
/* 		fprintf(stderr, "%s <filename>\n", argv[0]); */
/* 		exit(1); */
/* 	} */

/* 	/\* MPI_Init(&argc, &argv); *\/ */
/* 	/\* farc::DDT_Init(); *\/ */
/* 	// HRT_INIT(1, g_timerfreq); */

/* 	yyin = fopen(argv[1], "r"); */
/* 	if (yyin == NULL) { */
/* 		fprintf(stderr, "Could not open file %s\n", argv[1]); */
/* 	} */

/* 	while ((token = yylex()) != 0) {} */

/* 	fclose(yyin); */
/* } */
