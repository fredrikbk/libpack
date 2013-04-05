%{
#include "parser.h"
#define YY_NO_UNPUT

%}

digit         [0-9]
letter        [a-zA-Z]

%%
"byte"               { return BYTE; }
"char"               { return CHAR; }
"int"                { return INT; }
"double"             { return DOUBLE; }
"float"              { return FLOAT; }

"contiguous"         { return CONTIGUOUS; }
"vector"             { return VECTOR; }
"hvector"            { return HVECTOR; }
"hindexed"           { return HINDEXED; }

[-]                  { return ELEM; }

{digit}+|-{digit}+   { yylval.val = atoi(yytext); return NUM; }

[\[]                 { return LB; }
[\]]                 { return RB; }

[ \t\n\r]                  /* skip whitespace */
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
