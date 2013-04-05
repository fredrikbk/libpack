%{
#include <stdio.h>
int yylex (void);
void yyerror(const char *);
%}

%union {
  int val; 
};

%token <val> NUM
%token <sym> UNKNOWN SUBTYPE ELEM LB RB
%token <sym> BYTE CHAR INT DOUBLE FLOAT
%token <sym> CONTIGUOUS VECTOR HVECTOR HINDEXED STRUCT

%start input

%%

input:
| input datatype
;

datatype:
primitive |
derived
;

primitive:
  BYTE
| CHAR
| INT
| DOUBLE { printf("double\n") }
;

derived:
  contiguous
| vector
| hvector
;

contiguous:
 CONTIGUOUS NUM LB datatype RB      {printf("contiguous(%d)\n", $2); }
;

vector:
  VECTOR NUM NUM NUM LB datatype RB  { printf("vector(%d, %d, %d)\n", $2, $3, $4); }
;

hvector:
  HVECTOR NUM NUM NUM LB datatype RB { printf("hvector(%d, %d, %d)\n", $2, $3, $4); }
;


%%

void yyerror(const char *s) {
	fprintf (stderr, "Error: %s\n", s);
}
