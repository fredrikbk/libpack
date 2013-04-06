%{
#include <ddt_jit.hpp>
#include <hrtimer.h>


#include <stdio.h>
#include <assert.h>
#include <queue>

using namespace std;
using namespace farc;

extern "C" int yylex (void);
extern "C" void yyerror(const char *);

extern unsigned long long g_timerfreq;
HRT_TIMESTAMP_T start, stop;
uint64_t duration;

queue<Datatype*> farcDDTs;

%}

%union {
  int val; 
};

%token <val> NUM
%token <sym> UNKNOWN SUBTYPE ELEM LB RB
%token <sym> BYTE_ CHAR_ INT_ DOUBLE_ FLOAT_
%token <sym> CONTIGUOUS VECTOR HVECTOR HINDEXED STRUCT

%start input

%%

input:
| input topdatatype
;

topdatatype:
datatype {
	Datatype *datatype = farcDDTs.front();
	farcDDTs.pop();
	assert(farcDDTs.empty());

	datatype->print("");

	
	HRT_GET_TIMESTAMP(start);
	DDT_Commit(datatype);
	HRT_GET_TIMESTAMP(stop);
	HRT_GET_ELAPSED_TICKS(start, stop, &duration);
	printf("Time: %10.3lf us\n", HRT_GET_USEC(duration)); 

	


}
;

datatype:
primitive |
derived
;

primitive:
BYTE_
| CHAR_
| INT_
| DOUBLE_ {
	printf("double\n");
	farcDDTs.push(new PrimitiveDatatype(PrimitiveDatatype::DOUBLE));
}
;

derived:
contiguous
| vector
| hvector
;

contiguous:
CONTIGUOUS NUM LB datatype RB {
	printf("contiguous(%d)\n", $2);
}
;

vector:
VECTOR NUM NUM NUM LB datatype RB { 
	printf("vector(%d, %d, %d)\n", $2, $3, $4);
	Datatype *base = farcDDTs.front();
	farcDDTs.pop();
	farcDDTs.push(new VectorDatatype(base, $2, $3, $4));
}
;

hvector:
HVECTOR NUM NUM NUM LB datatype RB { printf("hvector(%d, %d, %d)\n", $2, $3, $4); }
;


%%

void yyerror(const char *s) {
	fprintf (stderr, "Error: %s\n", s);
}
