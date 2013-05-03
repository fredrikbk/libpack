#include <list>
#include <mpi.h>
#include <ddt_jit.hpp>

extern unsigned long long g_timerfreq;
extern std::vector<struct Datatype> datatypes;

extern "C" {
extern FILE *yyin;
int yyparse();
int yylex (void);
void yyerror(const char *);
}

struct Datatype {
	farc::Datatype *farc;
	MPI_Datatype    mpi;
};

struct Datatypes {
    std::list<struct Datatype> types;
};


