#include <mpi.h>
#include <ddt_jit.hpp>
#include <hrtimer.h>

using namespace std;

unsigned long long g_timerfreq;

extern FILE * yyin;
extern "C" int yylex();

extern "C" int yyparse();

int main(int argc, char **argv) {
	int token;

	if (argc < 2) {
		fprintf(stderr, "%s <filename>\n", argv[0]);
		exit(1);
	}

	/* MPI_Init(&argc, &argv); */
	/* farc::DDT_Init(); */
	// HRT_INIT(1, g_timerfreq);

	yyin = fopen(argv[1], "r");

	yyparse();

	// if (yyin == NULL) {
		// fprintf(stderr, "Could not open file %s\n", argv[1]);
	// }
	// while ((token = yylex()) != 0) {}

	fclose(yyin);
}

