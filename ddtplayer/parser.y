%{
#include <ddt_jit.hpp>
#include <hrtimer.h>

#include <mpi.h>
#include <stdio.h>
#include <assert.h>

#include <iostream>
#include <iomanip> 
#include <vector>
#include <queue>
#include <algorithm>

using namespace std;

void alloc_buffer(size_t size, void** buffer, int alignment);
void init_buffer(size_t size, void* buf, bool pattern);
int compare_buffers(size_t size, void *buf1, void *buf2);
void free_buffers(void *inbuf, void *outbuf);

extern FILE * yyin;
extern "C" int yylex (void);
void yyerror(const char *);

unsigned long long g_timerfreq;
HRT_TIMESTAMP_T start, stop;

struct Datatype {
	farc::Datatype *farc;
	MPI_Datatype    mpi;
};

vector<Datatype*> datatypes;

%}

%union {
	int val;
	struct Datatype *datatype;
};

%token <val> NUM
%token <sym> UNKNOWN SUBTYPE ELEM LP RP LB RB LC RC
%token <sym> CONTIGUOUS VECTOR HVECTOR HINDEXED STRUCT
%token <sym> BYTE_ CHAR_ INT_ DOUBLE_ FLOAT_

%type <datatype> datatype primitive derived contiguous vector hvector 

%start input

%%

input:
| input topdatatype 
;

topdatatype:
datatype {
	datatypes.push_back($1);

	Datatype *datatype = $1;


}
;

datatype:
primitive |
derived
;

primitive:
BYTE_ {
	$$ = new Datatype;
	$$->farc = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::BYTE);
	$$->mpi  = MPI_BYTE;
}
| CHAR_ {
	$$ = new Datatype;
	$$->farc = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::CHAR);
	$$->mpi  = MPI_CHAR;
}
| INT_ {
	$$ = new Datatype;
	$$->farc = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::INT);
	$$->mpi  = MPI_INT;
}
| FLOAT_ {
	$$ = new Datatype;
	$$->farc = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::FLOAT);
	$$->mpi  = MPI_FLOAT;
}
| DOUBLE_ {
	$$ = new Datatype;
	$$->farc = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::DOUBLE);
	$$->mpi  = MPI_DOUBLE;
}
;

derived:
contiguous
| vector
| hvector
/*| hindexed*/
;

contiguous:
CONTIGUOUS LP NUM RP LB datatype RB {
	$$ = $6;
	$$->farc = new farc::ContiguousDatatype($6->farc, $3);
	MPI_Type_contiguous($3, $6->mpi, &($$->mpi));
}
;

vector:
VECTOR LP NUM NUM NUM RP LB datatype RB { 
	$$ = $8;
	$$->farc = new farc::VectorDatatype($8->farc, $3, $4, $5);
	MPI_Type_vector($3, $4, $5, $8->mpi, &($$->mpi));
}
;

hvector:
HVECTOR LP NUM NUM NUM RP LB datatype RB {
	$$ = $8;
	$$->farc = new farc::HVectorDatatype($8->farc, $3, $4, $5);
	MPI_Type_hvector($3, $4, $5, $8->mpi, &($$->mpi));
}
;

/*
parenlist:
hindexed:
HINDEXED LP NUM COL parenlist RP LB datatype RB
;*/

%%

void yyerror(const char *s) {
	fprintf (stderr, "Error: %s\n", s);
}

void alloc_buffer(size_t size, void** buffer, int alignment) {
	if (alignment == 1) {
		*buffer = malloc(size);
	}
	else {
		posix_memalign(reinterpret_cast<void**>(buffer), 16, size);
	}
	assert(buffer != NULL);
}

void init_buffer(size_t size, void* buf, bool pattern) {
	if (pattern) {
		for (size_t i=0; i<size; i++) {
			((char*)buf)[i] = i+1;
		}
	}
	else {
		for (size_t i=0; i<size; i++) {
			((char*)buf)[i] = 0;
		}
	}
}

int compare_buffers(size_t size, void *mpi, void *farc) {
	int ret = 0;
    for (size_t i=0; i<size; i++) {
        if (((char*)farc)[i] != ((char*)farc)[i]) {
			printf("%3lu: mpi(%3i) != farc(%3i)\n", i, ((char*)mpi)[i], ((char*)farc)[i]);
			ret = -1;
		}
    }   
    return ret;
}

#define ALIGNMENT 1
#define WARMUP  5
#define NUMRUNS 10
void produce_report() {

	// Find text widths
	int name_w        = 0;
	int size_w        = 9;
	int mpi_commit_w  = 12;
	int farc_commit_w = 13;
	int mpi_pack_w    = 10;
	int farc_pack_w   = 11;
	int mpi_unpack_w  = 12;
	int farc_unpack_w = 13;

	for (unsigned int i=0; i<datatypes.size(); i++) {
		int textSize = datatypes[i]->farc->toString().size();
		if (textSize > name_w) {
			name_w = textSize;
		}
	}

	cout << setw(name_w)        << ""
		 << setw(size_w)        << "size"
		 << setw(mpi_commit_w)  << "mpi_commit"
		 << setw(farc_commit_w) << "farc_commit"
		 << setw(mpi_pack_w)    << "mpi_pack"
		 << setw(farc_pack_w)   << "farc_pack"
		 << setw(mpi_unpack_w)  << "mpi_unpack"
		 << setw(farc_unpack_w) << "farc_unpack"
		 << endl;
	
	// Produce report
	for (unsigned int i=0; i<datatypes.size(); i++) {
		Datatype *datatype = datatypes[i];

		double mpi_commit_time   = 0.0;
		double farc_commit_time  = 0.0;
		double mpi_pack_time     = 0.0;
		double farc_pack_time    = 0.0;
		double mpi_unpack_time   = 0.0;
		double farc_unpack_time  = 0.0;

		int size = datatype->farc->getSize();
		int extent = datatype->farc->getExtent();

		void *mpi_bigbuf, *mpi_smallbuf;
		alloc_buffer(size, &mpi_smallbuf, ALIGNMENT);
		alloc_buffer(extent, &mpi_bigbuf, ALIGNMENT);

		void *farc_bigbuf, *farc_smallbuf;
		alloc_buffer(size, &farc_smallbuf, ALIGNMENT);
		alloc_buffer(extent, &farc_bigbuf, ALIGNMENT);

		std::vector<uint64_t> times (NUMRUNS, 0);


		// mpi_commit
		for (unsigned int i=0; i<WARMUP; i++) {
			MPI_Type_commit(&(datatype->mpi));
		}
		for (unsigned int i=0; i<NUMRUNS; i++) {
			HRT_GET_TIMESTAMP(start);
			MPI_Type_commit(&(datatype->mpi));
			HRT_GET_TIMESTAMP(stop);
			HRT_GET_ELAPSED_TICKS(start, stop, &mpi_commit_time);
			HRT_GET_ELAPSED_TICKS(start, stop, &times[i]);
		}
		mpi_commit_time = HRT_GET_USEC(times[NUMRUNS/2]);

		// farc_commit
		for (unsigned int i=0; i<WARMUP; i++) {
			DDT_Commit(datatype->farc);
		}
		for (unsigned int i=0; i<NUMRUNS; i++) {
			HRT_GET_TIMESTAMP(start);
			DDT_Commit(datatype->farc);
			HRT_GET_TIMESTAMP(stop);
			HRT_GET_ELAPSED_TICKS(start, stop, &times[i]);
		}
		farc_commit_time = HRT_GET_USEC(times[NUMRUNS/2]);


		// mpi_pack
		init_buffer(extent, mpi_bigbuf, true);
		init_buffer(size, mpi_smallbuf, false);
		for (unsigned int i=0; i<WARMUP; i++) {
			int pos = 0;
			MPI_Pack(mpi_bigbuf, 1, datatype->mpi, mpi_smallbuf, size, &pos, MPI_COMM_WORLD);
		}
		for (unsigned int i=0; i<NUMRUNS; i++) {
			HRT_GET_TIMESTAMP(start);
			int pos = 0;
			MPI_Pack(mpi_bigbuf, 1, datatype->mpi, mpi_smallbuf, size, &pos, MPI_COMM_WORLD);
			HRT_GET_TIMESTAMP(stop);
			HRT_GET_ELAPSED_TICKS(start, stop, &times[i]);
		}
		mpi_pack_time = HRT_GET_USEC(times[NUMRUNS/2]);

		// farc pack
		init_buffer(extent, farc_bigbuf, true);
		init_buffer(size, farc_smallbuf, false);
		for (unsigned int i=0; i<WARMUP; i++) {
			DDT_Pack(farc_bigbuf, farc_smallbuf, datatype->farc, 1);
	   	}
		for (unsigned int i=0; i<NUMRUNS; i++) {
			HRT_GET_TIMESTAMP(start);
			DDT_Pack(farc_bigbuf, farc_smallbuf, datatype->farc, 1);
			HRT_GET_TIMESTAMP(stop);
			HRT_GET_ELAPSED_TICKS(start, stop, &times[i]);
		}
		farc_pack_time = HRT_GET_USEC(times[NUMRUNS/2]);

		// verify
		if (compare_buffers(extent, mpi_bigbuf, farc_bigbuf) != 0) {
			fprintf(stderr, "Error: %s: MPI and FARC input buffers differ after packing\n", 
				datatype->farc->toString().c_str());
		}
		if (compare_buffers(size, mpi_smallbuf, farc_smallbuf) != 0) {
			fprintf(stderr, "Error: %s: MPI and FARC output buffers differ after packing\n", 
				datatype->farc->toString().c_str());
		}


		// mpi_unpack
		init_buffer(size, mpi_smallbuf, true);
		init_buffer(extent, mpi_bigbuf, false);
		for (unsigned int i=0; i<WARMUP; i++) {
			int pos = 0;
			MPI_Unpack(mpi_smallbuf, size, &pos, mpi_bigbuf, 1, datatype->mpi, MPI_COMM_WORLD);
		}
		for (unsigned int i=0; i<NUMRUNS; i++) {
			HRT_GET_TIMESTAMP(start);
			int pos = 0;
			MPI_Unpack(mpi_smallbuf, size, &pos, mpi_bigbuf, 1, datatype->mpi, MPI_COMM_WORLD);
			HRT_GET_TIMESTAMP(stop);
			HRT_GET_ELAPSED_TICKS(start, stop, &times[i]);
		}
		mpi_unpack_time = HRT_GET_USEC(times[NUMRUNS/2]);

		// farc unpack
		init_buffer(size, farc_smallbuf, true);
		init_buffer(extent, farc_bigbuf, false);
		for (unsigned int i=0; i<WARMUP; i++) {
			DDT_Unpack(farc_smallbuf, farc_bigbuf, datatype->farc, 1);
		}
		for (unsigned int i=0; i<NUMRUNS; i++) {
			HRT_GET_TIMESTAMP(start);
			DDT_Unpack(farc_smallbuf, farc_bigbuf, datatype->farc, 1);
			HRT_GET_TIMESTAMP(stop);
			HRT_GET_ELAPSED_TICKS(start, stop, &times[i]);
		}
		farc_unpack_time = HRT_GET_USEC(times[NUMRUNS/2]);

		// verify
		if (compare_buffers(size, mpi_smallbuf, farc_smallbuf) != 0) {
			fprintf(stderr, "Error: %s: MPI and FARC input buffers differ after unpacking\n", 
				datatype->farc->toString().c_str());
		}
		if (compare_buffers(extent, mpi_bigbuf, farc_bigbuf) != 0) {
			fprintf(stderr, "Error: %s: MPI and FARC output buffers differ after unpacking\n", 
				datatype->farc->toString().c_str());
		}


		// free buffers
		free(mpi_bigbuf);
		free(mpi_smallbuf);
		free(farc_bigbuf);
		free(farc_smallbuf);


		// output
		cout << setw(name_w)        << datatype->farc->toString().c_str()
			 << setw(size_w)        << size
			 << setw(mpi_commit_w)  << setiosflags(ios::fixed) << setprecision(2) << mpi_commit_time
			 << setw(farc_commit_w) << setiosflags(ios::fixed) << setprecision(2) << farc_commit_time
			 << setw(mpi_pack_w)    << setiosflags(ios::fixed) << setprecision(3) << mpi_pack_time
			 << setw(farc_pack_w)   << setiosflags(ios::fixed) << setprecision(3) << farc_pack_time
			 << setw(mpi_unpack_w)  << setiosflags(ios::fixed) << setprecision(3) << mpi_unpack_time
			 << setw(farc_unpack_w) << setiosflags(ios::fixed) << setprecision(3) << farc_unpack_time
			 << endl;

	}
}

int main(int argc, char **argv) {
	int token;

	if (argc < 2) {
		fprintf(stderr, "%s <filename>\n", argv[0]);
		exit(1);
	}

	MPI_Init(&argc, &argv);
	farc::DDT_Init();
	HRT_INIT(1, g_timerfreq);

	yyin = fopen(argv[1], "r");
	if (yyin == NULL) {
		fprintf(stderr, "Error: could not open file %s\n", argv[1]);
		exit(1);
	}

	if (yyparse() == 0) {
		produce_report();
		printf("\n");
	}

	fclose(yyin);
}

