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
#include <list>
#include <algorithm>

using namespace std;

int calc_num(int start, int stride, int stop);

extern FILE * yyin;
extern "C" int yylex (void);
void yyerror(const char *);

unsigned long long g_timerfreq;

struct Datatype {
	farc::Datatype *farc;
	MPI_Datatype    mpi;
};

struct Datatypes {
	list<struct Datatype> types;
};

struct Index {
	int displ;
	int blocklen;
};

struct Indices {
	vector<struct Index *> indices;
};

vector<struct Datatype> datatypes;

%}

%union {
	int val;
	struct Datatypes *types;

	struct Index *index;
	struct Indices *indices;

	struct {
		int start;
		int stop;
		int stride;
	} range;
};

%token <val> NUM
%token <sym> UNKNOWN SUBTYPE ELEM
%token <sym> CONTIGUOUS VECTOR HVECTOR HINDEXED STRUCT
%token <sym> BYTE_ CHAR_ INT_ DOUBLE_ FLOAT_

%type <types> datatype primitive derived contiguous vector hvector hindexed
%type <indices> idxentries
%type <index>   idxentry
%type <range>   range

%start input

%%

input:
| input topdatatype 
;

topdatatype:
datatype {
	datatypes.insert( datatypes.end(), $1->types.begin(), $1->types.end() );
	free($1);
}
;

datatype:
primitive |
derived
;

primitive:
BYTE_ {
	$$ = new Datatypes;
	struct Datatype datatype;
	datatype.farc = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::BYTE);
	datatype.mpi  = MPI_BYTE;
	$$->types.push_back(datatype);
}
| CHAR_ {
	$$ = new Datatypes;
	struct Datatype datatype;
	datatype.farc = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::CHAR);
	datatype.mpi  = MPI_CHAR;
	$$->types.push_back(datatype);
}
| INT_ {
	$$ = new Datatypes;
	struct Datatype datatype;
	datatype.farc = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::INT);
	datatype.mpi  = MPI_INT;
	$$->types.push_back(datatype);
}
| FLOAT_ {
	$$ = new Datatypes;
	struct Datatype datatype;
	datatype.farc = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::FLOAT);
	datatype.mpi  = MPI_FLOAT;
	$$->types.push_back(datatype);
}
| DOUBLE_ {
	$$ = new Datatypes;
	struct Datatype datatype;
	datatype.farc = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::DOUBLE);
	datatype.mpi  = MPI_DOUBLE;
	$$->types.push_back(datatype);
}
;

derived:
contiguous
| vector
| hvector
| hindexed
;

range:
NUM {
	$$.start = $1;
	$$.stride = 1;
	$$.stop = $1;
}
| NUM ':' NUM ':' NUM {
	$$.start = $1;
	$$.stride = $3;
	$$.stop = $5;
}



contiguous:
CONTIGUOUS '(' range ')' '[' datatype ']' {
	$$ = new Datatypes;

	list<struct Datatype> *subtypes = &($6->types);
	list<struct Datatype> *types = &($$->types);

	for(list<struct Datatype>::iterator subtype = subtypes->begin();
		subtype != subtypes->end(); subtype++) {
		for (int count=$3.start; count <= $3.stop; count += $3.stride) {
			Datatype type;
			type.farc = new farc::ContiguousDatatype(count, subtype->farc);
			MPI_Type_contiguous(count, subtype->mpi, &(type.mpi));
			types->push_back(type);
		}
	}

	free($6);
}
;

vector:
VECTOR '(' range range range ')' '[' datatype ']' { 
	$$ = new Datatypes;

	list<struct Datatype> *subtypes = &($8->types);
	list<struct Datatype> *types = &($$->types);

	for(list<struct Datatype>::iterator subtype = subtypes->begin();
		subtype != subtypes->end(); subtype++) {
		for (int count=$3.start; count <= $3.stop; count += $3.stride) {
			for (int blocklen=$4.start; blocklen <= $4.stop; blocklen += $4.stride) {
				for (int stride=$5.start; stride <= $5.stop; stride += $5.stride) {
					Datatype type;
					type.farc = new farc::VectorDatatype(count, blocklen, stride, subtype->farc);
					MPI_Type_vector(count, blocklen, stride, subtype->mpi, &(type.mpi));
					types->push_back(type);
				}
			}
		}
	}

	free($8);
}
;

hvector:
HVECTOR '(' range range range ')' '[' datatype ']' {
	$$ = new Datatypes;

	list<struct Datatype> *subtypes = &($8->types);
	list<struct Datatype> *types = &($$->types);

	for(list<struct Datatype>::iterator subtype = subtypes->begin();
		subtype != subtypes->end(); subtype++) {
		for (int count=$3.start; count <= $3.stop; count += $3.stride) {
			for (int blocklen=$4.start; blocklen <= $4.stop; blocklen += $4.stride) {
				for (int stride=$5.start; stride <= $5.stop; stride += $5.stride) {
					Datatype type;
					type.farc = new farc::HVectorDatatype(count, blocklen, stride, subtype->farc);
					MPI_Type_hvector(count, blocklen, stride, subtype->mpi, &(type.mpi));
					types->push_back(type);
				}
			}
		}
	}

	free($8);
}
;

idxentry:
NUM ',' NUM {
	$$ = new Index;
	$$->displ = $1;
	$$->blocklen = $3;
}
;

idxentries:
/* empty rule */ {
	$$ = new Indices;
}
| idxentries idxentry {
	$$ = $1;
	$$->indices.push_back($2);
}
;

hindexed:
HINDEXED '(' idxentries ')' '[' datatype ']' {
	$$ = new Datatypes;

	list<struct Datatype> *subtypes = &($6->types);
	list<struct Datatype> *types = &($$->types);

	unsigned int num = $3->indices.size();
	long *displs = (long*)malloc(num * sizeof(long));
	int *blocklens = (int*)malloc(num * sizeof(int));

	for(int i=0; i<num; i++) {
		displs[i] = $3->indices[i]->displ;
		blocklens[i] = $3->indices[i]->blocklen;
		free($3->indices[i]);
	}
	free($3);

	for(list<struct Datatype>::iterator subtype = subtypes->begin();
		subtype != subtypes->end(); subtype++) {
		Datatype type;
		type.farc = new farc::HIndexedDatatype(num, blocklens, displs, subtype->farc);
		MPI_Type_hindexed(num, blocklens, displs, subtype->mpi, &(type.mpi));
		types->push_back(type);
	}

	free(displs);
	free(blocklens);
	free($6);
}
;

%%

void yyerror(const char *s) {
	fprintf (stderr, "Error: %s\n", s);
}

int calc_num(int start, int stride, int stop) {
	return 1+(stop-start)/stride;
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

#define WARMUP  5
#define NUMRUNS 10
#define TIME_HOT(code, median)                         \
do {                                                   \
	HRT_TIMESTAMP_T start, stop;                       \
	std::vector<uint64_t> times (NUMRUNS, 0);          \
	for (unsigned int i=0; i<WARMUP; i++) {            \
		code;                                          \
	}                                                  \
	for (unsigned int i=0; i<NUMRUNS; i++) {           \
		HRT_GET_TIMESTAMP(start);                      \
		code;                                          \
		HRT_GET_TIMESTAMP(stop);                       \
		HRT_GET_ELAPSED_TICKS(start, stop, &times[i]); \
	}                                                  \
	std::sort(times.begin(), times.end());             \
	median = HRT_GET_USEC(times[NUMRUNS/2]);           \
} while(0)

#define ALIGNMENT 1
void produce_report() {
	// Find text widths
	int name_w         = 0;
	int size_w         = 9;
	int mpi_commit_w   = 12;
	int farc_commit_w  = 13;
	int mpi_pack_w     = 10;
	int farc_pack_w    = 11;
	int pack_spdup_w   = 12;
	int mpi_unpack_w   = 12;
	int farc_unpack_w  = 13;
	int unpack_spdup_w = 14;

	for (unsigned int i=0; i<datatypes.size(); i++) {
		int textSize = datatypes[i].farc->toString().size();
		if (textSize > name_w) {
			name_w = textSize;
		}
	}

	cout << setw(name_w)         << ""
		 << setw(size_w)         << "size"
		 << setw(mpi_commit_w)   << "mpi_commit"
		 << setw(farc_commit_w)  << "farc_commit"
		 << setw(mpi_pack_w)     << "mpi_pack"
		 << setw(farc_pack_w)    << "farc_pack"
		 << setw(mpi_unpack_w)   << "mpi_unpack"
		 << setw(farc_unpack_w)  << "farc_unpack"
		 << setw(pack_spdup_w)   << "pack_spdup"
		 << setw(unpack_spdup_w) << "unpack_spdup"
		 << endl;
	
	// Produce report
	for (unsigned int i=0; i<datatypes.size(); i++) {
		HRT_TIMESTAMP_T start, stop;

		Datatype datatype = datatypes[i];

		double mpi_commit_time   = 0.0;
		double farc_commit_time  = 0.0;
		double mpi_pack_time     = 0.0;
		double farc_pack_time    = 0.0;
		double mpi_unpack_time   = 0.0;
		double farc_unpack_time  = 0.0;

		int size = datatype.farc->getSize();
		int size_mpi;
        MPI_Type_size(datatype.mpi, &size_mpi);
        if (size != size_mpi) {
            printf("SIZE MISSMATCH: MPI size: %i FARC size: %i\n", size_mpi, size);
		    cout << setw(name_w)        << datatype.farc->toString().c_str() << std::endl;
            exit(EXIT_FAILURE);
        }
		int extent = datatype.farc->getExtent();
        MPI_Aint extent_mpi;
        MPI_Type_extent(datatype.mpi, &extent_mpi);
        if (extent != extent_mpi) {
            cerr << "EXTENT MISSMATCH: MPI extent: " << extent_mpi
                 << " FARC extent: " << extent
                 << setw(name_w)
                 << datatype.farc->toString().c_str()
                 << endl;
            exit(EXIT_FAILURE);
        }
       

		void *mpi_bigbuf, *mpi_smallbuf, *mpi_bigbuf_centered;
		alloc_buffer(size, &mpi_smallbuf, ALIGNMENT);
        //TODO the *2 is used to provide "centered" buffer in case of negative strides
		alloc_buffer(extent*2, &mpi_bigbuf, ALIGNMENT);
        mpi_bigbuf_centered = ((char*) mpi_bigbuf) + extent;

		void *farc_bigbuf, *farc_smallbuf, *farc_bigbuf_centered;
		alloc_buffer(size, &farc_smallbuf, ALIGNMENT);
        //TODO the *2 is used to provide "centered" buffer in case of negative strides
		alloc_buffer(extent*2, &farc_bigbuf, ALIGNMENT);
        farc_bigbuf_centered = ((char*) farc_bigbuf) + extent;


		// mpi_commit
		TIME_HOT( MPI_Type_commit(&(datatype.mpi)), mpi_commit_time );

		// farc_commit
		TIME_HOT( DDT_Commit(datatype.farc), farc_commit_time );

		// mpi_pack
        
		init_buffer(extent*2, mpi_bigbuf, true);
		init_buffer(size, mpi_smallbuf, false);
		TIME_HOT( {int pos=0; MPI_Pack(mpi_bigbuf, 1, datatype.mpi, mpi_smallbuf, size, &pos, MPI_COMM_WORLD);}, mpi_pack_time );

		// farc pack
        //TODO the *2 is used to provide "centered" buffer in case of negative strides
		init_buffer(extent*2, farc_bigbuf, true);
		init_buffer(size, farc_smallbuf, false);
		TIME_HOT( DDT_Pack(farc_bigbuf_centered, farc_smallbuf, datatype.farc, 1), farc_pack_time);

		// verify
        //TODO the *2 is used to provide "centered" buffer in case of negative strides
		if (compare_buffers(extent*2, mpi_bigbuf, farc_bigbuf) != 0) {
			fprintf(stderr, "Error: %s: MPI and FARC input buffers differ after packing\n", 
				datatype.farc->toString().c_str());
		}
		if (compare_buffers(size, mpi_smallbuf, farc_smallbuf) != 0) {
			fprintf(stderr, "Error: %s: MPI and FARC output buffers differ after packing\n", 
				datatype.farc->toString().c_str());
		}


		// mpi_unpack
		init_buffer(size, mpi_smallbuf, true);
		init_buffer(extent, mpi_bigbuf, false);
		TIME_HOT( {int pos=0; MPI_Unpack(mpi_smallbuf, size, &pos, mpi_bigbuf, 1, datatype.mpi, MPI_COMM_WORLD);}, mpi_unpack_time);

		// farc unpack
		init_buffer(size, farc_smallbuf, true);
		init_buffer(extent, farc_bigbuf, false);
		TIME_HOT(DDT_Unpack(farc_smallbuf, farc_bigbuf_centered, datatype.farc, 1), farc_unpack_time);

		// verify
		if (compare_buffers(size, mpi_smallbuf, farc_smallbuf) != 0) {
			fprintf(stderr, "Error: %s: MPI and FARC input buffers differ after unpacking\n", 
				datatype.farc->toString().c_str());
		}
        //TODO the *2 is used to provide "centered" buffer in case of negative strides
		if (compare_buffers(extent*2, mpi_bigbuf, farc_bigbuf) != 0) {
			fprintf(stderr, "Error: %s: MPI and FARC output buffers differ after unpacking\n", 
				datatype.farc->toString().c_str());
		}

		double pack_speedup   = ((mpi_pack_time/farc_pack_time)-1)*100;
		double unpack_speedup = ((mpi_unpack_time/farc_unpack_time)-1)*100;

		// output
		cout.flags(std::ios::left);
		cout << setw(name_w)        << datatype.farc->toString().c_str();

		cout.flags(ios::right);
		cout.flags(ios::fixed);
		cout << setw(size_w)           << size
			 << setw(mpi_commit_w)     << setprecision(2) << mpi_commit_time
			 << setw(farc_commit_w)    << setprecision(2) << farc_commit_time
			 << setw(mpi_pack_w)       << setprecision(3) << mpi_pack_time
			 << setw(farc_pack_w)      << setprecision(3) << farc_pack_time
			 << setw(mpi_unpack_w)     << setprecision(3) << mpi_unpack_time
			 << setw(farc_unpack_w)    << setprecision(3) << farc_unpack_time
			 << setw(pack_spdup_w-1)   << setprecision(1) << pack_speedup     << "%"
			 << setw(unpack_spdup_w-1) << setprecision(1) << unpack_speedup   << "%"
			 << endl;

		// free buffers
		free(mpi_bigbuf);
		free(mpi_smallbuf);
		free(farc_bigbuf);
		free(farc_smallbuf);
	}

	// free datatypes
	for (unsigned int i=0; i<datatypes.size(); i++) {
		Datatype datatype = datatypes[i];

		// TODO: Also delete children
		DDT_Free(datatype.farc);
		MPI_Type_free(&(datatype.mpi));
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
	HRT_INIT(0, g_timerfreq);

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
