#include <mpi.h>
#include <cstdio>
#include <iomanip> 
#include <cstdlib>
#include <assert.h>
#include <iostream>
#include <hrtimer.h>
#include <algorithm>
#include <ddt_jit.hpp>

#include "cmdline.h"
#include "ddtplayer.hpp"

using namespace std;

gengetopt_args_info args_info;
unsigned long long g_timerfreq;
vector<struct Datatype> datatypes;

#define ALIGNMENT 1
#define WARMUP  5
#define NUMRUNS 10
#define CACHE_SIZE 16e6

#define CLEAR_CACHE                                    \
do {                                                   \
    uint8_t* A = (uint8_t*) malloc(CACHE_SIZE);        \
    for (size_t i=0; i<CACHE_SIZE; i++) {              \
        A[i] = rand() % 23;                            \
        A[i] += A[i / 2];                              \
    }                                                  \
    free(A);                                           \
} while(0)


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

#define TIME_COLD(code, median)                        \
do {                                                   \
	HRT_TIMESTAMP_T start, stop;                       \
	std::vector<uint64_t> times (NUMRUNS, 0);          \
	for (unsigned int i=0; i<WARMUP; i++) {            \
        CLEAR_CACHE;                                   \
		code;                                          \
	}                                                  \
	for (unsigned int i=0; i<NUMRUNS; i++) {           \
        CLEAR_CACHE;                                   \
		HRT_GET_TIMESTAMP(start);                      \
		code;                                          \
		HRT_GET_TIMESTAMP(stop);                       \
		HRT_GET_ELAPSED_TICKS(start, stop, &times[i]); \
	}                                                  \
	std::sort(times.begin(), times.end());             \
	median = HRT_GET_USEC(times[NUMRUNS/2]);           \
} while(0)


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
        if (((char*)farc)[i] != ((char*)mpi)[i]) {
			/* printf("%3lu: mpi(%3i) != farc(%3i)\n", i, ((char*)mpi)[i], ((char*)farc)[i]); */
			ret = -1;
		}
    }   
    return ret;
}

void produce_report() {

    int name_w = 0;

	for (unsigned int i=0; i<datatypes.size(); i++) {
		int textSize = datatypes[i].farc->toString(true).size();
		if (textSize > name_w) {
			name_w = textSize;
		}
	}

	cout << setw(0)         << "";
    cout << setw(9)         << "size";
    if (args_info.time_create_given) {
        cout << setw(26)   << "mpi_commit";
        cout << setw(26)  << "farc_commit";
    }
    if (args_info.time_hot_given) {
        cout << setw(26) << "mpi_pack_time_hot";
        cout << setw(26) << "farc_pack_time_hot";
        cout << setw(26) << "mpi_unpack_time_hot";
        cout << setw(26) << "farc_unpack_time_hot";
        cout << setw(26) << "pack_time_spdup_hot";
        cout << setw(26) << "unpack_time_spdup_hot";
    }
    if (args_info.time_cold_given) {
        cout << setw(26) << "mpi_pack_time_cold";
        cout << setw(26) << "farc_pack_time_cold";
        cout << setw(26) << "mpi_unpack_time_cold";
        cout << setw(26) << "farc_unpack_time_cold";
        cout << setw(26) << "pack_time_spdup_cold";
        cout << setw(26) << "unpack_time_spdup_cold";
    }

    cout << endl;
	
	// Produce report
	for (unsigned int i=0; i<datatypes.size(); i++) {
		HRT_TIMESTAMP_T start, stop;

		Datatype datatype = datatypes[i];

		double mpi_commit_time      = 0.0;
		double farc_commit_time     = 0.0;

        double mpi_pack_time_hot    = 0.0;
		double farc_pack_time_hot   = 0.0;
		double mpi_unpack_time_hot  = 0.0;
		double farc_unpack_time_hot = 0.0;

        double mpi_pack_time_cold    = 0.0;
		double farc_pack_time_cold   = 0.0;
		double mpi_unpack_time_cold  = 0.0;
		double farc_unpack_time_cold = 0.0;

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
        if (datatype.farc->getExtent() != extent_mpi) {
            cerr << "EXTENT MISSMATCH: MPI extent: " << extent_mpi
                 << " FARC extent: " << extent
                 << setw(name_w)
                 << datatype.farc->toString().c_str()
                 << endl;
            exit(EXIT_FAILURE);
        }

        // TODO Remove this. The *4 is used to provide "centered"
        // buffer in case of negative strides and in case of lower
        // bounds higher than 0.  It is a workaround until farc
        // supports these concepts.
        extent = extent * 4;
       
		void *mpi_bigbuf, *mpi_smallbuf;
		alloc_buffer(size, &mpi_smallbuf, ALIGNMENT);
		alloc_buffer(extent, &mpi_bigbuf, ALIGNMENT);

		void *farc_bigbuf, *farc_smallbuf;
		alloc_buffer(size, &farc_smallbuf, ALIGNMENT);
		alloc_buffer(extent, &farc_bigbuf, ALIGNMENT);

        // TODO Remove this.
        void *mpi_bigbuf_centered  = ((char*) mpi_bigbuf)  + (extent / 2);
        void *farc_bigbuf_centered = ((char*) farc_bigbuf) + (extent / 2);

		// mpi_commit
		TIME_HOT( MPI_Type_commit(&(datatype.mpi)), mpi_commit_time );

		// farc_commit
		TIME_HOT( DDT_Commit(datatype.farc), farc_commit_time );

		// mpi_pack
        init_buffer(extent, mpi_bigbuf, true);
		init_buffer(size, mpi_smallbuf, false);
		TIME_HOT( {int pos=0; MPI_Pack(mpi_bigbuf_centered, 1, datatype.mpi,
                    mpi_smallbuf, size, &pos, MPI_COMM_WORLD);}, 
                   mpi_pack_time_hot );
		TIME_COLD( {int pos=0; MPI_Pack(mpi_bigbuf_centered, 1, datatype.mpi,
                    mpi_smallbuf, size, &pos, MPI_COMM_WORLD);}, 
                   mpi_pack_time_cold );

		// farc pack
		init_buffer(extent, farc_bigbuf, true);
		init_buffer(size, farc_smallbuf, false);
		TIME_HOT( DDT_Pack(farc_bigbuf_centered, farc_smallbuf, datatype.farc, 1), farc_pack_time_hot);
		TIME_COLD( DDT_Pack(farc_bigbuf_centered, farc_smallbuf, datatype.farc, 1), farc_pack_time_cold);

		// verify
		if (compare_buffers(extent, mpi_bigbuf, farc_bigbuf) != 0) {
			cerr <<  "Error: " << datatype.farc->toString().c_str() 
                 << ": MPI and FARC input buffers differ after packing" << endl; 
				
		}
		if (compare_buffers(size, mpi_smallbuf, farc_smallbuf) != 0) {
			cerr <<  "Error: " << datatype.farc->toString().c_str() 
                 << ": MPI and FARC output buffers differ after packing" << endl; 
		}


		// mpi_unpack
		init_buffer(size, mpi_smallbuf, true);
		init_buffer(extent, mpi_bigbuf, false);
		TIME_HOT( {int pos=0; MPI_Unpack(mpi_smallbuf, size, &pos,
                   mpi_bigbuf_centered, 1, datatype.mpi, MPI_COMM_WORLD);},
                  mpi_unpack_time_hot);
		TIME_COLD( {int pos=0; MPI_Unpack(mpi_smallbuf, size, &pos,
                    mpi_bigbuf_centered, 1, datatype.mpi, MPI_COMM_WORLD);}, 
                   mpi_unpack_time_cold);

		// farc unpack
		init_buffer(size, farc_smallbuf, true);
		init_buffer(extent, farc_bigbuf, false);
		TIME_HOT(DDT_Unpack(farc_smallbuf, farc_bigbuf_centered, datatype.farc, 1), farc_unpack_time_hot);
		TIME_COLD(DDT_Unpack(farc_smallbuf, farc_bigbuf_centered, datatype.farc, 1), farc_unpack_time_cold);

		// verify
		if (compare_buffers(size, mpi_smallbuf, farc_smallbuf) != 0) {
			cerr <<  "Error: " << datatype.farc->toString().c_str() 
                 << ": MPI and FARC input buffers differ after unpacking" << endl; 
		}
        //TODO the *2 is used to provide "centered" buffer in case of negative strides
		if (compare_buffers(extent, mpi_bigbuf, farc_bigbuf) != 0) {
			cerr <<  "Error: " << datatype.farc->toString().c_str() 
                 << ": MPI and FARC output buffers differ after unpacking" << endl; 
		}

		double pack_speedup_hot   = ((mpi_pack_time_hot / farc_pack_time_hot)-1);
		double unpack_speedup_hot = ((mpi_unpack_time_hot / farc_unpack_time_hot)-1);
		double pack_speedup_cold   = ((mpi_pack_time_cold / farc_pack_time_cold)-1);
		double unpack_speedup_cold = ((mpi_unpack_time_cold / farc_unpack_time_cold)-1);


		// output
		cout.flags(std::ios::left);
		cout << setw(name_w)        << datatype.farc->toString(true).c_str();

		cout.flags(ios::right);
		cout.flags(ios::fixed);
		cout << setw(9)  << size;
        if (args_info.time_create_given) {
	        cout << setw(26) << setprecision(2) << mpi_commit_time;
		    cout << setw(26) << setprecision(2) << farc_commit_time;
        }
        if (args_info.time_hot_given) {
		    cout << setw(26) << setprecision(3) << mpi_pack_time_hot;
		    cout << setw(26) << setprecision(3) << farc_pack_time_hot;
		    cout << setw(26) << setprecision(3) << mpi_unpack_time_hot;
		    cout << setw(26) << setprecision(3) << farc_unpack_time_hot;
		    cout << setw(26) << setprecision(1) << pack_speedup_hot;
		    cout << setw(26) << setprecision(1) << unpack_speedup_hot;
        }
        if (args_info.time_hot_given) {
		    cout << setw(26) << setprecision(3) << mpi_pack_time_cold;
		    cout << setw(26) << setprecision(3) << farc_pack_time_cold;
		    cout << setw(26) << setprecision(3) << mpi_unpack_time_cold;
		    cout << setw(26) << setprecision(3) << farc_unpack_time_cold;
		    cout << setw(26) << setprecision(1) << pack_speedup_cold;
		    cout << setw(26) << setprecision(1) << unpack_speedup_cold;
        } 

		cout << endl;

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

    if (cmdline_parser(argc, argv, &args_info) != 0) {
        fprintf(stderr, "Could not parse commandline arguments\n");
        exit(1);
    }

    srand(time(NULL));

	MPI_Init(&argc, &argv);
	farc::DDT_Init();
	HRT_INIT(0, g_timerfreq);

	yyin = fopen(args_info.inputfile_arg, "r");
	if (yyin == NULL) {
		fprintf(stderr, "Error: could not open file %s\n", argv[1]);
		exit(1);
	}

	if (yyparse() == 0) {
		produce_report();
		printf("\n");
	}

    cmdline_parser_free(&args_info);
	fclose(yyin);
}

