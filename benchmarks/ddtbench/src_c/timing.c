#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "../hrtimer/hrtimer.h"

#include "mpi.h"

//! handles all the time measurement

//! we need this for the hrtimer stuff
unsigned long long g_timerfreq;
double g_hrttimer_startvalue;


#if TEST_TYPE != 2
  static HRT_TIMESTAMP_T timings[1025];
#endif

static int ids[1025];
static int epoch_used[5];
static int counter;

static char testname[50], method[50];
static char epoch_name[5][20];
static int bytes;
static MPI_File filehandle_values;

static int max_tests, current_tests;

/* PAPI relevant variables */
#if TEST_TYPE > 1
  #include "papi.h"
  static int PapiNumberEvents; // max 2
  static int PapiEventSet;
  static int PapiEventCode[2];
#endif

static long long PapiEventCounts[2050];
static char PapiEventName1[20] = "NA\0";
static char PapiEventName2[20] = "NA\0";

//! sets some init parameter like the testname, method and number of bytes
//! also sets static the epoch names, the epoch_used array and the initial timing with a special id (zero)

void timing_init( char* ptestname, char* pmethod, int pbytes ) {

  int i;

  counter = 1;

  strncpy( testname, ptestname, 50 );
  strncpy( method, pmethod, 50 );
  bytes = pbytes;

  snprintf( &epoch_name[0][0], 20, "ddt_create_overhead");
  snprintf( &epoch_name[1][0], 20, "pack");
  snprintf( &epoch_name[2][0], 20, "communication");
  snprintf( &epoch_name[3][0], 20, "unpack");
  snprintf( &epoch_name[4][0], 20, "ddt_free_overhead");

  for( i = 0 ; i < 5 ; i++) {
    epoch_used[i] = 0;
  }

  ids[0] = 0;

#if TEST_TYPE != 2
  HRT_GET_TIMESTAMP(timings[0]);
#endif

#if TEST_TYPE > 1
  /* start counter */
  PAPI_start( PapiEventSet);
#endif
}


void timing_hrt_init() {
  HRT_INIT(1, g_timerfreq);
}

//! writes the timing values to the output file
//! if the last parameter is set, then all the unused epochs are also printed

void timing_print( int last ) {

  char line[256];
  int i;

  char comparand[20];
  snprintf( comparand, 20, "NA" );

  char eventcount1[20];
  char eventcount2[20];
  char time_str[20];

  for( i=1; i<counter ; i++ ) {

//! calculate the time difference in usecs between timestamp i and the one before that
#if TEST_TYPE != 2
    uint64_t ticks_diff;
    HRT_GET_ELAPSED_TICKS(timings[i-1], timings[i], &ticks_diff);
    double tdiff = HRT_GET_USEC(ticks_diff);
    snprintf( time_str, 20, "%20.4f", tdiff );
#else
    snprintf( time_str, 20, "NA" );
#endif
// select the right values for PAPI columns
    snprintf( eventcount1, 20, "%lli",  PapiEventCounts[2*i] );
    if (strncmp(PapiEventName1, comparand, 20) == 0 ) { /* handling of event counter that are not set */
      snprintf( eventcount1, 20, "NA");
    }
    snprintf( eventcount2, 20, "%lli",  PapiEventCounts[2*i+1] );
    if (strncmp(PapiEventName2, comparand, 20) == 0 ) { /* handling of event counter that are not set */
      snprintf( eventcount2, 20, "NA");
    }

//! if a epoch id is not in 1..5 is found, only the id is printed (instead of
//! the the epoch name)
    if ( (ids[i] > 5) || (ids[i] < 1) ) {
      snprintf(line, 256, "%30s%30s%15i%20i%20s%20s%20s%20s%20s\n", testname, method, bytes, ids[i], time_str, PapiEventName1, eventcount1, PapiEventName2, eventcount2 );
    } else {
      snprintf(line, 256, "%30s%30s%15i%20s%20s%20s%20s%20s%20s\n", testname, method, bytes, epoch_name[ids[i]-1], time_str, PapiEventName1, eventcount1, PapiEventName2, eventcount2 );
//! keep tracking if a epoch is used or not
      if (epoch_used[ids[i]-1] == 0) {
        epoch_used[ids[i]-1] = 1;
      }
    }
    MPI_File_write(filehandle_values, line, strlen(line), MPI_CHAR, MPI_STATUS_IGNORE);
  }

  if (last == 1) {
#if TEST_TYPE > 1
    long long PapiValue[PapiNumberEvents];
    //stop the PAPI counter
    PAPI_stop( PapiEventSet, PapiValue );
    PAPI_reset( PapiEventSet );
#endif
//! each epoch must appear for all tests to make the analysis afterwards easier
//! at the end of the test, we print for each unused epoch a line with a timing of zero micro seconds
    for( i=0; i<5 ; i++ ) {
      if (epoch_used[i] == 0) {
        
        /* handling of event counter that are not set */
        if (strncmp(PapiEventName1, comparand, 20) == 0 ) {
          snprintf( eventcount1, 20, "NA");
        } else {
          snprintf( eventcount1, 20, "0");
        }
        if (strncmp(PapiEventName2, comparand, 20) == 0 ) {
          snprintf( eventcount2, 20, "NA");
        } else {
          snprintf( eventcount2, 20, "0");
        }

#if TEST_TYPE != 2 
        snprintf( time_str, 20, "%20.4f", 0.0 ); 
#else
        snprintf( time_str, 20, "NA" ); // don't record time
#endif

        snprintf(line, 256, "%30s%30s%15i%20s%20s%20s%20s%20s%20s\n", testname, method, bytes, epoch_name[i], time_str, PapiEventName1, eventcount1, PapiEventName2, eventcount2);
        MPI_File_write(filehandle_values, line, strlen(line), MPI_CHAR, MPI_STATUS_IGNORE);
        current_tests++;
      }
    }
  }

//! prints some progression information to stdout
  current_tests = current_tests + counter - 1;
  printf("Finished %7.3f%% of all tests\n", (double) current_tests/max_tests*100);
}

//! records the timing for the given epoch id
//! if number of recorded timings exceed 1024, then the function empties
//! the buffer and they will be printed

void timing_record( int id ) {

#if TEST_TYPE > 1
  PAPI_stop( PapiEventSet, &PapiEventCounts[2*counter] );
  PAPI_reset( PapiEventSet );
#endif

#if TEST_TYPE != 2
  HRT_GET_TIMESTAMP(timings[counter]);
#endif

  ids[counter++] = id;
  if ( counter > 1024 ) {
    timing_print(0);
    ids[0] = 0;
#if TEST_TYPE != 2
    HRT_GET_TIMESTAMP(timings[0]);
#endif
    counter = 1;
  }
#if TEST_TYPE > 1
  PAPI_start( PapiEventSet );
#endif
}

//! opens a file handle, to where the timing values are written

void timing_open_file( char* filename ) {
      
  char line[256];
  int ier;

  ier = MPI_File_open( MPI_COMM_SELF, filename, MPI_MODE_EXCL + MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, &filehandle_values );

  if ( ier != MPI_SUCCESS ) {
    printf("Error at open file %s for writing the timing values. The file probably already exists.\nWill now abort.\n", filename);
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }
  snprintf(line, 256, "%30s%30s%15s%20s%20s%20s%20s%20s%20s\n", "testname", "method", "bytes", "id", "time", "papi_evt1_type", "papi_evt1_val", "papi_evt2_type", "papi_evt2_val");
  MPI_File_write(filehandle_values, line, strlen(line), MPI_CHAR, MPI_STATUS_IGNORE);
}

#if TEST_TYPE > 1
void init_papi() {

  int retval;

  retval = PAPI_library_init(PAPI_VER_CURRENT);
  if (retval != PAPI_VER_CURRENT && retval > 0) {
    fprintf(stderr,"PAPI library version mismatch!\n");
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }

  if (retval < 0) {
    fprintf(stderr,"PAPI init error!\n");
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }

  retval = PAPI_is_initialized();
  if (retval != PAPI_LOW_LEVEL_INITED) {
    fprintf(stderr,"PAPI library not initialized!\n");
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }

  PapiEventSet = PAPI_NULL;
  if (PAPI_create_eventset(&PapiEventSet) != PAPI_OK) {
    fprintf(stderr,"PAPI could not create event set!\n");
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }

  //get the PAPI events
  char* papi_evt1 = getenv("PAPI_EVT1");
  char* papi_evt2 = getenv("PAPI_EVT2");

  PapiNumberEvents = 0;

  if ( (papi_evt1 != NULL) && (strlen(papi_evt1) != 0)) {
    retval = PAPI_event_name_to_code( papi_evt1, &PapiEventCode[0] );
    if ( retval != PAPI_OK ) {
      switch( retval ) {
        case PAPI_ENOTPRESET:   fprintf(stderr,"The event code provided in the environmental variable PAPI_EVT1 is not a valid PAPI preset!\n"); break;
        case PAPI_ENOEVNT:      fprintf(stderr,"The event provided in the environmental variable PAPI_EVT1 is not available on the underlying hardware!\n"); break;
        default:                fprintf(stderr,"The event code provided in the environmental variable PAPI_EVT1 is not correct!\n"); 
      }
      MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    strncpy( PapiEventName1, papi_evt1, 20 );

    PapiNumberEvents++;
    
    if (PAPI_add_event(PapiEventSet, PapiEventCode[0]) != PAPI_OK) {
      fprintf(stderr,"PAPI could not add counter to event set!\n");
      MPI_Abort( MPI_COMM_WORLD, 1 );
    }
  }

  if ( (papi_evt2 != NULL) && (strlen(papi_evt2) != 0) ) {
    retval = PAPI_event_name_to_code( papi_evt2, &PapiEventCode[1] );
    if ( retval != PAPI_OK ) {
      switch( retval ) {
        case PAPI_ENOTPRESET:   fprintf(stderr,"The event code provided in the environmental variable PAPI_EVT2 is not a valid PAPI preset!\n"); break;
        case PAPI_ENOEVNT:      fprintf(stderr,"The event provided in the environmental variable PAPI_EVT2 is not available on the underlying hardware!\n"); break;
        default:                fprintf(stderr,"The event code provided in the environmental variable PAPI_EVT2 is not correct!\n"); 
      }
      MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    strncpy( PapiEventName2, papi_evt2, 20 );

    PapiNumberEvents++;
  
    if ( PAPI_add_event(PapiEventSet, PapiEventCode[1]) != PAPI_OK) {
      fprintf(stderr,"PAPI could not add counter to event set!\n");
      MPI_Abort( MPI_COMM_WORLD, 1 );
    }
  } 

  if ( PapiNumberEvents == 0 ) {
    fprintf(stderr, "Althrough the TEST_TYPE parameter in Makefile.inc indicated a measurement with PAPI counters, the environmental variables PAPI_EVT1 and PAPI_EVT2 are not set.\n");
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }
}

void cleanup_papi() {
  if ( PAPI_cleanup_eventset(PapiEventSet) != PAPI_OK ) {
    fprintf(stderr,"PAPI could not clean up event set!\n");
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }
  if ( PAPI_destroy_eventset(&PapiEventSet) != PAPI_OK ) {
    fprintf(stderr,"PAPI could not destroy event set!\n");
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }
}
#endif

//! closes the above mentioned filehandle
void timing_close_file() {
  MPI_File_close( &filehandle_values );
}

//! sets the initial parameter for the progression bar
void timing_set_max_tests( int value ) {

  max_tests = value;
  current_tests = 0;

}

// fortran wrappers (one underscore)

void timing_init_( char* ptestname, char* pmethod, int* pbytes, int ptestname_len, int pmethod_len ) {

	int i;

	char* c_ptestname = malloc(ptestname_len+1);
	memcpy(c_ptestname, ptestname, ptestname_len);
	c_ptestname[ptestname_len] = '\0';

	for (i=ptestname_len-1; i>0; i--) {
		if (c_ptestname[i] == ' ') {
			c_ptestname[i] = '\0';
		}
		else {
			break;
		}
	}

	char* c_pmethod = malloc(pmethod_len+1);
	memcpy(c_pmethod, pmethod, pmethod_len);
	c_pmethod[pmethod_len] = '\0';

	for (i=pmethod_len-1; i>0; i--) {
		if (c_pmethod[i] == ' ') {
			c_pmethod[i] = '\0';
		}
		else {
			break;
		}
	}

	timing_init(c_ptestname, c_pmethod, *pbytes);

	free(c_ptestname);
	free(c_pmethod);

}

void timing_print_( int* last ) {
	timing_print(*last);
}

void timing_record_( int* id ) {
	timing_record(*id);
}

void timing_open_file_( char* filename, int len ) {

	int i;
	char* c_filename = malloc(len+1);
	memcpy(c_filename, filename, len);
	c_filename[len] = '\0';

	for (i=len-1; i>0; i--) {
		if (c_filename[i] == ' ') {
			c_filename[i] = '\0';
		}
		else {
			break;
		}
	}

	timing_open_file(c_filename);

	free(c_filename);

}

void timing_close_file_() {
	timing_close_file();
}

void timing_set_max_tests_( int* value ) {
	timing_set_max_tests(*value);
}

void timing_hrt_init_() {
  timing_hrt_init();
}

#if TEST_TYPE > 1
void init_papi_() {
  init_papi();
}

void cleanup_papi_() {
  cleanup_papi();
}
#endif

// fortran wrappers (two underscores)

void timing_init__( char* ptestname, int ptestname_len, char* pmethod, int pmethod_len, int* pbytes ) {

	int i;

	char* c_ptestname = malloc(ptestname_len+1);
	memcpy(c_ptestname, ptestname, ptestname_len);
	c_ptestname[ptestname_len] = '\0';

	for (i=ptestname_len-1; i>0; i--) {
		if (c_ptestname[i] == ' ') {
			c_ptestname[i] = '\0';
		}
		else {
			break;
		}
	}

	char* c_pmethod = malloc(pmethod_len+1);
	memcpy(c_pmethod, pmethod, pmethod_len);
	c_pmethod[pmethod_len] = '\0';

	for (i=pmethod_len-1; i>0; i--) {
		if (c_pmethod[i] == ' ') {
			c_pmethod[i] = '\0';
		}
		else {
			break;
		}
	}

	timing_init(c_ptestname, c_pmethod, *pbytes);

	free(c_ptestname);
	free(c_pmethod);
}

void timing_print__( int* last ) {
	timing_print(*last);
}

void timing_record__( int* id ) {
	timing_record(*id);
}

void timing_open_file__( char* filename, int len ) {

	int i;
	char* c_filename = malloc(len+1);
	memcpy(c_filename, filename, len);
	c_filename[len] = '\0';

	for (i=len-1; i>0; i--) {
		if (c_filename[i] == ' ') {
			c_filename[i] = '\0';
		}
		else {
			break;
		}
	}

	timing_open_file(c_filename);

	free(c_filename);

}

void timing_close_file__() {
	timing_close_file();
}

void timing_set_max_tests__( int* value ) {
	timing_set_max_tests(*value);
}

void timing_hrt_init__() {
  timing_hrt_init();
}
#if TEST_TYPE > 1
void init_papi__() {
  init_papi();
}

void cleanup_papi__() {
  cleanup_papi();
}
#endif
