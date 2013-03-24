#ifndef FARC_INCLUDED
#define FARC_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

typedef void* FARC_Datatype;
typedef int   FARC_Primitive;

FARC_Datatype FARC_CHAR; 
FARC_Datatype FARC_SIGNED_CHAR;
FARC_Datatype FARC_UNSIGNED_CHAR;
FARC_Datatype FARC_BYTE;
FARC_Datatype FARC_WCHAR;
FARC_Datatype FARC_SHORT;
FARC_Datatype FARC_UNSIGNED_SHORT;
FARC_Datatype FARC_INT;
FARC_Datatype FARC_UNSIGNED;
FARC_Datatype FARC_LONG;
FARC_Datatype FARC_UNSIGNED_LONG;
FARC_Datatype FARC_FLOAT;
FARC_Datatype FARC_DOUBLE;
FARC_Datatype FARC_LONG_DOUBLE;
FARC_Datatype FARC_LONG_LONG_INT;
FARC_Datatype FARC_UNSIGNED_LONG_LONG;
FARC_Datatype FARC_LONG_LONG;

typedef long FARC_Aint;


/* Functions */
int FARC_Init();
int FARC_Finalize();

int FARC_Contiguous(int count, FARC_Datatype old_type, FARC_Datatype *newtype_p);
int FARC_Vector(int count, int blocklen, int stride, FARC_Datatype oldtype, FARC_Datatype *newtype_p);
int FARC_Hvector(int count, int blocklen, FARC_Aint stride, FARC_Datatype oldtype, FARC_Datatype *newtype_p);
int FARC_Indexed_block(int count, int blocklen, int displacements[], FARC_Datatype oldtype, FARC_Datatype *newtype_p);
int FARC_Hindexed(int count, int blocklens[], FARC_Aint displacements[], FARC_Datatype oldtype, FARC_Datatype *newtype_p);
int FARC_Struct(int count, int blocklens[], FARC_Aint displacements[], FARC_Datatype oldtypes[], FARC_Datatype *newtype_p);

int FARC_Free(FARC_Datatype *ddt);

int FARC_Commit(FARC_Datatype *ddt);
int FARC_Commit_pack(FARC_Datatype *ddt);
int FARC_Commit_unpack(FARC_Datatype *ddt);

int FARC_Pack(void* inbuf, void* outbuf, FARC_Datatype type, int count);
int FARC_Unpack(void* inbuf, void* outbuf, FARC_Datatype type, int count);

int FARC_Get_extent(FARC_Datatype datatype, FARC_Aint *lb, FARC_Aint *extent);
int FARC_Get_size(FARC_Datatype datatype, int *size);

#if defined(__cplusplus)
}
#endif
#endif
