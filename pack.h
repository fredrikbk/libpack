#ifndef FARC_INCLUDED
#define FARC_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

typedef void* FARC_Datatype;
typedef int   FARC_Primitivetype;

#define FARC_CHAR                0 
#define FARC_SIGNED_CHAR         1
#define FARC_UNSIGNED_CHAR       2
#define FARC_BYTE                3
#define FARC_WCHAR               4
#define FARC_SHORT               5
#define FARC_UNSIGNED_SHORT      6
#define FARC_INT                 7
#define FARC_UNSIGNED            8
#define FARC_LONG                9
#define FARC_UNSIGNED_LONG      10
#define FARC_FLOAT              11
#define FARC_DOUBLE             12
#define FARC_LONG_DOUBLE        13
#define FARC_LONG_LONG_INT      14
#define FARC_UNSIGNED_LONG_LONG 15
#define FARC_LONG_LONG          16

typedef long FARC_Aint;


/* Functions */
int FARC_Init();
int FARC_Finalize();

int FARC_Primitive(FARC_Primitivetype primitivetype, FARC_Datatype *newtype_p);
int FARC_Contiguous(int count, FARC_Datatype oldtype, FARC_Datatype *newtype_p);
int FARC_Vector(int count, int blocklen, int stride, FARC_Datatype oldtype, FARC_Datatype *newtype_p);
int FARC_Hvector(int count, int blocklen, FARC_Aint stride, FARC_Datatype oldtype, FARC_Datatype *newtype_p);
int FARC_Indexed_block(int count, int blocklen, int displacements[], FARC_Datatype oldtype, FARC_Datatype *newtype_p);
int FARC_Hindexed(int count, int blocklens[], FARC_Aint displacements[], FARC_Datatype oldtype, FARC_Datatype *newtype_p);
int FARC_Struct(int count, int blocklens[], FARC_Aint displacements[], FARC_Datatype oldtypes[], FARC_Datatype *newtype_p);

int FARC_Free(FARC_Datatype *ddt);

int FARC_Commit(FARC_Datatype *ddt);
int FARC_Commit_pack(FARC_Datatype *ddt);
int FARC_Commit_unpack(FARC_Datatype *ddt);

int FARC_Pack(void* inbuf, int incount, FARC_Datatype intype, void* outbuf);
int FARC_Unpack(void* inbuf, void* outbuf, int outcount, FARC_Datatype outtype);

int FARC_Get_extent(FARC_Datatype datatype, FARC_Aint *lb, FARC_Aint *extent);
int FARC_Get_size(FARC_Datatype datatype, int *size);

#if defined(__cplusplus)
}
#endif
#endif
