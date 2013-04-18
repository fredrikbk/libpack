// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

#ifndef LPK_INCLUDED
#define LPK_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

typedef void* LPK_Datatype;
typedef int   LPK_Primitivetype;

#define LPK_CHAR                0
#define LPK_SIGNED_CHAR         1
#define LPK_UNSIGNED_CHAR       2
#define LPK_BYTE                3
#define LPK_WCHAR               4
#define LPK_SHORT               5
#define LPK_UNSIGNED_SHORT      6
#define LPK_INT                 7
#define LPK_UNSIGNED            8
#define LPK_LONG                9
#define LPK_UNSIGNED_LONG      10
#define LPK_FLOAT              11
#define LPK_DOUBLE             12
#define LPK_LONG_DOUBLE        13
#define LPK_LONG_LONG_INT      14
#define LPK_UNSIGNED_LONG_LONG 15
#define LPK_LONG_LONG          16

typedef long LPK_Aint;


/* Functions */
int LPK_Init();
int LPK_Finalize();

int LPK_Primitive(LPK_Primitivetype primitivetype, LPK_Datatype *newtype_p);
int LPK_Contiguous(int count, LPK_Datatype oldtype, LPK_Datatype *newtype_p);
int LPK_Vector(int count, int blocklen, int stride, LPK_Datatype oldtype, LPK_Datatype *newtype_p);
int LPK_Hvector(int count, int blocklen, LPK_Aint stride, LPK_Datatype oldtype, LPK_Datatype *newtype_p);
int LPK_Indexed_block(int count, int blocklen, int displacements[], LPK_Datatype oldtype, LPK_Datatype *newtype_p);
int LPK_Hindexed(int count, int blocklens[], LPK_Aint displacements[], LPK_Datatype oldtype, LPK_Datatype *newtype_p);
int LPK_Struct(int count, int blocklens[], LPK_Aint displacements[], LPK_Datatype oldtypes[], LPK_Datatype *newtype_p);

int LPK_Free(LPK_Datatype *ddt);

int LPK_Compile(LPK_Datatype *ddt);
int LPK_Compile_pack(LPK_Datatype *ddt);
int LPK_Compile_unpack(LPK_Datatype *ddt);

int LPK_Pack(void* inbuf, int incount, LPK_Datatype intype, void* outbuf);
int LPK_Unpack(void* inbuf, void* outbuf, int outcount, LPK_Datatype outtype);

int LPK_Get_extent(LPK_Datatype datatype, LPK_Aint *lb, LPK_Aint *extent);
int LPK_Get_size(LPK_Datatype datatype, int *size);

#if defined(__cplusplus)
}
#endif
#endif
