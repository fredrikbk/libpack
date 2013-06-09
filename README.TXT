LibPack: Runtime-compiled pack functions (not only) for MPI Datatypes
=====================================================================

Introduction:

    The libpack library lets the user build a description of an arbitrary
    data layout. The API is very similar to that of MPI derived datatypes
    (libpack also includes a wrapper which interposes MPI calls and calls
    libpack functions instead). After the data layout description is complete,
    it has to be compiled by calling a commit-function. Data can only be
    packed with datatypes which have been committed. However, using a datatype
    in the construction of another datatype does not require it to be
    committed.

    During the commit function call the libpack library generates native
    machine code which is capable of packing and unpacking data with the
    specified layout. This is different to what the MPI implementations in use
    today do: they interpret the datatype instead. This means they built an
    optimized in-memory representation of the data-layout, which is fed to an
    interpreter to actually pack / unpack data with the specified layout. This
    means that they do not have to pay much overhead at the commit time for
    code generation, on the other hand, interpretation is often slower than
    native execution (to mitigate this many interpreted languages also
    runtime-compile "hot" functions, i.e., functions which dominate the
    execution time).

    Our library is meant to be a research vehicle to experiment with such
    techniques. In order to allow other to seamlessly integrate it into their
    own work and modify it according to their needs we will describe the
    current architecture of our library in the following section.

Contents:

    The library consists of several parts:

    - a core, which is responsible for generating pack functions. The core is
      written in C++.

      The core is implemented in the files ddt_jit.cpp and the native
      code generation for the different derived datatypes defined by the MPI
      standard is implemented in the codegen_* files.

    - a set of wrappers, which allow others to use this library without
      modifying their code to much. The wrappers use some common
      functionality, for example to track the usage of temporary buffers and
      datatype handles, which are implemented in interposer_common.cpp

        - a C interface (which can also be called from FORTRAN and other
          languages, as it has only primitive types as function arguments and
          return values). The C interface is defined and implemented in pack.h
          and pack.cpp.

        - an MPI wrapper for most of the datatype related MPI functions. This
          allows to use an unmodified MPI program and substitute MPIs
          datatype processing with the one offered by libpack. When doing so
          please make sure that the functions used by the MPI program are
          actually implemented in libpack. Currently, we do, for example,
          not provide wrappers for MPI one-sided. To use the MPI wrapper,
          simply link libpack.a into your program. Be aware that the order in
          which you link in libraries matters. When using GNU ld, you should
          put libfarc.a at the end of your linker flags. Note that using this
          wrapper is not as efficient as the use of libpack could be, as the
          interface defined by MPI forces us to do some expensive lookups and
          checks: for example when we wrap a non-blocking receive with a
          derived datatype we have to dynamically allocate a temporary
          receive / unpack buffer, and subsequently check in every MPI Test(),
          Wait(), Waitany(), etc. call if this receive is finished now, and if
          yes, free the temporary resources This overhead would be much 
          smaller if MPI would allow the user to attach private data to MPI
          requests. However, this is not done to keep requests as small as
          possible. The MPI wrappers are implemented in interposer.c and of
          course they do not need their own header file, as they use the
          function declarations of mpi.h.

    - a set of tests, including a test harness which makes it easy to add your
      own tests. The tests can be found in the directories "tests", where
      tests for the C and C++ interface are located, and in the directory
      pmpi-tests where tests for the MPI wrappers are put in.

    - a small datatype benchmarking program called "ddtplayer". DDTPlayer is
      capable of parsing textual descriptions of datatypes, such as 

                hidx(0,1 17952,1)[vec(34:10:64 1 34)[double]]

      and construct the described types in libpack as well as their MPI DDT
      equivalent. It then benchmarks the packing and unpacking performance of
      both and prints the absolute times as well as the speedup of libpack
      over MPI. The line above describes the following datatypes:

      An hindexed datatype with displacements 0, 17952 and counts 1,1 with
      different vectors as basetypes: The count of the vectors starts at 34
      and ends at 64 and increases in steps of size 10. the basetype of the
      vectors is a double.

Implementation:

    In this section we will describe the implementation of the core in detail.
    A derived datatype is constructed as a tree of C++ classes. The leaves of
    such a datatype tree are primitive datatypes, such as int, doubles, floats,
    etc. They can be constructed using libpack by calling 

    farc::Datatype* t1;
    t1 = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::INT);

    which would make t1 a primitive datatype of the integer type. Building a
    vector of those integer types can now be done with the line

    farc::Datatype* t2 = new farc::VectorDatatype(2, 3, 5, t1);

    which creates a vector with count=2, blocklen=3 and stride=5. For the
    semantics of the different datatypes, as well as the meaning and order of
    their arguments, please consult the MPI standard. Internally, this would
    create the following datastructure:

    VectorDatatype
        count = 2
        blocklen = 3
        stride = 5
        lower_bound = 0
        upper_bound = 32
        size = 24;
        basetype = *PrimitiveDatatype
                      size = 4
                      lower_bound = 0
                      upper_bound = 4

    Note that the constructed datatype retains all the information provided
    during the construction of the MPI datatype. This is necessary to
    implement equivalents to the MPI functions MPI_Type_get_envelope and
    MPI_Type_get_contents.

    Each Datatype class contains a clone() function which deep-copies a
    datatype This is necessary for the construction of nested datatypes,
    as when constructing the vector in the example above, we can not simply
    use the provided pointer to the primitive type and store it as the vectors
    basetype, as it is possible that the user frees the primitive datatype
    before committing the vector. Therefore, all datatypes provided as
    arguments to the constructors of other datatypes are cloned, and the
    cloned version is used internally.

    Each Datatype object needs to implement the packCodegen and unpackCodegen
    method. Those methods emit the LLVM IR for packing the datatype. An
    explanation of LLVM IR is out of the scope of this document. The
    interested reader is referred to the document "Kaleidoscope: Implementing
    a Language with LLVM" which comes with LLVM and explains how to implement
    a JIT-compiled language with the same method we use for libpack. The
    function generated by codegenPack() looks like this (in pseudocode):

    pack(int8_t* inbuf, int count, int8_t* outbuf) {

        // the generated code may not modify its arguments so that it
        // is fully composable with other generated code parts

        out1 = outbuf;
        in1 = inbuf;

        // outer loop over datatype count (as supplied to pack or mpi_send)
        for (int i=0; i<count; i++) {

            nextout1 = out1 + this->size;
            nextin1 = NULL;

            // inner loop over vector count
            for (;;) {

                // here we generate code to pack basetype
                // basetype->packCodegen(in2, this->blocklen, out2)
            
                out2 = out2 + this->basetype->size * this->blocklen;
                in2 = in2 + this->basetype->extent * this->stride;

                if (out2 == out1) break;
    
            }
            nextin1 = in1 + this->extent;
            out1 = nextout1;
            in1 = nextin1;

        }

    }

    This function is generated by in codegenVector() in codegen_vector.cpp.
    Note that the LLVM IR which is used to generate pack functions uses single
    static assignment notation. This means there are no mutable variables.
    Each variable is assigned exactly once in its lifetime. However, where
    mutable variables are needed in traditional programs (i.e, counter
    variables in loops), phi-nodes can be used. A phi node can have a
    different value, depending on the control flow, i.e., the out1 variable in
    the above pseudocode is such a phi node - if it is reached through the
    loop header (initial iteration of the loop) its value is outbuf, but if it
    is reached through the body of the loop its value is nextout1.

    Note that the pack and unpack functions for most datatypes are similar
    enough that we can generate the actual pack/unpack code with the same
    function, using a few if-statements in the common generation function. The
    same is true for some datatypes: The indexed and hindexed datatype for
    example use a common function to generate their packing and unpacking
    functions.

    When the pack/unpack functions have been generated, LLVM JIT compiler
    gives us a normal function pointer, which is attached to the C++ object
    which represents the datatype.

Future Work:

    There are still some areas where libpack could be improved.
    We are currently researching how to optimize pipelining of packing and
    sending, as well as how to generate an optimal instruction sequence for
    the packing. Other ideas are to improve the code generation, perhaps by
    including our own IR for datatypes and to include other scenarios where
    packing / serializing data is of interest (i.e., host to accelerator
    communication, C++ object serialization, etc.).

Citation:

    If you use libpack for your work, please cite

    --------------------------------------------------
    T. Schneider, F. Kjolstad and T. Hoefler
    MPI Datatype Processing using Runtime Compilation
    In the proceeding of EuroMPI'13, Madrid, Spain
    Springer, LNCS, 2013
    --------------------------------------------------
