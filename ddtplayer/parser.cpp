/* A Bison parser, made by GNU Bison 2.5.  */

/* Bison implementation for Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2011 Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.
   
   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.5"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Copy the first part of user declarations.  */

/* Line 268 of yacc.c  */
#line 1 "parser.y"

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



/* Line 268 of yacc.c  */
#line 119 "parser.cpp"

/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     NUM = 258,
     UNKNOWN = 259,
     SUBTYPE = 260,
     ELEM = 261,
     CONTIGUOUS = 262,
     VECTOR = 263,
     HVECTOR = 264,
     HINDEXED = 265,
     STRUCT = 266,
     BYTE_ = 267,
     CHAR_ = 268,
     INT_ = 269,
     DOUBLE_ = 270,
     FLOAT_ = 271
   };
#endif



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 293 of yacc.c  */
#line 48 "parser.y"

	int val;
	struct Datatypes *types;

	struct Index *index;
	struct Indices *indices;

	struct {
		int start;
		int stop;
		int stride;
	} range;



/* Line 293 of yacc.c  */
#line 187 "parser.cpp"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 343 of yacc.c  */
#line 199 "parser.cpp"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   53

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  23
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  13
/* YYNRULES -- Number of rules.  */
#define YYNRULES  24
/* YYNRULES -- Number of states.  */
#define YYNSTATES  57

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   271

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      18,    19,     2,     2,    22,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    17,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    20,     2,    21,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,     4,     7,     9,    11,    13,    15,    17,
      19,    21,    23,    25,    27,    29,    31,    33,    39,    47,
      57,    67,    71,    72,    75
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      24,     0,    -1,    -1,    24,    25,    -1,    26,    -1,    27,
      -1,    28,    -1,    12,    -1,    13,    -1,    14,    -1,    16,
      -1,    15,    -1,    30,    -1,    31,    -1,    32,    -1,    35,
      -1,     3,    -1,     3,    17,     3,    17,     3,    -1,     7,
      18,    29,    19,    20,    26,    21,    -1,     8,    18,    29,
      29,    29,    19,    20,    26,    21,    -1,     9,    18,    29,
      29,    29,    19,    20,    26,    21,    -1,     3,    22,     3,
      -1,    -1,    34,    33,    -1,    10,    18,    34,    19,    20,
      26,    21,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint8 yyrline[] =
{
       0,    76,    76,    77,    81,    88,    89,    93,   100,   107,
     114,   121,   131,   132,   133,   134,   138,   143,   152,   173,
     198,   223,   231,   234,   241
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "NUM", "UNKNOWN", "SUBTYPE", "ELEM",
  "CONTIGUOUS", "VECTOR", "HVECTOR", "HINDEXED", "STRUCT", "BYTE_",
  "CHAR_", "INT_", "DOUBLE_", "FLOAT_", "':'", "'('", "')'", "'['", "']'",
  "','", "$accept", "input", "topdatatype", "datatype", "primitive",
  "derived", "range", "contiguous", "vector", "hvector", "idxentry",
  "idxentries", "hindexed", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,    58,    40,    41,
      91,    93,    44
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    23,    24,    24,    25,    26,    26,    27,    27,    27,
      27,    27,    28,    28,    28,    28,    29,    29,    30,    31,
      32,    33,    34,    34,    35
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     5,     7,     9,
       9,     3,     0,     2,     7
};

/* YYDEFACT[STATE-NAME] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,     0,     0,     0,     0,     7,     8,     9,
      11,    10,     3,     4,     5,     6,    12,    13,    14,    15,
       0,     0,     0,    22,    16,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    23,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    21,     0,    17,    18,
       0,     0,    24,     0,     0,    19,    20
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     1,    12,    13,    14,    15,    25,    16,    17,    18,
      35,    28,    19
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -20
static const yytype_int8 yypact[] =
{
     -20,     0,   -20,   -16,   -15,   -14,   -13,   -20,   -20,   -20,
     -20,   -20,   -20,   -20,   -20,   -20,   -20,   -20,   -20,   -20,
       3,     3,     3,   -20,    -6,     1,     3,     3,    -2,    16,
       5,     3,     3,    -1,     6,   -20,    10,    28,    11,    20,
      42,    28,    43,    26,    29,    30,   -20,    27,   -20,   -20,
      28,    28,   -20,    31,    32,   -20,   -20
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -20,   -20,   -20,   -19,   -20,   -20,     2,   -20,   -20,   -20,
     -20,   -20,   -20
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
       2,    33,    20,    21,    22,    23,    24,     3,     4,     5,
       6,    29,     7,     8,     9,    10,    11,    34,    43,    36,
      30,    40,    47,    26,    27,    37,    41,    42,    31,    32,
      44,    53,    54,    38,    39,     3,     4,     5,     6,    45,
       7,     8,     9,    10,    11,    46,    48,    49,    52,    50,
      51,     0,    55,    56
};

#define yypact_value_is_default(yystate) \
  ((yystate) == (-20))

#define yytable_value_is_error(yytable_value) \
  YYID (0)

static const yytype_int8 yycheck[] =
{
       0,     3,    18,    18,    18,    18,     3,     7,     8,     9,
      10,    17,    12,    13,    14,    15,    16,    19,    37,     3,
      19,    22,    41,    21,    22,    20,    20,    17,    26,    27,
      19,    50,    51,    31,    32,     7,     8,     9,    10,    19,
      12,    13,    14,    15,    16,     3,     3,    21,    21,    20,
      20,    -1,    21,    21
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    24,     0,     7,     8,     9,    10,    12,    13,    14,
      15,    16,    25,    26,    27,    28,    30,    31,    32,    35,
      18,    18,    18,    18,     3,    29,    29,    29,    34,    17,
      19,    29,    29,     3,    19,    33,     3,    20,    29,    29,
      22,    20,    17,    26,    19,    19,     3,    26,     3,    21,
      20,    20,    21,    26,    26,    21,    21
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  However,
   YYFAIL appears to be in use.  Nevertheless, it is formally deprecated
   in Bison 2.4.2's NEWS entry, where a plan to phase it out is
   discussed.  */

#define YYFAIL		goto yyerrlab
#if defined YYFAIL
  /* This is here to suppress warnings from the GCC cpp's
     -Wunused-macros.  Normally we don't worry about that warning, but
     some users do, and we want to make it easy for users to remove
     YYFAIL uses, which will produce warnings from Bison 2.5.  */
#endif

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* This macro is provided for backward compatibility. */

#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (0, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  YYSIZE_T yysize1;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = 0;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - Assume YYFAIL is not used.  It's too flawed to consider.  See
       <http://lists.gnu.org/archive/html/bison-patches/2009-12/msg00024.html>
       for details.  YYERROR is fine as it does not invoke this
       function.
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                yysize1 = yysize + yytnamerr (0, yytname[yyx]);
                if (! (yysize <= yysize1
                       && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                  return 2;
                yysize = yysize1;
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  yysize1 = yysize + yystrlen (yyformat);
  if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
    return 2;
  yysize = yysize1;

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */
#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */


/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks thru separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yytoken = 0;
  yyss = yyssa;
  yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */
  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 4:

/* Line 1806 of yacc.c  */
#line 81 "parser.y"
    {
	datatypes.insert( datatypes.end(), (yyvsp[(1) - (1)].types)->types.begin(), (yyvsp[(1) - (1)].types)->types.end() );
	free((yyvsp[(1) - (1)].types));
}
    break;

  case 7:

/* Line 1806 of yacc.c  */
#line 93 "parser.y"
    {
	(yyval.types) = new Datatypes;
	struct Datatype datatype;
	datatype.farc = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::BYTE);
	datatype.mpi  = MPI_BYTE;
	(yyval.types)->types.push_back(datatype);
}
    break;

  case 8:

/* Line 1806 of yacc.c  */
#line 100 "parser.y"
    {
	(yyval.types) = new Datatypes;
	struct Datatype datatype;
	datatype.farc = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::CHAR);
	datatype.mpi  = MPI_CHAR;
	(yyval.types)->types.push_back(datatype);
}
    break;

  case 9:

/* Line 1806 of yacc.c  */
#line 107 "parser.y"
    {
	(yyval.types) = new Datatypes;
	struct Datatype datatype;
	datatype.farc = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::INT);
	datatype.mpi  = MPI_INT;
	(yyval.types)->types.push_back(datatype);
}
    break;

  case 10:

/* Line 1806 of yacc.c  */
#line 114 "parser.y"
    {
	(yyval.types) = new Datatypes;
	struct Datatype datatype;
	datatype.farc = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::FLOAT);
	datatype.mpi  = MPI_FLOAT;
	(yyval.types)->types.push_back(datatype);
}
    break;

  case 11:

/* Line 1806 of yacc.c  */
#line 121 "parser.y"
    {
	(yyval.types) = new Datatypes;
	struct Datatype datatype;
	datatype.farc = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::DOUBLE);
	datatype.mpi  = MPI_DOUBLE;
	(yyval.types)->types.push_back(datatype);
}
    break;

  case 16:

/* Line 1806 of yacc.c  */
#line 138 "parser.y"
    {
	(yyval.range).start = (yyvsp[(1) - (1)].val);
	(yyval.range).stride = 1;
	(yyval.range).stop = (yyvsp[(1) - (1)].val);
}
    break;

  case 17:

/* Line 1806 of yacc.c  */
#line 143 "parser.y"
    {
	(yyval.range).start = (yyvsp[(1) - (5)].val);
	(yyval.range).stride = (yyvsp[(3) - (5)].val);
	(yyval.range).stop = (yyvsp[(5) - (5)].val);
}
    break;

  case 18:

/* Line 1806 of yacc.c  */
#line 152 "parser.y"
    {
	(yyval.types) = new Datatypes;

	list<struct Datatype> *subtypes = &((yyvsp[(6) - (7)].types)->types);
	list<struct Datatype> *types = &((yyval.types)->types);

	for(list<struct Datatype>::iterator subtype = subtypes->begin();
		subtype != subtypes->end(); subtype++) {
		for (int count=(yyvsp[(3) - (7)].range).start; count <= (yyvsp[(3) - (7)].range).stop; count += (yyvsp[(3) - (7)].range).stride) {
			Datatype type;
			type.farc = new farc::ContiguousDatatype(subtype->farc, count);
			MPI_Type_contiguous(count, subtype->mpi, &(type.mpi));
			types->push_back(type);
		}
	}

	free((yyvsp[(6) - (7)].types));
}
    break;

  case 19:

/* Line 1806 of yacc.c  */
#line 173 "parser.y"
    { 
	(yyval.types) = new Datatypes;

	list<struct Datatype> *subtypes = &((yyvsp[(8) - (9)].types)->types);
	list<struct Datatype> *types = &((yyval.types)->types);

	for(list<struct Datatype>::iterator subtype = subtypes->begin();
		subtype != subtypes->end(); subtype++) {
		for (int count=(yyvsp[(3) - (9)].range).start; count <= (yyvsp[(3) - (9)].range).stop; count += (yyvsp[(3) - (9)].range).stride) {
			for (int blocklen=(yyvsp[(4) - (9)].range).start; blocklen <= (yyvsp[(4) - (9)].range).stop; blocklen += (yyvsp[(4) - (9)].range).stride) {
				for (int stride=(yyvsp[(5) - (9)].range).start; stride <= (yyvsp[(5) - (9)].range).stop; stride += (yyvsp[(5) - (9)].range).stride) {
					Datatype type;
					type.farc = new farc::VectorDatatype(subtype->farc, count, blocklen, stride);
					MPI_Type_vector(count, blocklen, stride, subtype->mpi, &(type.mpi));
					types->push_back(type);
				}
			}
		}
	}

	free((yyvsp[(8) - (9)].types));
}
    break;

  case 20:

/* Line 1806 of yacc.c  */
#line 198 "parser.y"
    {
	(yyval.types) = new Datatypes;

	list<struct Datatype> *subtypes = &((yyvsp[(8) - (9)].types)->types);
	list<struct Datatype> *types = &((yyval.types)->types);

	for(list<struct Datatype>::iterator subtype = subtypes->begin();
		subtype != subtypes->end(); subtype++) {
		for (int count=(yyvsp[(3) - (9)].range).start; count <= (yyvsp[(3) - (9)].range).stop; count += (yyvsp[(3) - (9)].range).stride) {
			for (int blocklen=(yyvsp[(4) - (9)].range).start; blocklen <= (yyvsp[(4) - (9)].range).stop; blocklen += (yyvsp[(4) - (9)].range).stride) {
				for (int stride=(yyvsp[(5) - (9)].range).start; stride <= (yyvsp[(5) - (9)].range).stop; stride += (yyvsp[(5) - (9)].range).stride) {
					Datatype type;
					type.farc = new farc::HVectorDatatype(subtype->farc, count, blocklen, stride);
					MPI_Type_hvector(count, blocklen, stride, subtype->mpi, &(type.mpi));
					types->push_back(type);
				}
			}
		}
	}

	free((yyvsp[(8) - (9)].types));
}
    break;

  case 21:

/* Line 1806 of yacc.c  */
#line 223 "parser.y"
    {
	(yyval.index) = new Index;
	(yyval.index)->displ = (yyvsp[(1) - (3)].val);
	(yyval.index)->blocklen = (yyvsp[(3) - (3)].val);
}
    break;

  case 22:

/* Line 1806 of yacc.c  */
#line 231 "parser.y"
    {
	(yyval.indices) = new Indices;
}
    break;

  case 23:

/* Line 1806 of yacc.c  */
#line 234 "parser.y"
    {
	(yyval.indices) = (yyvsp[(1) - (2)].indices);
	(yyval.indices)->indices.push_back((yyvsp[(2) - (2)].index));
}
    break;

  case 24:

/* Line 1806 of yacc.c  */
#line 241 "parser.y"
    {
	(yyval.types) = new Datatypes;

	list<struct Datatype> *subtypes = &((yyvsp[(6) - (7)].types)->types);
	list<struct Datatype> *types = &((yyval.types)->types);

	unsigned int num = (yyvsp[(3) - (7)].indices)->indices.size();
	long *displs = (long*)malloc(num * sizeof(long));
	int *blocklens = (int*)malloc(num * sizeof(int));

	for(int i=0; i<num; i++) {
		displs[i] = (yyvsp[(3) - (7)].indices)->indices[i]->displ;
		blocklens[i] = (yyvsp[(3) - (7)].indices)->indices[i]->blocklen;
		free((yyvsp[(3) - (7)].indices)->indices[i]);
	}
	free((yyvsp[(3) - (7)].indices));

	for(list<struct Datatype>::iterator subtype = subtypes->begin();
		subtype != subtypes->end(); subtype++) {
		Datatype type;
		type.farc = new farc::HIndexedDatatype(num, blocklens, displs, subtype->farc);
		MPI_Type_hindexed(num, blocklens, displs, subtype->mpi, &(type.mpi));
		types->push_back(type);
	}

	free(displs);
	free(blocklens);
	free((yyvsp[(6) - (7)].types));
}
    break;



/* Line 1806 of yacc.c  */
#line 1701 "parser.cpp"
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined(yyoverflow) || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}



/* Line 2067 of yacc.c  */
#line 272 "parser.y"


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
            printf("EXTENT MISSMATCH: MPI extent: %i FARC extent: %i\n", extent_mpi, extent);
		    cout << setw(name_w)        << datatype.farc->toString().c_str() << std::endl;
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

