/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

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
/* Tokens.  */
#define NUM 258
#define UNKNOWN 259
#define SUBTYPE 260
#define ELEM 261
#define CONTIGUOUS 262
#define VECTOR 263
#define HVECTOR 264
#define HINDEXED 265
#define STRUCT 266
#define BYTE_ 267
#define CHAR_ 268
#define INT_ 269
#define DOUBLE_ 270
#define FLOAT_ 271




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 37 "parser.y"
{
	int val;
	struct Datatype *datatype;
}
/* Line 1529 of yacc.c.  */
#line 86 "parser.hpp"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

