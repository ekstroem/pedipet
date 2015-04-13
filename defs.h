/**********************************************
 *
 * This file sets the appropriate global
 * constants for pedipet.
 *
 * All constants should be preceded by a comment
 * giving a short description of the constant
 *
 * (C) Claus Ekstrøm 2001 
 *
 **********************************************/


#ifndef _DEFS_H_
#define _DEFS_H_

// Version number
#define  VERSION               "0.61a"

//
#define  LINESIZE              79

// The number of characters that can be used for names
// (id, pedigree names etc).
// Curently this also influences the number of markers that
// can be handles (dur to sloppy programming). SHOULD BE NO LESS
// THAN MAXMARKERS
#define  NAMESIZE              251   

//
#define  MAXMARKERS            251

//
#define  MAXALLELES            66        // Maximum number of alleles pr. locus

//
#define  MAXPERSONS            50        // Maximum number of persons in a single pedigree

#define  MAXTRAITS             20
#define  BUFFERSIZE            2048




// Not all compilers have defined pi
#ifndef PI
#define PI 3.141592654
#endif

#define ZERO  0
#define NaN               (0.0f / 0.0f)
#define POSITIVE_INFINITY (1.0f / 0.0f)
#define NEGATIVE_INFINITY (-1.0f / 0.0f)


#endif
