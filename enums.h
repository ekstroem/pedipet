/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
 * This file sets the appropriate global
 * enumerations for pedipet.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 */

#ifndef _ENUMS_H_
#define _ENUMS_H_


#ifdef __cplusplus
extern "C" {
#endif

//
// Possible sex-values
//
enum
{
  S_UNKNOWN,
  S_MALE,
  S_FEMALE,
  MAXSEX,
};


//
// Mapping functions
//
enum
{
  M_HALDANE,
  M_KOSAMBI,
  MAXMAP,
};


enum
{
  COL_BLACK,
  COL_RED,
  COL_GREEN,
  COL_YELLOW,
  COL_BLUE,
  COL_MAGENTA,
  COL_CYAN,
  COL_WHITE,
  MAXCOL,
};



// Possible option keywords
enum
{
  O_COLOUR,
  O_CONVERTLINKAGEID,
  O_DOMINANCE,
  O_IDCASESENSITIVE,
  O_ITERATIONINFO,
  O_NEWERRORFILE,
  O_NOBOUND,
  O_NOSEED,
  O_REDUCEALLELES,
  O_SAVEDIR,
  O_SEED,
  O_USEERRORFILE,
  O_XLINKED,
  MAXOPTIONS,
};



#ifdef __cplusplus
}
#endif


#endif
