/*
 * PEDIPET
 * 
 * cestring.c
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
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


#include <ctype.h>

#include "cestring.h"

// Makes a copy of source without white spaces
void RemoveWhiteSpace(char *dest, char *source)
{
  while(*source) {
    if (!isspace(*source))
      *dest++ = *source;
    source++;
  }
  *dest++ = '\0';
}

// Makes a copy of source with single white spaces
void SimplifyWhiteSpace(char *dest, char *source)
{
  int space = 0;

  while(*source) {
    if (!isspace(*source)) {
      space = 0;
      *dest++ = *source;
    }
    else {
      if (!space)
        *dest++ = *source;
      space = 1;
    }
    source++;
  }
  *dest++ = '\0';
}

// Removes trailing white space
//void RemoveTrailingWhiteSpace(char *dest, char *source)
//{
//}
//*/

