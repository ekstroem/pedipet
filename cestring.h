/*
 * PEDIPET
 * 
 * String rutines for plain C
 *
 * cestring.h
 * 
 * Copyright (C) 1996--2004 Claus Ekstrøm
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


#ifndef CESTRING
#define CESTRING

void RemoveWhiteSpace(char *dest, char *source);
// Copies source to dest and removes all white space
// Dest must be large enough to contain source

void SimplifyWhiteSpace(char *dest, char *source);
// Reduces two or more white spaces to a single white space (ASCII 32)
// Dest must be large enough to contain source


#endif
