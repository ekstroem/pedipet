/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
 * Generic list rutines
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


#ifndef _LISTS_H_
#define _LISTS_H_

// Simple list structure
typedef struct list
{
	struct list *next;
} list;


// Adds a single element p1 to the list l1
void addlist (void *l1, void *p1);

// Adds a list of elements (starting at) p1 to the list l1
void addtolist (void *l1, void *p1);


// Removes an element and frees its memory
void removelist (void *l,void *p);

// Removes an element from a list, but does not free memory
void choplist (list **l,list *p);

// Removes element p from l1 and moves it to list l2
void translist (void *l1,void *l2,void *p);

// Removes an entire list
void freelist (void *p1);

// Returns the number of elements in the list
int listlen (void *l);

#endif
