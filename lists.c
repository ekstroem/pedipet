/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
 * List rutines
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

#include <assert.h>
#include <stdlib.h>

#include "lists.h"

void addlist (void *l1,void *p1)
{
	list **l;
	list *p,*q;

	l = (list **) l1;
	p = (list *) p1;

	p->next = 0;

	if (*l)
	{
		for (q = *l; q->next; q = q->next);
		q->next = p;
	}
	else
		*l = p;
}


void addtolist (void *l1,void *p1)
{
	list **l;
	list *p,*q;

	l = (list **) l1;
	p = (list *) p1;

	if (*l)
	{
		for (q = *l; q->next; q = q->next);
		q->next = p;
	}
	else
		*l = p;
}


void choplist (list **l,list *p)
{
	list *q;

	if (*l == p)
		*l = p->next;
	else
	{
		for (q = *l; q->next != p; q = q->next)
			assert (q);
		q->next = p->next;
	}
}

void translist (void *l1,void *l2,void *p)
{
	choplist (l1,p);
	addlist (l2,p);
}

void removelist (void *l,void *p)
{
	choplist (l,p);
	free (p);
}

void freelist (void *p1)
{
	list *p,*p2;

	p = p1;

	while (p) {
	  p2 = p->next;
	  free (p);
	  p = p2;
	}
	// Modified below from p = 0 (????)
        p1 = NULL;
}

int listlen (void *l)
{
	int i=0;
	list *p;

	for (p = l, i = 0; p; p = p->next, i++)
	  ; // Do nothing
	return i;
}
