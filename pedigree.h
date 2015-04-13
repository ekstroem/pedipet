/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
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



#ifndef _PEDIGREE_H_
#ifndef _PEDIPET_H_
#define _PEDIGREE_H_

#include "defs.h"


static double sqrarg;
#define	 addlist2(l,p)	       (*l = p, l = &p->next)
#define	 min(a,b)              ((a) < (b) ? (a) : (b))
#define	 max(a,b)              ((a) > (b) ? (a) : (b))
#define  absolut(a)            ((a) < 0.0 ? (-a) : (a))
#define	 xisdigit(c)	       ((c) == '-' || ((c) >= '0' && (c) <= '9'))
#define	 addptr(p,i)	       ((void *)(((char *)p) + i))
#define  SQR(a)                ((sqrarg=(a))==0.0 ? 0.0 : sqrarg*sqrarg)
#define  strcmpl(s1,s2)        (strcasecmp(s1,s2))
#define  forind                for (ind = individuals; ind; ind = ind->next)

#define  MISSING               -1

#define  OPTIONSFILE           "options.pp"
#define  ERRORFILE             "errors.pp"
#define  DESCFILE              "desc.pp"
#define  HEFILE                "haselst.s"
#define  SYSMISS               "."       // Symbol for the system missing value
#define  MATCHKEYWORD          3  // Number of letters used for command mathcing
#include "enums.h"
#include "structure.h"


#endif
#endif

