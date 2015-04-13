/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
 *  Various functions for importing and exporting data in a variety of formats
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


#include "cestring.h"
#include "pedipet.h"


//---------------
// Export rutines
//---------------

// Solar
int ExportSolarPedigree (individual *indlist);
int ExportSolarMarker (individual *indlist);
int ExportSolarMap (void);
int ExportSolarFrequency (void);
int ExportSolarPhenotype (individual *indlist);
void MakeSolarIBDFile(int marker);

// GeneHunter
void MakeGHdatafile (void);

// CRI-MAP
int ExportCRIMAPGen(individual *indlist, int xlinked);
int ExportCRIMAPPar(individual *indlist);

// Linkage
void importlinkagepedigree (void);
void ReadLinkageParameterFile(void);
void ExportLinkagePedigree (individual *indlist, int aff);
void ExportLinkageParameter (individual *indlist);

// MENDEL
void ExportMendelPedigree (individual *indlist);
void ExportMendelLocus(individual *indlist);
void ExportMendelMap();

// RelPair
int ExportRelpairControl();
int ExportRelpairLocus(int chromosome);
int ExportRelpairPedigree(individual *indlist);

// Merlin
int ExportMerlinPedigree (individual *indlist);
int ExportMerlinDataFile ();
int ExportMerlinMap (void);
int ExportMerlinFreq (void);

// General file io
int FileExists(char *name);

