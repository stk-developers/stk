/***************************************************************************
 *   copyright            : (C) 2004 by Lukas Burget,UPGM,FIT,VUT,Brno     *
 *   email                : burget@fit.vutbr.cz                            *
 ***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef FILEIO_H
#define FILEIO_H

#include <string>
#include <stdio.h>
#include "common.h"

using namespace std;

typedef struct {
    INT_32 nSamples;              /* Structure for HTK header */
    INT_32 sampPeriod;
    INT_16 sampSize;
    INT_16 sampKind;
}
HTK_Header;


#ifdef __cplusplus
  extern "C" {
#endif


int WriteHTKHeader (FILE * fp_out, HTK_Header header, bool swap);
int WriteHTKFeature (FILE * fp_out, FLOAT *out, size_t fea_len, bool swap);
int ReadHTKHeader (FILE * fp_in, HTK_Header *header, bool swap);

int ReadHTKFeature (FILE * fp_in, FLOAT *in, size_t fea_len, bool swap,
                    bool decompress, FLOAT *A, FLOAT *B);

int Mkdir4File(const char *file_name);

typedef struct {
  char *last_file_name;
  char *last_cmn_file;
  char *last_cvn_file;
  char *last_cvg_file;
  FILE *fp;
  FLOAT *cmn;
  FLOAT *cvn;
  FLOAT *cvg;
  HTK_Header last_header;
  FLOAT *A;
  FLOAT *B;
} RHFBuffer;

FLOAT *ReadHTKFeatures(
  char *file_name,
  bool swap,
  int extLeft,
  int extRight,
  int targetKind,
  int derivOrder,
  int *derivWinLen,
  HTK_Header *header,
  const char *cmn_file,
  const char *cvn_file,
  const char *cvg_file,
  RHFBuffer *buff);

#ifdef __cplusplus
}
#endif




////////////////////////////////////////////////////////////////////////////////
/// File class
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// File class is an interface to common file I/O operation in STK. It gives
/// virtual methods to be later redefined separately for input and output.
////////////////////////////////////////////////////////////////////////////////
class STK_FStream
{
protected:
	/// Holds physical file name
	string fName;

	/// Pointer to the C FILE structure
	FILE* fStream;

	/// Holds info whether the file is open or not
	bool isopen;

public:
	/// A constructor
	STK_FStream();

	/// Abstract method to open the file
	virtual int open(const string& filename) = 0;

	/// Closes the file
	virtual void close();

	/// Returns whether the file is open or not
	bool is_open() const {return isopen;};

	/// @brief Returns a pointer to the FILE* structure
  ///
	/// It is not safe to give direc access to the FILE structure
	/// so const function is given
	FILE* fp() {return fStream;}
};


class mystream : public ifstream
{
public:
  string fName;
};


////////////////////////////////////////////////////////////////////////////////
/// IFile class
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// IFile class encapsulates file input operations.
////////////////////////////////////////////////////////////////////////////////
class STK_IFStream : public STK_FStream
{
public:
	~STK_IFStream();

	/// Opens the file for input
	int open(const string& filename);
};



////////////////////////////////////////////////////////////////////////////////
/// OFile class
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// OFile class encapsulates file output operations.
////////////////////////////////////////////////////////////////////////////////
class STK_OFStream : public STK_FStream
{
public:
	~STK_OFStream();

	/// Opens the file for output
	int open(const string& filename);
};



#endif // FILEIO_H
