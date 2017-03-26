//@HEADER
// ************************************************************************
//
//               Epetra: Linear Algebra Services Package
//                 Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER


#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif
#include "Epetra_SerialComm.h"
#include "../epetra_test_err.h"
#include "Epetra_Version.h"

// prototypes

int check(Epetra_SerialDenseSolver & solver, double * A1, int LDA,
	  int N1, int NRHS1, double OneNorm1,
	  double * B1, int LDB1,
	  double * X1, int LDX1,
	  bool Transpose, bool verbose);

void GenerateHilbert(double *A, int LDA, int N);

bool Residual( int N, int NRHS, double * A, int LDA, bool Transpose,
	       double * X, int LDX, double * B, int LDB, double * resid);

int matrixCpyCtr(bool verbose, bool debug);
int matrixAssignment(bool verbose, bool debug);
void printHeading(const char* heading);
double* getRandArray(int length);
void printMat(const char* name, Epetra_SerialDenseMatrix& matrix);
void printArray(double* array, int length);
bool identicalSignatures(Epetra_SerialDenseMatrix& a, Epetra_SerialDenseMatrix& b, bool testLDA = true);
bool seperateData(Epetra_SerialDenseMatrix& a, Epetra_SerialDenseMatrix& b);


int main(int argc, char *argv[])
{
  int ierr = 0, i, j, k;
  bool debug = false;
  return ierr ;
}

int check(Epetra_SerialDenseSolver &solver, double * A1, int LDA1,

  return(0);
}

 void GenerateHilbert(double *A, int LDA, int N) {
   for (int j=0; j<N; j++)
     for (int i=0; i<N; i++)
       A[i+j*LDA] = 1.0/((double)(i+j+1));
   return;
 }

bool Residual( int N, int NRHS, double * A, int LDA, bool Transpose,
	       double * X, int LDX, double * B, int LDB, double * resid) {

  Epetra_BLAS Blas;
  char Transa = 'N';
  if (Transpose) Transa = 'T';
  Blas.GEMM(Transa, 'N', N, NRHS, N, -1.0, A, LDA,
	    X, LDX, 1.0, B, LDB);
  bool OK = true;
  for (int i=0; i<NRHS; i++) {
    resid[i] = Blas.NRM2(N, B+i*LDB);
    if (resid[i]>1.0E-7) OK = false;
  }

  return(OK);
}


//=========================================================================
// test matrix operator= (copy & view)

//=========================================================================
// test matrix copy constructor (copy & view)
int matrixCpyCtr(bool verbose, bool debug) {
	const int m1rows = 5;
	const int m1cols = 4;
	const int m2rows = 2;
	const int m2cols = 6;

	int ierr = 0;
	int returnierr = 0;
	if(verbose) printHeading("Testing matrix copy constructors");

	if(verbose) cout << "checking copy constructor (view)" << endl;
	double* m1rand = getRandArray(m1rows * m1cols);
	if(debug) printArray(m1rand, m1rows * m1cols);
	Epetra_SerialDenseMatrix m1(View, m1rand, m1rows, m1rows, m1cols);
	if(debug) {
		cout << "original matrix:" << endl;
		printMat("m1",m1);
	}
	Epetra_SerialDenseMatrix m1clone(m1);
	if(debug) {
		cout << "clone matrix:" << endl;
		printMat("m1clone",m1clone);
	}
	if(verbose) cout << "making sure signatures match" << endl;
	EPETRA_TEST_ERR(!identicalSignatures(m1, m1clone), ierr);
	delete[] m1rand;
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;
	
	if(verbose) cout << "\nchecking copy constructor (copy)" << endl;
	double* m2rand = getRandArray(m2rows * m2cols);
	if(debug) printArray(m2rand, m2rows * m2cols);
	Epetra_SerialDenseMatrix m2(Copy, m2rand, m2rows, m2rows, m2cols);
	if(debug) {
		cout << "original matrix:" << endl;
		printMat("m2",m2);
	}
	Epetra_SerialDenseMatrix m2clone(m2);
	if(debug) {
		cout << "clone matrix:" << endl;
		printMat("m2clone",m2clone);
	}
	if(verbose) cout << "checking that signatures match" << endl;
	EPETRA_TEST_ERR(!identicalSignatures(m2, m2clone), ierr);
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;

	if(verbose) cout << "\nmodifying entry in m2, m2clone should be unchanged" << endl;
	EPETRA_TEST_ERR(!seperateData(m2, m2clone), ierr);
	if(debug) {
		printArray(m2rand, m2rows * m2cols);
		cout << "orig:" << endl;
		printMat("m2",m2);
		cout << "clone:" << endl;
		printMat("m2clone",m2clone);
	}
	delete[] m2rand;
	returnierr += ierr;
	if(ierr == 0)
		if(verbose) cout << "Checked OK." << endl;
	ierr = 0;
	
	return(returnierr);
}

//=========================================================================
// prints section heading with spacers/formatting
void printHeading(const char* heading) {
	cout << "\n==================================================================\n";
	cout << heading << endl;
	cout << "==================================================================\n";
}

//=========================================================================
// prints SerialDenseMatrix/Vector with formatting
void printMat(const char* name, Epetra_SerialDenseMatrix& matrix) {
	//cout << "--------------------" << endl;
	cout << "*** " << name << " ***" << endl;
	cout << matrix;
	//cout << "--------------------" << endl;
}

//=========================================================================
// prints double* array with formatting
void printArray(double* array, int length) {
	cout << "user array (size " << length << "): ";
	for(int i = 0; i < length; i++)
		cout << array[i] << "  ";
	cout << endl;
}

//=========================================================================
// returns a double* array of a given length, with random values on interval (-1,1).
// this is the same generator used in SerialDenseMatrix
double* getRandArray(int length) {
  const double a = 16807.0;
	const double BigInt = 2147483647.0;
	const double DbleOne = 1.0;
	const double DbleTwo = 2.0;
	double seed = rand();

	double* array = new double[length];

	for(int i = 0; i < length; i++) {
		seed = fmod(a * seed, BigInt);
		array[i] = DbleTwo * (seed / BigInt) - DbleOne;
	}

	return(array);
}

//=========================================================================
// checks the signatures of two matrices
bool identicalSignatures(Epetra_SerialDenseMatrix& a, Epetra_SerialDenseMatrix& b, bool testLDA) {

	if((a.M()  != b.M()  )|| // check properties first
		 (a.N()  != b.N()  )||
		 (a.CV() != b.CV() ))
		return(false);

	if(testLDA == true)      // if we are coming from op= c->c #2 (have enough space)
		if(a.LDA() != b.LDA()) // then we don't check LDA (but we do check it in the test function)
			return(false);

	if(a.CV() == View) { // if we're still here, we need to check the data
		if(a.A() != b.A()) // for a view, this just means checking the pointers
			return(false);   // for a copy, this means checking each element
	}
	else { // CV == Copy
		const int m = a.M();
		const int n = a.N();
		for(int i = 0; i < m; i++)
			for(int j = 0; j < n; j++) {
				if(a(i,j) != b(i,j))
					return(false);
			}
	}

	return(true); // if we're still here, signatures are identical
}

//=========================================================================
// checks if two matrices are independent or not
bool seperateData(Epetra_SerialDenseMatrix& a, Epetra_SerialDenseMatrix& b) {
	bool seperate;

	int r = EPETRA_MIN(a.M(),b.M()) / 2; // ensures (r,c) is valid
	int c = EPETRA_MIN(a.N(),b.N()) / 2; // in both matrices

	double orig_a = a(r,c);
	double new_value = a(r,c) + 1;
	if(b(r,c) == new_value) // there's a chance b could be independent, but
		new_value++;          // already have new_value in (r,c).
	
	a(r,c) = new_value;
	if(b(r,c) == new_value)
		seperate = false;
	else
		seperate = true;

	a(r,c) = orig_a; // undo change we made to a

	return(seperate);
}
