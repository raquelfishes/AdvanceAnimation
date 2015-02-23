#pragma once

#include "Vec3.h"

class CGSolver
{
protected:
   Matrix3 *diagA;
   Matrix3 *sparseA;
   int *sparserow, *sparsecol;
   Vec3 *b;
   int nDiag, nSparse;
   Vec3 *rsd, *dir, *Adir;
   int maxIter;
   float tol;

public:
   CGSolver(void);
   CGSolver(int nodes, int springs);
   ~CGSolver(void);

   void Init(void);

   void SetTolerance(float t)
   {
      tol=t;
   }

   void SetMaxIter(int m)
   {
      maxIter=m;
   }

   //Set values in 'A' and 'b'
   void ResetSystem(void);

   void AddToDiagBlock(const Matrix3& val, int i);

   void AddToSparseBlock(const Matrix3& val, int i);

   void InitSparseBlockIndices(int i, int r, int c);

   void AddToVectorBlock(const Vec3& val, int i);

   //Solve for 'x' in 'A * x = b'
   //The solver needs the initial values x0
   void SolveCG(Vec3 *x, const Vec3 *x0);

protected:
   void InitializeCG(Vec3 *x, const Vec3 *x0);

   //Product y = M * x
   void SparseMult(const Vec3 *x, Vec3 *y) const;

   //return Product x^T * y
   float InnerProduct(const Vec3 *x, const Vec3 *y) const;

   //Sum r = a + sc*b
   void VectorSum(const Vec3 *a, const Vec3 *b, float sc, Vec3 *r) const;

};

