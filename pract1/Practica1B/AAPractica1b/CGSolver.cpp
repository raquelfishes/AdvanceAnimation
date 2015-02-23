// CGSolver.cpp : Defines the entry point for the console application.
//

#include "CGSolver.h"

CGSolver::CGSolver(void)
{
   nDiag=0;
   nSparse=0;
}

CGSolver::CGSolver(int nodes, int springs)
{
   nDiag=nodes;
   nSparse=springs;
}

CGSolver::~CGSolver(void)
{
   if(!nDiag)
   {
      return;
   }

   delete[] diagA;
   delete[] sparseA;
   delete[] sparserow;
   delete[] sparsecol;
   delete[] b;
   delete[] rsd;
   delete[] dir;
   delete[] Adir;
}

void CGSolver::Init(void)
{
   if(!nDiag)
   {
      return;
   }

   diagA=new Matrix3[nDiag];
   sparseA=new Matrix3[nSparse];
   sparserow=new int[nSparse];
   sparsecol=new int[nSparse];
   b=new Vec3[nDiag];
   rsd=new Vec3[nDiag];
   dir=new Vec3[nDiag];
   Adir=new Vec3[nDiag];
}

//Solve for 'x' in 'A * x = b'
//The solver needs the initial values x0
void CGSolver::SolveCG(Vec3 *x, const Vec3 *x0)
{
   float rrNew=2.0f*tol;
   float rrOld=1.0;

   InitializeCG(x, x0);

   for(int i=0; i<maxIter && rrNew>tol; i++)
   {
      //Compute A * d
      SparseMult(dir, Adir);

      //Compute old r^T * r
      rrOld=InnerProduct(rsd, rsd);

      //Compute d^T * A * d
      float dAd=InnerProduct(dir, Adir);

      //Compute alpha
      float alpha=(fabs(dAd)<1e-32f) ? 0.0f : rrOld/dAd;

      //Compute new iteration value
      VectorSum(x, dir, alpha, x);

      //Recompute residual
      VectorSum(rsd, Adir, -alpha, rsd);

      //Compute new r^T * r
      rrNew=InnerProduct(rsd, rsd);

      //Recompute direction
      VectorSum(rsd, dir, rrNew/rrOld, dir);
   }

}

void CGSolver::InitializeCG(Vec3 *x, const Vec3 *x0)
{
   SparseMult(x0, rsd);
   for(int i=0; i<nDiag; i++)
   {
      dir[i]=rsd[i]=b[i]-rsd[i];
   }
   for(int i=0; i<nDiag; i++)
   {
      x[i]=x0[i];
   }
}

//Product y = M * x
void CGSolver::SparseMult(const Vec3 *x, Vec3 *y) const
{
   int i;

   //Obtain y=M*x via sparse matrix multiplication
   for(i=0; i<nDiag; i++)
   {
      y[i]=diagA[i]*x[i];
   }

   for(i=0; i<nSparse; i++)
   {
      y[sparserow[i]]+=sparseA[i]*x[sparsecol[i]];
      y[sparsecol[i]]+=sparseA[i]^x[sparserow[i]];
   }
}

//return Product x^T * y
float CGSolver::InnerProduct(const Vec3 *x, const Vec3 *y) const
{
   float res=0.0f;
   for(int i=0; i<nDiag; i++)
   {
      res+=x[i].dot(y[i]);
   }
   return res;
}

//Sum r = a + sc*b
void CGSolver::VectorSum(const Vec3 *a, const Vec3 *b, float sc, Vec3 *r) const
{
   for(int i=0; i<nDiag; i++)
   {
      r[i]=a[i]+sc*b[i];
   }
}

void CGSolver::ResetSystem(void)
{
   for(int i=0; i<nDiag; i++)
   {
      diagA[i] = Matrix3::ZERO;
      b[i] = Vec3::ZERO;
   }
   for(int i=0; i<nSparse; i++)
   {
      sparseA[i] = Matrix3::ZERO;
   }
}

void CGSolver::AddToDiagBlock(const Matrix3& val, int i)
{
   diagA[i]+=val;
}

void CGSolver::AddToSparseBlock(const Matrix3& val, int i)
{
   sparseA[i]+=val;
}

void CGSolver::InitSparseBlockIndices(int i, int r, int c)
{
   sparserow[i]=r;
   sparsecol[i]=c;
}

void CGSolver::AddToVectorBlock(const Vec3& val, int i)
{
   b[i]+=val;
}


//Sample Usage
CGSolver *solver;

void Example(void)
{
   int numberOfMovingNodes=10;
   int numberOfSpringsWithOneOrTwoMovingNodes=20; //All springs that are attached to at least 1 moving node
   int numberOfSpringsBetweenTwoMovingNodes=10; //Only the springs that are attached to two moving nodes
   int maxIterations=10; //This is probably too low in terms of accuracy, but you should set it in terms of stability
   float tolerance=1e-10f; //This value is meaningless; ideally one should account for the number of nodes, etc.

   //Initialization
   solver = new CGSolver(numberOfMovingNodes, numberOfSpringsBetweenTwoMovingNodes);
   solver->Init();
   solver->SetMaxIter(maxIterations);
   solver->SetTolerance(tolerance);
      for(int i=0; i<numberOfSpringsWithOneOrTwoMovingNodes; i++)
      {
         int indexOfNodeA=0;
         int indexOfNodeB=1;
         bool nodeAMoves=indexOfNodeA>=0;
         bool nodeBMoves=indexOfNodeB>=0;
         int indexOfSpring=(nodeAMoves && nodeBMoves) ? 0 : -1; //Assign indices only to springs attached to two moving nodes
                                                                //because this index is used for traversing the array of blocks
                                                                //outside the diagonal
         //Initialize node indices
         if(nodeBMoves && nodeAMoves)
         {
            solver->InitSparseBlockIndices(indexOfSpring, indexOfNodeA, indexOfNodeB);
         }
      }


   //Simulation loop
   while(true)
   {
      //Reset linear system A * x = b
      solver->ResetSystem(); //This sets A and b to 0.
      //If you have some constant values, you may want to modify this and do the reset such that the constant part is loaded

      //Collect A matrix and b vector from points
      //The terms of the A matrix probably include mass, etc.
      for(int i=0; i<numberOfMovingNodes; i++)
      {
         int indexOfNode=0; //Assign indices only to moving nodes,
                            //because this index indicates the location of the block in the big matrix
         Matrix3 matrixBlock;
         Vec3 vectorBlock;

         //Add corresponding block to the diagonal of A
         solver->AddToDiagBlock(matrixBlock, indexOfNode);

         //Add corresponding block to b
         solver->AddToVectorBlock(vectorBlock, indexOfNode);
      }

      //Collect A matrix from springs
      //These are typically the stiffness terms
      for(int i=0; i<numberOfSpringsWithOneOrTwoMovingNodes; i++)
      {
         int indexOfNodeA=0;
         int indexOfNodeB=1;
         bool nodeAMoves=indexOfNodeA>=0;
         bool nodeBMoves=indexOfNodeB>=0;
         int indexOfSpring=(nodeAMoves && nodeBMoves) ? 0 : -1; //Assign indices only to springs attached to two moving nodes
                                                                //because this index is used for traversing the array of blocks
                                                                //outside the diagonal

         Matrix3 matrixBlockDiagonalNodeA, matrixBlockDiagonalNodeB;
         Matrix3 matrixBlockNonDiagonal;

         //Add corresponding blocks to the diagonal of A
         if(nodeAMoves)
         {
            solver->AddToDiagBlock(matrixBlockDiagonalNodeA, indexOfNodeA);
         }
         if(nodeBMoves)
         {
            solver->AddToDiagBlock(matrixBlockDiagonalNodeB, indexOfNodeB);
         }
         if(nodeBMoves && nodeAMoves)
         {
            solver->AddToSparseBlock(matrixBlockNonDiagonal, indexOfSpring);
         }
      }

      //You may also want to add other terms, such as external forces

      Vec3 *v0Block=new Vec3[numberOfMovingNodes]; //Initial values for the velocities
                                                                           //(Maybe the velocities from the previous time step)
      Vec3 *vBlock=new Vec3[numberOfMovingNodes]; //Result

      //Solve linear system
      //Solve for the new velocities
      solver->SolveCG(vBlock, v0Block);

   }

}

