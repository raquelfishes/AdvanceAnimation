#ifndef VEC3_DEFINED
#define VEC3_DEFINED

#include <math.h>

// simple 3d vector class
//This class includes methods for adding and subtracting 3D vectors,
//dot and cross product, product with a scalar, and computing the length and length squared.
//If you need other operations, you should code them in other methods.
//Do not modify Vec3, because your program will not run in our test applications!
class Vec3
{
protected:
   float v[3];

public:
   Vec3(void){}
   Vec3(float x_, float y_, float z_){v[0]=x_; v[1]=y_; v[2]=z_;}
   ~Vec3(void){}

   const static Vec3 ZERO;

   float operator[](int i) const
   {
      return v[i];
   }
   float& operator[](int i)
   {
      return v[i];
   }

   void operator+=(const Vec3& v_)
   {
      v[0]+=v_[0]; v[1]+=v_[1]; v[2]+=v_[2];
   }
   void operator-=(const Vec3& v_)
   {
      v[0]-=v_[0]; v[1]-=v_[1]; v[2]-=v_[2];
   }

   Vec3 operator+(const Vec3& v_) const
   {
      return Vec3(v[0]+v_[0], v[1]+v_[1], v[2]+v_[2]);
   }
   Vec3 operator-(const Vec3& v_) const
   {
      return Vec3(v[0]-v_[0], v[1]-v_[1], v[2]-v_[2]);
   }

   friend Vec3 operator-(const Vec3& v)
   {
      return Vec3(-v[0], -v[1], -v[2]);
   }
   friend Vec3 operator*(float k, const Vec3& v)
   {
      return Vec3(k*v[0], k*v[1], k*v[2]);
   }

   float dot(const Vec3& v_) const
   {
      return v[0]*v_[0]+v[1]*v_[1]+v[2]*v_[2];
   }
   Vec3 cross(const Vec3& v_) const
   {
      return Vec3(v[1]*v_[2]-v[2]*v_[1], v[2]*v_[0]-v[0]*v_[2], v[0]*v_[1]-v[1]*v_[0]);
   }

   float length(void) const
   {
      return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
   }
   float length_sq(void) const
   {
      return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
   }
   
};

class Matrix3
{
protected:
   float v[3][3];

public:
   Matrix3(void){}
   Matrix3(float a, float b, float c, float d, float e, float f, float g, float h, float i)
   {
      v[0][0]=a; v[0][1]=b; v[0][2]=c;
      v[1][0]=d; v[1][1]=e; v[1][2]=f;
      v[2][0]=g; v[2][1]=h; v[2][2]=i;
   }
   ~Matrix3(void){}

   const static Matrix3 ZERO;

   const static Matrix3 IDENTITY;

   float* operator[](int i)
   {
      return v[i];
   }

   const float* operator[](int i) const
   {
      return v[i];
   }

   void operator+=(const Matrix3& m)
   {
      v[0][0]+=m[0][0]; v[0][1]+=m[0][1]; v[0][2]+=m[0][2]; 
      v[1][0]+=m[1][0]; v[1][1]+=m[1][1]; v[1][2]+=m[1][2];
      v[2][0]+=m[2][0]; v[2][1]+=m[2][1]; v[2][2]+=m[2][2];
   }
   void operator-=(const Matrix3& m)
   {
      v[0][0]-=m[0][0]; v[0][1]-=m[0][1]; v[0][2]-=m[0][2];
      v[1][0]-=m[1][0]; v[1][1]-=m[1][1]; v[1][2]-=m[1][2];
      v[2][0]-=m[2][0]; v[2][1]-=m[2][1]; v[2][2]-=m[2][2];
   }

   Matrix3 operator+(const Matrix3& m) const
   {
      return Matrix3(v[0][0]+m[0][0], v[0][1]+m[0][1], v[0][2]+m[0][2],
         v[1][0]+m[1][0], v[1][1]+m[1][1], v[1][2]+m[1][2],
         v[2][0]+m[2][0], v[2][1]+m[2][1], v[2][2]+m[2][2]);
   }
   Matrix3 operator-(const Matrix3& m) const
   {
      return Matrix3(v[0][0]-m[0][0], v[0][1]-m[0][1], v[0][2]-m[0][2],
         v[1][0]-m[1][0], v[1][1]-m[1][1], v[1][2]-m[1][2],
         v[2][0]-m[2][0], v[2][1]-m[2][1], v[2][2]-m[2][2]);
   }

   friend Matrix3 operator-(const Matrix3& m)
   {
      return Matrix3(-m[0][0], -m[0][1], -m[0][2],
         -m[1][0], -m[1][1], -m[1][2],
         -m[2][0], -m[2][1], -m[2][2]);
   }
   friend Matrix3 operator*(float k, const Matrix3& m)
   {
      return Matrix3(k*m[0][0], k*m[0][1], k*m[0][2],
         k*m[1][0], k*m[1][1], k*m[1][2],
         k*m[2][0], k*m[2][1], k*m[2][2]);
   }

   Matrix3 operator*(const Matrix3& m) const
   {
      Matrix3 res = Matrix3::ZERO;
      for (int i=0; i<3; i++)
      {
         for (int j=0; j<3; j++)
         {
            for (int k=0; k<3; k++)
            {
               res[i][j] += v[i][k] * m[k][j];
            }
         }
      }
      return res;
   }
   Vec3 operator*(const Vec3& m) const
   {
      Vec3 res = Vec3::ZERO;
      for (int i=0; i<3; i++)
      {
         for (int j=0; j<3; j++)
         {
            res[i] += v[i][j] * m[j];
         }
      }
      return res;
   }
   //Transpose multiplication
   Vec3 operator^(const Vec3& m) const
   {
      Vec3 res = Vec3::ZERO;
      for (int i=0; i<3; i++)
      {
         for (int j=0; j<3; j++)
         {
            res[i] += v[j][i] * m[j];
         }
      }
      return res;
   }
   //Outer product
   friend Matrix3 operator^(const Vec3& v1, const Vec3& v2);
};

inline Matrix3 operator^(const Vec3& v1, const Vec3& v2)
{
   return Matrix3(v1[0]*v2[0], v1[0]*v2[1], v1[0]*v2[2],
      v1[1]*v2[0], v1[1]*v2[1], v1[1]*v2[2],
      v1[2]*v2[0], v1[2]*v2[1], v1[2]*v2[2]);
}


#endif // VEC3_DEFINED
