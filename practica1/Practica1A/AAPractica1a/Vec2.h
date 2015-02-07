#ifndef VEC2_DEFINED
#define VEC2_DEFINED

#include <math.h>

// simple 2d vector class
//This class includes methods for adding and subtracting 2D vectors,
//dot and cross product, product with a scalar, and computing the length and length squared.
//If you need other operations, you should code them in other methods.
//Do not modify Vec2, because your program will not run in our test applications!
class Vec2
{
public:
   float x, y;

public:
   Vec2(void){}
   Vec2(float x_, float y_){x=x_; y=y_;}
   ~Vec2(void){}

   const static Vec2 ZERO;

   static Vec2 CrossProd(const Vec2& v)
   {
      return Vec2(-v.y, v.x);
   }

   void operator+=(const Vec2& v)
   {
      x+=v.x; y+=v.y;
   }
   void operator-=(const Vec2& v)
   {
      x-=v.x; y-=v.y;
   }

   Vec2 operator+(const Vec2& v) const
   {
      return Vec2(x+v.x, y+v.y);
   }
   Vec2 operator-(const Vec2& v) const
   {
      return Vec2(x-v.x, y-v.y);
   }

   friend Vec2 operator*(float k, const Vec2& v)
   {
      return Vec2(k*v.x, k*v.y);
   }

   float dot(const Vec2& v) const
   {
      return x*v.x+y*v.y;
   }
   float cross(const Vec2& v) const
   {
      return x*v.y-y*v.x;
   }
   Vec2 rotate90(void) const
   {
      return Vec2(-y, x);
   }

   float length(void) const
   {
      return sqrt(x*x+y*y);
   }
   float length_sq(void) const
   {
      return x*x+y*y;
   }
   
};

class Matrix2
{
public:
   float v[2][2];

   Matrix2(void){}
   Matrix2(float a, float b, float c, float d)
   {
      v[0][0]=a; v[0][1]=b; v[1][0]=c; v[1][1]=d;
   }
   ~Matrix2(void){}

   const static Matrix2 ZERO;

   const static Matrix2 IDENTITY;

   void operator+=(const Matrix2& m)
   {
      v[0][0]+=m.v[0][0]; v[0][1]+=m.v[0][1]; v[1][0]+=m.v[1][0]; v[1][1]+=m.v[1][1];
   }
   void operator-=(const Matrix2& m)
   {
      v[0][0]-=m.v[0][0]; v[0][1]-=m.v[0][1]; v[1][0]-=m.v[1][0]; v[1][1]-=m.v[1][1];
   }

   Matrix2 operator+(const Matrix2& m) const
   {
      return Matrix2(v[0][0]+m.v[0][0], v[0][1]+m.v[0][1], v[1][0]+m.v[1][0], v[1][1]+m.v[1][1]);
   }
   Matrix2 operator-(const Matrix2& m) const
   {
      return Matrix2(v[0][0]-m.v[0][0], v[0][1]-m.v[0][1], v[1][0]-m.v[1][0], v[1][1]-m.v[1][1]);
   }

   friend Matrix2 operator*(float k, const Matrix2& m)
   {
      return Matrix2(k*m.v[0][0], k*m.v[0][1], k*m.v[1][0], k*m.v[1][1]);
   }

   Matrix2 operator*(const Matrix2& m) const
   {
      return Matrix2(v[0][0]*m.v[0][0]+v[0][1]*m.v[1][0],
         v[0][0]*m.v[0][1]+v[0][1]*m.v[1][1],
         v[1][0]*m.v[0][0]+v[1][1]*m.v[1][0],
         v[1][0]*m.v[0][1]+v[1][1]*m.v[1][1]);
   }
   Vec2 operator*(const Vec2& m) const
   {
      return Vec2(v[0][0]*m.x+v[0][1]*m.y,
         v[1][0]*m.x+v[1][1]*m.y);
   }
   friend Matrix2 operator*(const Vec2& v1, const Vec2& v2);
};

inline Matrix2 operator*(const Vec2& v1, const Vec2& v2)
{
   return Matrix2(v1.x*v2.x, v1.x*v2.y, v1.y*v2.x, v1.y*v2.y);
}


#endif // VEC2_DEFINED
