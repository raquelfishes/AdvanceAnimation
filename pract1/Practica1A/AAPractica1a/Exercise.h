void AdvanceTimeStep1(float k, float m, float d, float L, float kc, float dt, int method,
   float p1, float v1, float& p2, float& v2, float& p2old, bool collision);
void AdvanceTimeStep2(float k, float m, float d, float L, float kA, float A, float dt, int method,
   const Vec2& p1, const Vec2& v1, Vec2& p2, Vec2& v2, Vec2& p3, Vec2& v3, Vec2& p4, Vec2& v4, bool springs, bool area);

void eulerExplicit(float k, float m, float d, float L, float kc, float dt, float p1, float v1, float& p2, float& v2, bool collision);
void eulerImplicit(float k, float m, float d, float L, float kc, float dt, float p1, float v1, float& p2, float& v2, bool collision);
void eulerSymplectic(float k, float m, float d, float L, float kc, float dt, float p1, float v1, float& p2, float& v2, bool collision);
void midPoint(float k, float m, float d, float L, float kc, float dt, float p1, float v1, float& p2, float& v2, bool collision);
void verlet(float k, float m, float d, float L, float kc, float dt, float p1, float v1, float& p2, float& v2, float& p2old, bool collision);
void eulerSymplectic2(float k, float m, float d, float L, float kA, float A, float dt, const Vec2& p1, const Vec2& v1, Vec2& p2, Vec2& v2, Vec2& p3, Vec2& v3, Vec2& p4, Vec2& v4, bool springs, bool area);
void eulerImplicit2(float k, float m, float d, float L, float kA, float A, float dt, const Vec2& p1, const Vec2& v1, Vec2& p2, Vec2& v2, Vec2& p3, Vec2& v3, Vec2& p4, Vec2& v4, bool springs, bool area);
