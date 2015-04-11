
#include "scene.h"
#include "pcg_solver.h"

namespace
{
	//////////////////////////////////////////////
	// Add any custom classes or functions here! //
	//////////////////////////////////////////////

	inline int clamp(int x, int a, int b){
		return (x<a) ? a : ((x>b) ? b : x);
	}

	inline float bilerp(Vec2 pos, Array2<float> ink, Index2 size){
		const Vec2 posMin(floorf(pos.x), floorf(pos.y));
		const Vec2 posMax(ceilf(pos.x), ceilf(pos.y));
		const Vec2 t(pos - posMin);

		const Index2 id1(clamp((int)posMin.x, 0, size.x - 1), clamp((int)posMin.y, 0, size.y - 1));
		const Index2 id2(clamp((int)posMax.x, 0, size.x - 1), clamp((int)posMin.y, 0, size.y - 1));
		const Index2 id3(clamp((int)posMin.x, 0, size.x - 1), clamp((int)posMax.y, 0, size.y - 1));
		const Index2 id4(clamp((int)posMax.x, 0, size.x - 1), clamp((int)posMax.y, 0, size.y - 1));

		///// Bilineal interpolation - General case
		/// alfa = x - i
		/// beta = y - j
		/// U_alfa_j =  (1 - alfa) * u_ij + alfa * u_i+1j
		/// U_alfa_j+1 = (1 - alfa) * u_ij+1 + alfa * u_i+1j+1
		/// I = (1 - beta)u_alfa_j + beta * u_alfa_j+1

		const float y1 = (1.0f - t.x) * ink[id1] + t.x * ink[id2];
		const float y2 = (1.0f - t.x) * ink[id3] + t.x * ink[id4];
		const float value = (1 - t.y) * y1 + t.y * y2;
		return value;

	}

    
}

// advection
void Fluid2::fluidAdvection( const float dt )
{
    // ink advection
    {
		Array2<float> inkAux(ink);
		for (unsigned int i = 0; i < ink.getSize().x; ++i){
			for (unsigned int j = 0; j < ink.getSize().y; ++j){
				const Index2 id(i, j);

				/// Calculate global coordinates
				const Vec2 pos(grid.getCellPos(id));
				/// Calculate advencion method
				const Vec2 vel((velocityX[id] + velocityX[Index2(i + 1, j)]) * 0.5f, (velocityY[id] + velocityY[Index2(i, j + 1)]) * 0.5f);
				const Vec2 endpos(pos - dt * vel);
				/// Bilineal interpolation
				ink[id] = bilerp(grid.getCellIndex(endpos),inkAux,inkAux.getSize());;
			}
		}
    }

    // velocity advection
    {
		Array2< float > vXaux(velocityX);
		Array2< float > vYaux(velocityY);
		for (unsigned int i = 0; i < velocityX.getSize().x; ++i){
			for (unsigned int j = 0; j < velocityY.getSize().y; ++j){

			}
		}
    }
}

// emission
void Fluid2::fluidEmission()
{
	if( Scene::testcase >= Scene::SMOKE )
	{
        // emit source ink
        {
        }
        // emit source velocity
        {
        }
	}
}

// volume forces
void Fluid2::fluidVolumeForces( const float dt )
{
	if( Scene::testcase >= Scene::SMOKE )
	{
        // gravity
	}
}

// viscosity
void Fluid2::fluidViscosity( const float dt )
{
	if( Scene::testcase >= Scene::SMOKE )
	{
        // viscosity
	}
}

// pressure
void Fluid2::fluidPressureProjection( const float dt )
{
	if( Scene::testcase >= Scene::SMOKE )
	{
        // pressure
	}
}
