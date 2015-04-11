
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

	inline float bilerp(Vec2 pos, Array2<float> field, Index2 size){
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

		const float y1 = (1.0f - t.x) * field[id1] + t.x * field[id2];
		const float y2 = (1.0f - t.x) * field[id3] + t.x * field[id4];
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
				/// Calculate index of the velocitis
				const Index2 idvXaux(clamp(i + 1, 0, velocityX.getSize().x - 1), clamp(j, 0, velocityX.getSize().y - 1));
				const Index2 idvYaux(clamp(i, 0, velocityY.getSize().x - 1), clamp(j+1, 0, velocityY.getSize().y - 1));
				/// Calculate advencion method
				const Vec2 vel((velocityX[id] + velocityX[idvXaux]) * 0.5f, (velocityY[id] + velocityY[idvYaux]) * 0.5f);
				const Vec2 endpos(pos - dt * vel);
				/// Bilineal interpolation with index
				ink[id] = bilerp(grid.getCellIndex(endpos),inkAux,inkAux.getSize());;
			}
		}
    }
	
    // velocity advection
    {
		Array2< float > uAux(velocityX);
		Array2< float > vAux(velocityY);
		const Index2 sizeU = velocityX.getSize();
		const Index2 sizeV = velocityY.getSize();
		for (unsigned int i = 0; i < sizeU.x; ++i){
			for (unsigned int j = 0; j < sizeU.y; ++j){
				const Index2 id(i, j);
				/// Calculate global coordinates VelocityX
				const Vec2 pos(grid.getFaceXPos(id));
				/// Calculate index of the previus velocities in X to compute the velocityX
				const Index2 id1(clamp(i - 1, 0, sizeV.x - 1), clamp(j, 0, sizeV.y - 1));
				const Index2 id2(clamp(i, 0, sizeV.x - 1), clamp(j, 0, sizeV.y - 1));
				const Index2 id3(clamp(i - 1, 0, sizeV.x - 1), clamp(j + 1, 0, sizeV.y - 1));
				const Index2 id4(clamp(i,   0, sizeV.x-1), clamp(j+1, 0, sizeV.y-1));
				/// Compute advection methond
				const Vec2 vel(uAux[id], (vAux[id1] + vAux[id2] + vAux[id3] + vAux[id4])*0.25f);
				const Vec2 endpos(pos - dt*vel);
				/// Bilineal interpolation with index, face index in axis 0, horizontal
				velocityX[id] = bilerp(grid.getFaceIndex(pos, 0), uAux, uAux.getSize());
			}
		}
		for (unsigned int i = 0; i < sizeV.x; ++i){
			for (unsigned int j = 0; j < sizeV.y; ++j){
				const Index2 id(i, j);
				/// Calculate global coordinates VelocityY
				const Vec2 pos(grid.getFaceYPos(id));
				/// Calculate index of the previus velocities in X to compute the velocityX
				const Index2 id1(clamp(i, 0, sizeU.x - 1), clamp(j - 1, 0, sizeU.y - 1));
				const Index2 id2(clamp(i, 0, sizeU.x - 1), clamp(j - 1, 0, sizeU.y - 1));
				const Index2 id3(clamp(i + 1, 0, sizeU.x - 1), clamp(j, 0, sizeU.y - 1));
				const Index2 id4(clamp(i + 1, 0, sizeU.x - 1), clamp(j, 0, sizeU.y - 1));
				/// Compute advection methond
				const Vec2 vel((uAux[id1] + uAux[id2] + uAux[id3] + uAux[id4])*0.25f,vAux[id]);
				const Vec2 endpos(pos - dt*vel);
				/// Bilineal interpolation with index, face index in axis 0, horizontal
				velocityY[id] = bilerp(grid.getFaceIndex(pos, 1), vAux, vAux.getSize());
			}
		}
    }
	
}

// emission
void Fluid2::fluidEmission()
{
	if( Scene::testcase >= Scene::SMOKE )
	{
		Bbox2 source(-0.18f, -1.9f, 0.18f, -1.7f);
		Vec2 ismin = grid.getCellIndex(source.minPosition);
		Vec2 ismax = grid.getCellIndex(source.maxPosition);
		Index2 cornermin((int)floor(ismin.x), (int)floor(ismin.y));
		Index2 cornermax((int)ceil(ismax.x), (int)ceil(ismax.y));
        // emit source ink
        {
			for (unsigned int i = cornermin.x; i <= cornermax.x; ++i)
				for (unsigned int j = cornermin.y; j <= cornermax.y; ++j)
					ink[Index2(i, j)] = 1.0f;
        }
        // emit source velocity
        {
			for (unsigned int i = cornermin.x; i <= cornermax.x + 1; ++i)
				for (unsigned int j = cornermin.y; j <= cornermax.y; ++j)
					velocityX[Index2(i, j)] = 0.0f;

			for (unsigned int i = cornermin.x; i <= cornermax.x; ++i)
				for (unsigned int j = cornermin.y; j <= cornermax.y + 1; ++j)
					velocityY[Index2(i, j)] = 2.0f;
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
		Array2< float > uAux(velocityX);
		Array2< float > vAux(velocityY);
		const Index2 sizeU = velocityX.getSize();
		const Index2 sizeV = velocityY.getSize();

		const Vec2 invDx(1.0f / pow(grid.getCellDx().x, 2), 1.0f / pow(grid.getCellDx().y, 2));
		const float dtVisDen = dt * Scene::kViscosity / Scene::kDensity;

		for (unsigned int i = 0; i < sizeU.x; ++i){
			for (unsigned int j = 0; j < sizeU.y; ++j){
				const Index2 id(i, j);
				/// Calculate index 
				const Index2 id1(clamp(i + 1, 0, sizeU.x - 1), j);
				const Index2 id2(clamp(i - 1, 0, sizeU.x - 1), j);
				const Index2 id3(i, clamp(j + 1, 0, sizeU.y - 1));
				const Index2 id4(i, clamp(j - 1, 0, sizeU.y - 1));
				/// Compute viscosity in U
				velocityX[id] = uAux[id] + dtVisDen * ((uAux[id1] - 2 * uAux[id] + uAux[id2])*invDx.x + (uAux[id3] - 2 * uAux[id] + uAux[id4])*invDx.y);
			}
		}

		for (unsigned int i = 0; i < sizeV.x; ++i){
			for (unsigned int j = 0; j < sizeV.y; ++j){
				const Index2 id(i, j);
				/// Calculate index 
				const Index2 id1(clamp(i + 1, 0, sizeV.x - 1), j);
				const Index2 id2(clamp(i - 1, 0, sizeV.x - 1), j);
				const Index2 id3(i, clamp(j + 1, 0, sizeV.y - 1));
				const Index2 id4(i, clamp(j - 1, 0, sizeV.y - 1));
				/// Compute viscosity in V
				velocityX[id] = vAux[id] + dtVisDen * ((vAux[id1] - 2 * vAux[id] + vAux[id2])*invDx.x + (vAux[id3] - 2 * vAux[id] + vAux[id4])*invDx.y);
			}
		}

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
