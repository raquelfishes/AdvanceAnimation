
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

	inline float bilerp(const Vec2 pos, Array2<float>& field, const Index2 size){
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
		return (1 - t.y) * y1 + t.y * y2;
	}
}

// advection
void Fluid2::fluidAdvection(const float dt)
{
	// ink advection
	{
		Array2<float> inkAux(ink);
		const Index2 sizeInk = ink.getSize();
		for (unsigned int i = 0; i < sizeInk.x; ++i){
			for (unsigned int j = 0; j < sizeInk.y; ++j){
				const Index2 id(i, j);
				/// Calculate global coordinates
				const Vec2 pos(grid.getCellPos(id));
				/// Calculate advencion method
				const Vec2 vel((velocityX[id] + velocityX[Index2(i+1,j)]) * 0.5f, (velocityY[id] + velocityY[Index2(i,j+1)]) * 0.5f);
				const Vec2 endpos(pos - dt * vel);
				/// Bilineal interpolation with index
				ink[id] = bilerp(grid.getCellIndex(endpos), inkAux, sizeInk);;
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
				/// Calculate index of the previus velocities in X to compute the velocityX
				const Index2 id1(clamp(i - 1, 0, sizeV.x - 1), clamp(j, 0, sizeV.y - 1));
				const Index2 id2(clamp(i, 0, sizeV.x - 1), clamp(j, 0, sizeV.y - 1));
				const Index2 id3(clamp(i - 1, 0, sizeV.x - 1), clamp(j + 1, 0, sizeV.y - 1));
				const Index2 id4(clamp(i, 0, sizeV.x - 1), clamp(j + 1, 0, sizeV.y - 1));
				/// Calculate global coordinates VelocityX
				const Vec2 pos(grid.getFaceXPos(id));
				/// Compute advection methond
				const Vec2 vel(uAux[id], (vAux[id1] + vAux[id2] + vAux[id3] + vAux[id4])*0.25f);
				const Vec2 endpos(pos - dt*vel);
				/// Bilineal interpolation with index, face index in axis 0, horizontal
				velocityX[id] = bilerp(grid.getFaceIndex(endpos, 0), uAux, sizeU);
			}
		}
		for (unsigned int i = 0; i < sizeV.x; ++i){
			for (unsigned int j = 0; j < sizeV.y; ++j){
				const Index2 id(i, j);
				/// Calculate index of the previus velocities in X to compute the velocityX
				const Index2 id1(clamp(i, 0, sizeU.x - 1), clamp(j - 1, 0, sizeU.y - 1));
				const Index2 id2(clamp(i, 0, sizeU.x - 1), clamp(j - 1, 0, sizeU.y - 1));
				const Index2 id3(clamp(i + 1, 0, sizeU.x - 1), clamp(j, 0, sizeU.y - 1));
				const Index2 id4(clamp(i + 1, 0, sizeU.x - 1), clamp(j, 0, sizeU.y - 1));
				/// Calculate global coordinates VelocityY
				const Vec2 pos(grid.getFaceYPos(id));
				/// Compute advection methond
				const Vec2 vel((uAux[id1] + uAux[id2] + uAux[id3] + uAux[id4])*0.25f, vAux[id]);
				const Vec2 endpos(pos - dt*vel);
				/// Bilineal interpolation with index, face index in axis 0, horizontal
				velocityY[id] = bilerp(grid.getFaceIndex(endpos, 1), vAux, sizeV);
			}
		}
	}

}

// emission
void Fluid2::fluidEmission()
{
	if (Scene::testcase >= Scene::SMOKE)
	{
		
		Bbox2 source1(-0.2f, -0.2f, 0.2f, 0.2f);

		Vec2 minCell1 = grid.getCellIndex(source1.minPosition);
		Vec2 maxCell1 = grid.getCellIndex(source1.maxPosition);
		Vec2 difCell1 = (maxCell1 - minCell1) / Vec2(2.0f, 2.0f);
		Vec2 middleCell1 = difCell1 + minCell1;

		Index2 cornerMin1((int)floor(minCell1.x), (int)floor(minCell1.y));
		Index2 cornerMax1((int)ceil(maxCell1.x), (int)ceil(maxCell1.y));
		Index2 middleMin1((int)floor(middleCell1.x), (int)floor(middleCell1.y));
		Index2 middleMax1((int)ceil(middleCell1.x), (int)ceil(middleCell1.y));

		for (unsigned int i = cornerMin1.x; i <= cornerMax1.x; ++i){
			for (unsigned int j = cornerMin1.y; j <= cornerMax1.y; ++j){
				const Index2 id(i, j);
				if (j > cornerMin1.y && j < cornerMax1.y){
					if (i < middleMin1.x){
						ink[id] = 1.0f;
						velocityX[id] = -5.0f;
						velocityY[id] = 0.0f;
					}
					if (i > middleMax1.x){
						ink[id] = 1.0f;
						velocityX[id] = 5.0f;
						velocityY[id] = 0.0f;
					}
				}
			}
		}
		if (rand() > RAND_MAX / 3){
			Bbox2 source2(-0.1f, -1.9f, 0.1f, -1.7f);

			Vec2 minCell2 = grid.getCellIndex(source2.minPosition);
			Vec2 maxCell2 = grid.getCellIndex(source2.maxPosition);

			Index2 cornerMin2((int)floor(minCell2.x), (int)floor(minCell2.y));
			Index2 cornerMax2((int)ceil(maxCell2.x), (int)ceil(maxCell2.y));

			for (unsigned int i = cornerMin2.x; i <= cornerMax2.x; ++i){
				for (unsigned int j = cornerMin2.y; j <= cornerMax2.y; ++j){
					const Index2 id(i, j);
					ink[id] = 1.0f;
					velocityX[id] = 0.0f;
					velocityY[id] = 10.0f;
				}
			}

			Bbox2 source3(-0.1f, 1.7f, 0.1f, 1.9f);

			Vec2 minCell3 = grid.getCellIndex(source3.minPosition);
			Vec2 maxCell3 = grid.getCellIndex(source3.maxPosition);

			Index2 cornerMin3((int)floor(minCell3.x), (int)floor(minCell3.y));
			Index2 cornerMax3((int)ceil(maxCell3.x), (int)ceil(maxCell3.y));

			for (unsigned int i = cornerMin3.x; i <= cornerMax3.x; ++i){
				for (unsigned int j = cornerMin3.y; j <= cornerMax3.y; ++j){
					const Index2 id(i, j);
					ink[id] = 1.0f;
					velocityX[id] = 0.0f;
					velocityY[id] = -10.0f;
				}
			}
		}
	}
}

// volume forces
void Fluid2::fluidVolumeForces(const float dt)
{
	if (Scene::testcase >= Scene::SMOKE)
	{
		// gravity
		const float gravity = dt * Scene::kGravity;
		const Index2 sizeV = velocityY.getSize();
		for (unsigned int i = 0; i < sizeV.x; ++i){
			for (unsigned int j = 0; j < sizeV.y; ++j){
				const Index2 id(i, j);
				velocityY[id] += gravity;
			}
		}
	}
}


// viscosity
void Fluid2::fluidViscosity(const float dt)
{
	if (Scene::testcase >= Scene::SMOKE)
	{
		// viscosity
		Array2< float > uAux(velocityX);
		Array2< float > vAux(velocityY);
		const Index2 sizeU = velocityX.getSize();
		const Index2 sizeV = velocityY.getSize();

		const Vec2 dx(grid.getCellDx());
		const Vec2 invDxPow(1.0f / pow(dx.x, 2), 1.0f / pow(dx.y, 2));
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
				velocityX[id] += dtVisDen * ((uAux[id1] - 2.0f * uAux[id] + uAux[id2])*invDxPow.x + (uAux[id3] - 2 * uAux[id] + uAux[id4])*invDxPow.y);
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
				velocityX[id] += dtVisDen * ((vAux[id1] - 2 * vAux[id] + vAux[id2])*invDxPow.x + (vAux[id3] - 2 * vAux[id] + vAux[id4])*invDxPow.y);
			}
		}

	}
}

// pressure
void Fluid2::fluidPressureProjection(const float dt)
{
	if (Scene::testcase >= Scene::SMOKE)
	{
		// pressure
		const Index2 sizeP = pressure.getSize();
		const Index2 sizeU = velocityX.getSize();
		const Index2 sizeV = velocityY.getSize();

		const Vec2 dx(grid.getCellDx());
		const Vec2 invDx(1.0f / dx.x, 1.0f / dx.y);
		const Vec2 invDxPow(1.0f / pow(dx.x, 2), 1.0f / pow(dx.y, 2));

		/// wall boundary conditions
		for (unsigned int j = 0; j < sizeU.y; ++j){
			/// Left wall as solid
			velocityX[Index2(0, j)] = 0.0f;
			/// Right wall as solid
			velocityX[Index2(sizeU.x - 1, j)] = 0.0f;
		}
		for (unsigned int i = 0; i < sizeV.x; ++i){
			/// Botton as solid
			velocityY[Index2(i, 0)] = 0.0f;
		}

		// b
		const float pDt = Scene::kDensity / dt;
		std::vector< double > b(sizeP.x * sizeP.y);
		for (unsigned int i = 0; i < sizeP.x; ++i){
			for (unsigned int j = 0; j < sizeP.y; ++j){
				const Index2 id(i, j);
				b[pressure.getLinearIndex(i, j)] = -pDt * ((velocityX[Index2(i + 1, j)] - velocityX[id]) * invDx.x + (velocityY[Index2(i, j + 1)] - velocityY[id]) * invDx.y);
			}
		}

		// A
		SparseMatrix< double > A(sizeP.x * sizeP.y, 5);
		for (unsigned int i = 0; i < sizeP.x; ++i){
			for (unsigned int j = 0; j < sizeP.y; ++j){
				const unsigned int id = pressure.getLinearIndex(i, j);
				if (i > 0) {
					const unsigned int id1 = pressure.getLinearIndex(i - 1, j);
					A.add_to_element(id, id, 1. * invDxPow.x);
					A.add_to_element(id, id1, -1. * invDxPow.x);
				}
				if (i < sizeP.x - 1) {
					const unsigned int id1 = pressure.getLinearIndex(i + 1, j);
					A.add_to_element(id, id, 1. * invDxPow.x);
					A.add_to_element(id, id1, -1. * invDxPow.x);
				}
				if (j > 0) {
					const unsigned int id1 = pressure.getLinearIndex(i, j - 1);
					A.add_to_element(id, id, 1. * invDxPow.y);
					A.add_to_element(id, id1, -1. * invDxPow.y);
				}
				/// TOP as air
				A.add_to_element(id, id, 1. * invDxPow.y);
				if (j < sizeP.y - 1) {
					const unsigned int id1 = pressure.getLinearIndex(i, j + 1);
					A.add_to_element(id, id1, -1. * invDxPow.y);
				}
			}
		}

		// pcg solver
		PCGSolver< double > solver;

		double residual_out;
		int iterations_out;
		std::vector< double > p(sizeP.x * sizeP.y);

		solver.set_solver_parameters(1e-6, 10000);
		solver.solve(A, b, p, residual_out, iterations_out);

		std::cout << "Pressure system result: res=" << residual_out << ", iter=" << iterations_out << std::endl;

		// set pressure
		for (unsigned int i = 0; i < sizeP.x; ++i){
			for (unsigned int j = 0; j < sizeP.y; ++j){
				const Index2 id(i, j);
				pressure[id] = (float)p[pressure.getLinearIndex(i, j)];
			}
		}

		// apply pressure gradient
		const float dtOverRho = dt / Scene::kDensity;
		for (unsigned int i = 1; i < sizeU.x - 1; ++i){
			for (unsigned int j = 0; j < sizeU.y; ++j){
				const Index2 id(i, j);
				const Index2 id1(i - 1, j);
				const float gradp = (pressure[id] - pressure[id1]) * invDx.x;
				velocityX[id] -= dtOverRho * gradp;
			}
		}
		for (unsigned int i = 0; i < sizeV.x; ++i){
			for (unsigned int j = 1; j < sizeV.y - 1; ++j){
				const Index2 id(i, j);
				const Index2 id1(i, j - 1);
				const float gradp = (pressure[id] - pressure[id1]) * invDx.y;
				velocityY[id] -= dtOverRho * gradp;
			}
		}
		// apply pressure gradient: TOP as air
		for (unsigned int i = 0; i < sizeV.x; ++i){
			const Index2 id(i, sizeV.y - 1);
			const Index2 id1(i, sizeP.y - 1);
			const float gradp = (0.0f - pressure[id1]) * invDx.y;
			velocityY[id] -= dtOverRho * gradp;
		}
	}
}
