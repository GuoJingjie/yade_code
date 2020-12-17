/*************************************************************************
*  Copyright (C) 2014 by Bruno Chareyre <bruno.chareyre@grenoble-inp.fr>     *
*  Copyright (C) 2013 by T. Sweijen (T.sweijen@uu.nl)                    *
*  Copyright (C) 2012 by Chao Yuan <chao.yuan@3sr-grenoble.fr>           *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifdef TWOPHASEFLOW
#include "TwoPhaseFlowEngine.hpp"
#include <boost/range/algorithm_ext/erase.hpp>

namespace yade { // Cannot have #include directive inside.

using math::max;
using math::min; // using inside .cpp file is ok.

#define INFT(cell0) solver->T[solver->currentTes].Triangulation().is_infinite(cell0)

void TwoPhaseFlowEngine::verifyCompatibilityBC()
{
	//This is merely a function to Real check boundary conditions to avoid ill-posed B.C.
	std::cerr << endl << "Boundary and initial conditions are set for: ";

	if (drainageFirst && primaryTPF) {
		std::cerr << "Primary Drainage";
		if (initialPC > -1 * waterBoundaryPressure) {
			std::cerr << endl << "Warning, initial capillary pressure larger than imposed capillary pressure, this may cause imbibition";
		}
	}
	if (drainageFirst && !primaryTPF) {
		std::cerr << "Secondary Drainage";
		if (initialPC > -1 * waterBoundaryPressure) {
			std::cerr << endl << "Warning, initial capillary pressure larger than imposed capillary pressure, this may cause imbibition";
		}
	}
	if (!drainageFirst && primaryTPF) {
		std::cerr << "Primary Imbibition";
		if (initialPC < -1 * waterBoundaryPressure) {
			std::cerr << endl << "Warning, initial capillary pressure smaller than imposed capillary pressure, this may cause drainage";
		}
	}
	if (!drainageFirst && !primaryTPF) {
		std::cerr << "Secondary Imbibition";
		if (initialPC < -1 * waterBoundaryPressure) {
			std::cerr << endl << "Warning, initial capillary pressure smaller than imposed capillary pressure, this may cause drainage";
		}
	}

	std::cout << endl << "Water pressure at: " << waterBoundaryPressure << " and air pressure at: " << airBoundaryPressure << " InitialPC: " << initialPC;
}

void TwoPhaseFlowEngine::setPoreNetwork()
{
	//Reorder cell id's
	RTriangulation&     tri     = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	unsigned int        i       = 0;
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (!cell->info().isFictious) {
			if (cell->info().poreId == -1) {
				cell->info().poreId = i;
				if (cell->info().mergedID > 0) {
					for (FiniteCellsIterator Mcell = tri.finite_cells_begin(); Mcell != cellEnd; Mcell++) {
						if (Mcell->info().mergedID == cell->info().mergedID) {
							Mcell->info().poreId = i;
						}
					}
				}
				i = i + 1;
			}
		}
	}
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (!cell->info().isFictious) {
			if (cell->info().poreId == -1) {
				std::cout << " cell -1 " << cell->info().id;
			}
		}
	}
	numberOfPores = i;
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (!cell->info().isFictious) {
			for (unsigned int k = 0; k < 4; k++) {
				if (!cell->neighbor(k)->info().isFictious) {
					if (cell->info().mergedID == 0
					    || (cell->neighbor(k)->info().mergedID != cell->info().mergedID && cell->info().mergedID != 0)) {
						cell->info().poreIdConnectivity[k] = cell->neighbor(k)->info().poreId; //FIXME: REDUNDANT
					} else
						cell->info().poreIdConnectivity[k] = -1;
				}
			}
		}
	}


	makeListOfPoresInCells(false);
}


void TwoPhaseFlowEngine::setListOfPores()
{
	RTriangulation&     tri     = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	bool                stop    = false;

	//set list of pores
	if ((deformation && remesh) || firstDynTPF) {
		listOfPores.clear();
		for (unsigned int j = 0; j < numberOfPores; j++) {
			for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
				if (cell->info().poreId == (int)j && !cell->info().isGhost && !stop) {
					listOfPores.push_back(cell);
					stop = true;
				}
			}
			stop = false;
		}


		for (unsigned int i = 0; i < numberOfPores; i++) {
			for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
				if (cell->info().poreId == (int)i && !listOfPores[i]->info().isNWRes) {
					for (unsigned int j = 0; j < 4; j++) {
						if (cell->neighbor(j)->info().isWRes) {
							listOfPores[i]->info().isWResInternal   = true;
							listOfPores[i]->info().conductivityWRes = cell->info().kNorm()[j];
						}
					}
				}
			}
		}

		for (unsigned int i = 0; i < numberOfPores; i++) {
			if (listOfPores[i]->info().isWResInternal) { //Important to track for centroidAverage water pressure
				for (unsigned int j = 0; j < listOfPores[i]->info().poreNeighbors.size(); j++) {
					if (!listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().isWResInternal) {
						listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().waterBC = true;
					}
				}
			}


			for (unsigned int j = 0; j < listOfPores[i]->info().poreNeighbors.size(); j++) {
				for (unsigned int k = 0; k < listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().poreNeighbors.size(); k++) {
					if (listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().poreNeighbors[k] == (int)i) {
						listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().listOfEntryPressure[k]
						        = listOfPores[i]->info().listOfEntryPressure[j];
						listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().listOfThroatArea[k]
						        = listOfPores[i]->info().listOfThroatArea[j];
						listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().listOfkNorm[k]
						        = listOfPores[i]->info().listOfkNorm[j];
						listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().listOfkNorm2[k]
						        = listOfPores[i]->info().listOfkNorm2[j];
					}
				}
			}
		}


		for (unsigned int i = 0; i < numberOfPores; i++) {
			listOfPores[i]->info().minSaturation = 1e-3;
		}
	}
	//update for deformation


	//reset knorm
	for (unsigned int i = 0; i < numberOfPores; i++) {
		listOfPores[i]->info().listOfkNorm.clear();

		listOfPores[i]->info().listOfkNorm = listOfPores[i]->info().listOfkNorm2;

		if (listOfPores[i]->info().saturation < 1.0) {
			for (unsigned int j = 0; j < listOfPores[i]->info().poreNeighbors.size(); j++) {
				if (listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().saturation < 1.0) {
					if (-1.0 * listOfPores[i]->info().p() > listOfPores[i]->info().listOfEntryPressure[j]) {
						Real radiusCurvature = 0.0;
						if (listOfPores[i]->info().saturation < 1.0
						    || listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().saturation == 1.0) {
							radiusCurvature = 2.0 * surfaceTension / (-1.0 * listOfPores[i]->info().p());
						}
						if (listOfPores[i]->info().saturation == 1.0
						    || listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().saturation < 1.0) {
							radiusCurvature = 2.0 * surfaceTension
							        / (-1.0 * listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().p());
						}
						if (listOfPores[i]->info().saturation < 1.0
						    || listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().saturation < 1.0) {
							radiusCurvature = 4.0 * surfaceTension
							        / (-1.0
							           * (listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().p()
							              + listOfPores[i]->info().p()));
						}

						Real areaWater = listOfPores[i]->info().listOfThroatArea[j] - 3.1415926535 * radiusCurvature * radiusCurvature;
						if (areaWater < 0.0) {
							areaWater = listOfPores[i]->info().listOfThroatArea[j];
						}

						Real hydraulicRad = 4.0 * surfaceTension
						        / (-1.0
						           * (listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().p() + listOfPores[i]->info().p()));

						Point&  p1                            = listOfPores[i]->info();
						Point&  p2                            = listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info();
						CVector l                             = p1 - p2;
						Real    distance                      = sqrt(l.squared_length());
						listOfPores[i]->info().listOfkNorm[j] = areaWater * hydraulicRad * hydraulicRad / (viscosity * distance);
						if (listOfPores[i]->info().listOfkNorm[j] < 1e-15) {
							listOfPores[i]->info().listOfkNorm[j] = 1e-15;
						}
						if (listOfPores[i]->info().listOfkNorm[j] < 0.0
						    || listOfPores[i]->info().listOfkNorm[j] != listOfPores[i]->info().listOfkNorm[j]) {
							std::cerr << " Error! " << listOfPores[i]->info().listOfkNorm[j];
						}


						for (unsigned int k = 0; k < listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().poreNeighbors.size();
						     k++) {
							if (listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().poreNeighbors[k] == (int)i) {
								listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().listOfkNorm[k]
								        = listOfPores[i]->info().listOfkNorm[j];
							}
						}
					}
				}
			}
		}
	}

	for (unsigned int i = 0; i < numberOfPores; i++) {
		for (unsigned int j = 0; j < listOfPores[i]->info().poreNeighbors.size(); j++) {
			for (unsigned int k = 0; k < listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().poreNeighbors.size(); k++) {
				if (listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().poreNeighbors[k] == (int)i) {
					listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().listOfkNorm[k] = listOfPores[i]->info().listOfkNorm[j];
				}
			}
		}
	}
}


void TwoPhaseFlowEngine::solvePressure()
{
	RTriangulation&     tri     = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	Real                oldDT   = 0.0;

	//Define matrix, triplet list, and linear solver
	tripletList.clear(); // tripletList.resize(T_nnz);
	VectorXr residualsList(numberOfPores);
	VectorXr pressuresList(numberOfPores); //Solve aMatrix * pressuresList = residualsList

	//define lists
	if ((deformation && remesh) || firstDynTPF) {
		aMatrix.resize(numberOfPores, numberOfPores);
		saturationList.assign(numberOfPores, 0.0);
		hasInterfaceList.assign(numberOfPores, false);
		listOfFlux.assign(numberOfPores, 0.0);
		listOfMergedVolume.assign(numberOfPores, 0.0); //NOTE CHANGED AFTER PUSH ON GIT
	}

	//reset various lists
	for (unsigned int i = 0; i < numberOfPores; i++) {
		residualsList[i]    = 0.0;
		pressuresList[i]    = 0.0;
		saturationList[i]   = listOfPores[i]->info().saturation;
		hasInterfaceList[i] = listOfPores[i]->info().hasInterface;
		listOfFlux[i]       = 0.0;
	}

	//Fill matrix
	for (unsigned int i = 0; i < numberOfPores; i++) {
		//Get diagonal coeff
		Real dsdp2  = 0.0;
		Real coeffA = 0.0, coeffA2 = 0.0;
		if (hasInterfaceList[i] && !firstDynTPF) {
			if (listOfPores[i]->info().p() == 0) {
				std::cout << endl << "Error, pressure = 0 " << listOfPores[i]->info().p() << listOfPores[i]->info().id;
				listOfPores[i]->info().p() = -1.0 * listOfPores[i]->info().thresholdPressure;
			}
			dsdp2  = dsdp(listOfPores[i], listOfPores[i]->info().p());
			coeffA = dsdp2
			        * ((listOfPores[i]->info().mergedVolume / scene->dt)
			           + (listOfPores[i]->info().accumulativeDV
			              - listOfPores[i]->info().accumulativeDVSwelling)); //Only consider the change in porosity due to particle movement
		}


		//fill matrix off-diagonals
		if (!listOfPores[i]->info().isWResInternal) {
			for (unsigned int j = 0; j < listOfPores[i]->info().poreNeighbors.size(); j++) {
				tripletList.push_back(ETriplet(i, listOfPores[i]->info().poreNeighbors[j], -1.0 * listOfPores[i]->info().listOfkNorm[j]));
				coeffA2 += listOfPores[i]->info().listOfkNorm[j];
			}
		}

		//Set boundary conditions
		if (listOfPores[i]->info().isWResInternal) {
			tripletList.push_back(ETriplet(i, i, 1.0));
			residualsList[i] = waterBoundaryPressure;
		}

		//Fill matrix diagonal
		if (!listOfPores[i]->info().isWResInternal) {
			if (hasInterfaceList[i]) {
				residualsList[i] += -1.0 * listOfPores[i]->info().saturation
				                * (listOfPores[i]->info().accumulativeDV - listOfPores[i]->info().accumulativeDVSwelling)
				        + coeffA * listOfPores[i]->info().p();
			}
			if (!hasInterfaceList[i] && deformation && listOfPores[i]->info().saturation > listOfPores[i]->info().minSaturation) {
				residualsList[i] += -1.0 * (listOfPores[i]->info().accumulativeDV - listOfPores[i]->info().accumulativeDVSwelling);
			}
			tripletList.push_back(ETriplet(i, i, coeffA + coeffA2));
		}
	}


	//Solve Matrix
	aMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
	eSolver.analyzePattern(aMatrix);
	eSolver.factorize(aMatrix);
	eSolver.compute(aMatrix);


	//Solve for pressure: FIXME: add check for quality of matrix, if problematic, skip all below.
	pressuresList = eSolver.solve(residualsList);

	//Compute flux
	Real flux                = 0.0;
	Real accumulativeDefFlux = 0.0;
	Real waterBefore = 0.0, waterAfter = 0.0;
	Real boundaryFlux = 0.0, lostVolume = 0.0;
	oldDT = scene->dt;

	//check water balance
	for (unsigned int i = 0; i < numberOfPores; i++) {
		waterBefore += listOfPores[i]->info().saturation * listOfPores[i]->info().mergedVolume;
	}

	//compute flux
	for (unsigned int i = 0; i < numberOfPores; i++) {
		flux = 0.0;
		for (unsigned int j = 0; j < listOfPores[i]->info().poreNeighbors.size(); j++) {
			flux += listOfPores[i]->info().listOfkNorm[j] * (pressuresList[i] - pressuresList[listOfPores[i]->info().poreNeighbors[j]]);
		}
		listOfFlux[i] = flux;
		if (listOfPores[i]->info().saturation > listOfPores[i]->info().minSaturation) {
			accumulativeDefFlux += listOfPores[i]->info().accumulativeDV;
		}
		if (listOfPores[i]->info().isWResInternal) {
			boundaryFlux += flux;
		}
		if (!listOfPores[i]->info().isWResInternal && !hasInterfaceList[i] && math::abs(listOfFlux[i]) > 1e-15 && !deformation) {
			std::cerr << " | Flux not 0.0" << listOfFlux[i] << " isNWRES:  " << listOfPores[i]->info().isNWRes
			          << " saturation: " << listOfPores[i]->info().saturation << " P:" << listOfPores[i]->info().p()
			          << " isNWef:" << listOfPores[i]->info().isNWResDef << "|";
			lostVolume += listOfFlux[i] * scene->dt;
		}
	}


	Real summFluxList  = 0.0;
	Real summFluxUnsat = 0.0;
	for (unsigned int i = 0; i < numberOfPores; i++) {
		if (hasInterfaceList[i]) {
			summFluxUnsat += listOfFlux[i];
		}
		summFluxList += listOfFlux[i];
	}

	//update saturation
	for (unsigned int i = 0; i < numberOfPores; i++) {
		if (!deformation && hasInterfaceList[i] && listOfFlux[i] != 0.0 && !listOfPores[i]->info().isWResInternal) {
			Real ds = -1.0 * scene->dt * (listOfFlux[i])
			        / (listOfPores[i]->info().mergedVolume + listOfPores[i]->info().accumulativeDV * scene->dt);
			saturationList[i] = ds
			        + (saturationList[i] * listOfPores[i]->info().mergedVolume
			           / (listOfPores[i]->info().mergedVolume + listOfPores[i]->info().accumulativeDV * scene->dt));
		}
		if (deformation && hasInterfaceList[i] && !listOfPores[i]->info().isWResInternal) {
			saturationList[i] = saturationList[i]
			        + (pressuresList[i] - listOfPores[i]->info().p()) * dsdp(listOfPores[i], listOfPores[i]->info().p())
			        + (scene->dt
			           / (listOfPores[i]->info().mergedVolume
			              + (listOfPores[i]->info().accumulativeDV - listOfPores[i]->info().accumulativeDVSwelling) * scene->dt))
			                * listOfPores[i]->info().accumulativeDVSwelling;
		}
	}

	for (unsigned int i = 0; i < numberOfPores; i++) {
		waterAfter += saturationList[i] * (listOfPores[i]->info().mergedVolume + listOfPores[i]->info().accumulativeDV * scene->dt);
	}


	accumulativeFlux += (summFluxList)*scene->dt;
	accumulativeDeformationFlux += accumulativeDefFlux * scene->dt;

	fluxInViaWBC += boundaryFlux * scene->dt;


	if (!deformation && math::abs(boundaryFlux * scene->dt + (waterBefore - waterAfter)) / math::abs(boundaryFlux * scene->dt) > 1e-3
	    && math::abs(boundaryFlux) > 1e-18) { //FIXME test has to optimized for deforming pore units
		std::cerr << endl
		          << "No volume balance! Flux balance: WBFlux:"
		          << math::abs(boundaryFlux * scene->dt + (waterBefore - waterAfter)) / math::abs(boundaryFlux * scene->dt) << " "
		          << boundaryFlux * scene->dt << "Flux: " << summFluxList * scene->dt << "deltaVolume: " << waterBefore - waterAfter
		          << "Flux in IFACE: " << summFluxUnsat * scene->dt << " lostVolume: " << lostVolume * scene->dt;
		// 	  stopSimulation = true;
	}
	// --------------------------------------find new dt -----------------------------------------------------
	Real dt = 0.0, finalDT = 1e6;
	int  saveID = -1;
	for (unsigned int i = 0; i < numberOfPores; i++) {
		//Time step for deforming pore units
		if (deformation) {
			dt = -1.0 * listOfPores[i]->info().mergedVolume
			        / (listOfPores[i]->info().accumulativeDV + listOfPores[i]->info().accumulativeDVSwelling); //Residence time total pore volume
			if (dt > deltaTimeTruncation && dt < finalDT) {
				finalDT = dt;
				saveID  = -1;
			}

			if (listOfPores[i]->info().accumulativeDVSwelling > 0.0
			    || listOfPores[i]->info().accumulativeDV > 0.0) { // Residence time during increase in pore size
				if (listOfPores[i]->info().accumulativeDVSwelling > listOfPores[i]->info().accumulativeDV) {
					dt = listOfPores[i]->info().mergedVolume * (1.0 - saturationList[i]) / listOfPores[i]->info().accumulativeDVSwelling;
				}
				if (listOfPores[i]->info().accumulativeDVSwelling <= listOfPores[i]->info().accumulativeDV) {
					dt = listOfPores[i]->info().mergedVolume * (1.0 - saturationList[i]) / listOfPores[i]->info().accumulativeDV;
				}
				if (dt > deltaTimeTruncation && dt < finalDT) {
					finalDT = dt;
					saveID  = -2;
				}
			}
			if (listOfPores[i]->info().accumulativeDVSwelling < 0.0 || listOfPores[i]->info().accumulativeDV < 0.0) {
				if (listOfPores[i]->info().accumulativeDVSwelling < listOfPores[i]->info().accumulativeDV) {
					dt = -1.0 * listOfPores[i]->info().mergedVolume * saturationList[i] / listOfPores[i]->info().accumulativeDVSwelling;
				}
				if (listOfPores[i]->info().accumulativeDVSwelling >= listOfPores[i]->info().accumulativeDV) {
					dt = -1.0 * listOfPores[i]->info().mergedVolume * saturationList[i] / listOfPores[i]->info().accumulativeDV;
				}
				if (dt > deltaTimeTruncation && dt < finalDT) {
					finalDT = dt;
					saveID  = -2;
				}
			}
		}

		//Time step for dynamic flow
		if (hasInterfaceList[i]) {
			//thresholdSaturation
			if (math::abs(listOfPores[i]->info().thresholdSaturation - saturationList[i]) > truncationPrecision) {
				dt = -1.0 * (listOfPores[i]->info().thresholdSaturation - saturationList[i]) * listOfPores[i]->info().mergedVolume
				        / listOfFlux[i];
				if (dt > deltaTimeTruncation && dt < finalDT) {
					finalDT = dt; /*saveID = 1;*/
				}
			}
			//Empty pore
			if (math::abs(0.0 - saturationList[i]) > truncationPrecision && listOfFlux[i] > 0.0) { //only for drainage
				dt = -1.0 * (0.0 - saturationList[i]) * listOfPores[i]->info().mergedVolume / listOfFlux[i];
				if (dt > deltaTimeTruncation && dt < finalDT) {
					finalDT = dt; /*saveID = 2;*/
				}
			}
			//Saturated pore
			if (math::abs(1.0 - saturationList[i]) > truncationPrecision && listOfFlux[i] < 0.0) { //only for imbibition
				dt = -1.0 * (1.0 - saturationList[i]) * listOfPores[i]->info().mergedVolume / listOfFlux[i];
				if (dt > deltaTimeTruncation && dt < finalDT) {
					finalDT = dt; /*saveID = 3;*/
				}
			}
		}
	}
	if (finalDT == 1e6) {
		finalDT = deltaTimeTruncation;
		saveID  = 5;
		if (!firstDynTPF && !remesh) {
			std::cout << endl << "NO dt found!";
			stopSimulation = true;
		}
	}
	scene->dt = finalDT * safetyFactorTimeStep;
	if (debugTPF) {
		std::cerr << endl << "Time step: " << finalDT << " Limiting process:" << saveID;
	}


	// --------------------------------------update cappilary pressure (need to correct for linearization of ds/dp)-----------------------------------------------------
	for (unsigned int i = 0; i < numberOfPores; i++) {
		if (hasInterfaceList[i] && !listOfPores[i]->info().isWResInternal && (!deformation || listOfPores[i]->info().saturation != 0.0)) {
			pressuresList[i] = porePressureFromPcS(listOfPores[i], saturationList[i]);
		}
	}


	// --------------------------------------Find invasion events-----------------------------------------------------
	for (unsigned int i = 0; i < numberOfPores; i++) {
		if (saturationList[i] > 1.0 - truncationPrecision && (listOfFlux[i] < 0.0 || deformation) && saturationList[i] != 1.0) {
			if (saturationList[i] > 1.0) {
				if (saturationList[listOfPores[i]->info().invadedFrom] >= 1.0) {
					waterVolumeTruncatedLost += (saturationList[i] - 1.0) * listOfPores[i]->info().mergedVolume;
					saturationList[i] = 1.0;
				}
				if (saturationList[listOfPores[i]->info().invadedFrom] < 1.0) {
					saturationList[listOfPores[i]->info().invadedFrom] += (saturationList[i] - 1.0) * listOfPores[i]->info().mergedVolume
					        / listOfPores[listOfPores[i]->info().invadedFrom]->info().mergedVolume;
					saturationList[i] = 1.0;
					if (saturationList[listOfPores[i]->info().invadedFrom] > 1.0) {
						waterVolumeTruncatedLost += (saturationList[listOfPores[i]->info().invadedFrom] - 1.0)
						        * listOfPores[listOfPores[i]->info().invadedFrom]->info().mergedVolume;
						saturationList[listOfPores[i]->info().invadedFrom] = 1.0;
					}
				}
			}
			saturationList[i]                 = 1.0;
			hasInterfaceList[i]               = false;
			listOfPores[i]->info().isNWResDef = true;
		}
	}

	//Check for drainage
	for (unsigned int i = 0; i < numberOfPores; i++) {
		if ((fractionMinSaturationInvasion == -1 && hasInterfaceList[i] && saturationList[i] < listOfPores[i]->info().thresholdSaturation)
		    || (fractionMinSaturationInvasion > 0.0 && saturationList[i] < fractionMinSaturationInvasion)) {
			for (unsigned int j = 0; j < listOfPores[i]->info().poreNeighbors.size(); j++) {
				if (airBoundaryPressure - pressuresList[listOfPores[i]->info().poreNeighbors[j]] > listOfPores[i]->info().listOfEntryPressure[j]
				    && !hasInterfaceList[listOfPores[i]->info().poreNeighbors[j]]
				    && !listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().isWResInternal
				    && saturationList[listOfPores[i]->info().poreNeighbors[j]]
				            > listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().minSaturation
				    && saturationList[listOfPores[i]->info().poreNeighbors[j]] <= 1.0) {
					hasInterfaceList[listOfPores[i]->info().poreNeighbors[j]] = true;
					saturationList[listOfPores[i]->info().poreNeighbors[j]]   = 1.0 - truncationPrecision;
					pressuresList[listOfPores[i]->info().poreNeighbors[j]] = porePressureFromPcS(listOfPores[i], 1.0 - truncationPrecision);
					listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().invadedFrom = i;
				}
			}
		}
	}


	//truncate saturation
	for (unsigned int i = 0; i < numberOfPores; i++) {
		if ((saturationList[i] < truncationPrecision || saturationList[i] <= listOfPores[i]->info().minSaturation)) {
			waterVolumeTruncatedLost -= (listOfPores[i]->info().minSaturation - saturationList[i]) * listOfPores[i]->info().mergedVolume;
			saturationList[i] = listOfPores[i]->info().minSaturation;
			// 	    hasInterfaceList[i] = false; // NOTE: in case of deactivation of empty cell, set hasInterfaceList[i] to false
			pressuresList[i]               = porePressureFromPcS(listOfPores[i], listOfPores[i]->info().minSaturation); //waterBoundaryPressure;
			listOfPores[i]->info().isNWRes = true;
			for (unsigned int j = 0; j < listOfPores[i]->info().poreNeighbors.size(); j++) {
				if (!hasInterfaceList[listOfPores[i]->info().poreNeighbors[j]]
				    && !listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().isNWRes
				    && listOfFlux[listOfPores[i]->info().poreNeighbors[j]] >= 0.0
				    && !listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().isWResInternal) {
					hasInterfaceList[listOfPores[i]->info().poreNeighbors[j]] = true;
					saturationList[listOfPores[i]->info().poreNeighbors[j]]   = 1.0 - truncationPrecision;
					pressuresList[listOfPores[i]->info().poreNeighbors[j]]
					        = porePressureFromPcS(listOfPores[listOfPores[i]->info().poreNeighbors[j]], 1.0 - truncationPrecision);
				}
			}
		}

		if (listOfPores[i]->info().isNWResDef && saturationList[i] < listOfPores[i]->info().thresholdSaturation) {
			listOfPores[i]->info().isNWResDef = false;
		}
	}


	if (deformation) {
		for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
			cell->info().poreBodyVolume += cell->info().dv() * oldDT;
		}
		//copyPoreDataToCells();  //NOTE: For two-way coupling this function should be activated, but it is a bit costly for computations.
	}
	for (unsigned int i = 0; i < numberOfPores; i++) {
		if (deformation) {
			listOfPores[i]->info().mergedVolume += listOfPores[i]->info().accumulativeDV * oldDT;
			listOfPores[i]->info().poreBodyRadius
			        = getChi(listOfPores[i]->info().numberFacets) * math::pow(listOfPores[i]->info().mergedVolume, (1. / 3.));
		}
		if (saturationList[i] > 1.0) {
			std::cerr << endl << "Error!, saturation larger than 1? ";
			saturationList[i] = 1.0; //NOTE ADDED AFTER TRUNK UPDATE should be 0.0?
			                         // 		 stopSimulation = true;
		}
		listOfPores[i]->info().saturation     = saturationList[i];
		listOfPores[i]->info().p()            = pressuresList[i];
		listOfPores[i]->info().hasInterface   = bool(hasInterfaceList[i]);
		listOfPores[i]->info().flux           = listOfFlux[i];
		listOfPores[i]->info().dv()           = 0.0; //NOTE ADDED AFTER TRUNK UPDATE
		listOfPores[i]->info().accumulativeDV = 0.0;
	}
}

void TwoPhaseFlowEngine::getQuantities()
{
	Real waterVolume = 0.0, pressureWaterVolume = 0.0, waterVolume_NHJ = 0.0, pressureWaterVolume_NHJ = 0.0, waterVolumeP = 0.0, YDimension = 0.0,
	     simplePressureAverage = 0.0;
	voidVolume                 = 0.0;
	for (unsigned int i = 0; i < numberOfPores; i++) {
		voidVolume += listOfPores[i]->info().mergedVolume;
		waterVolume += listOfPores[i]->info().mergedVolume * listOfPores[i]->info().saturation;
		YDimension += solver->cellBarycenter(listOfPores[i])[1] * listOfPores[i]->info().mergedVolume * listOfPores[i]->info().saturation;
		simplePressureAverage += listOfPores[i]->info().mergedVolume * listOfPores[i]->info().p();

		if (math::abs(listOfPores[i]->info().p()) < 1e10) {
			pressureWaterVolume += listOfPores[i]->info().mergedVolume * listOfPores[i]->info().saturation * listOfPores[i]->info().p();
			waterVolumeP += listOfPores[i]->info().mergedVolume * listOfPores[i]->info().saturation;
		}
		if (listOfPores[i]->info().saturation < 1.0) {
			waterVolume_NHJ += listOfPores[i]->info().mergedVolume * listOfPores[i]->info().saturation;
			pressureWaterVolume_NHJ += listOfPores[i]->info().mergedVolume * listOfPores[i]->info().saturation * listOfPores[i]->info().p();
		}
	}

	Real areaAveragedPressureAcc = 0.0, areaSphere = 0.0;
	airWaterInterfacialArea = 0.0;
	for (unsigned int i = 0; i < numberOfPores; i++) {
		if (listOfPores[i]->info().hasInterface) {
			if (listOfPores[i]->info().saturation < 1.0 && listOfPores[i]->info().saturation >= listOfPores[i]->info().thresholdSaturation) {
				areaSphere = 4.0 * 3.14159265359
				        * math::pow(getChi(listOfPores[i]->info().numberFacets)
				                            * math::pow(
				                                    listOfPores[i]->info().mergedVolume * (1.0 - listOfPores[i]->info().saturation), 0.3333),
				                    2);
			}
			if (listOfPores[i]->info().saturation < listOfPores[i]->info().thresholdSaturation && listOfPores[i]->info().saturation > 0.0
			    && listOfPores[i]->info().saturation > listOfPores[i]->info().minSaturation) { //FIXME FIXME FIXME 1 june 2016
				areaSphere = 4.0 * 3.14159265359 * math::pow((2.0 * surfaceTension / (-1.0 * listOfPores[i]->info().p())), 2.0)
				        + 2.0 * getN(listOfPores[i]->info().numberFacets)
				                * (listOfPores[i]->info().poreBodyRadius - (2.0 * surfaceTension / (-1.0 * listOfPores[i]->info().p())))
				                * (2.0 * surfaceTension / (-1.0 * listOfPores[i]->info().p()))
				                * (2.0 * 3.14159265359 - getDihedralAngle(listOfPores[i]->info().numberFacets));
			}
			areaAveragedPressureAcc += areaSphere * listOfPores[i]->info().p();
			airWaterInterfacialArea += areaSphere;
		}
	}

	areaAveragedPressure           = areaAveragedPressureAcc / airWaterInterfacialArea;
	waterSaturation                = waterVolume / voidVolume;
	waterPressure                  = pressureWaterVolume / waterVolumeP;
	waterPressurePartiallySatPores = pressureWaterVolume_NHJ / waterVolume_NHJ;
	simpleWaterPressure            = simplePressureAverage / voidVolume;
	totalWaterVolume               = waterVolume;


	if (!deformation) {
		Real volumeWaterVBC = 0.0, volumeWaterPressureBC = 0.0, volumeWaterBC = 0.0, volumeWaterAirBC = 0.0, volumeWaterPressureAirBC = 0.0,
		     volumeAirBC = 0.0, Ybottom = 0.0, Ytop = 0.0;
		for (unsigned int i = 0; i < numberOfPores; i++) {
			if (listOfPores[i]->info().waterBC) {
				volumeWaterVBC += listOfPores[i]->info().saturation * listOfPores[i]->info().mergedVolume;
				volumeWaterPressureBC += listOfPores[i]->info().saturation * listOfPores[i]->info().mergedVolume * listOfPores[i]->info().p();
				volumeWaterBC += listOfPores[i]->info().mergedVolume;
				Ybottom += solver->cellBarycenter(listOfPores[i])[1] * listOfPores[i]->info().mergedVolume;
			}
			if (listOfPores[i]->info().airBC) {
				volumeWaterAirBC += listOfPores[i]->info().saturation * listOfPores[i]->info().mergedVolume;
				volumeWaterPressureAirBC
				        += listOfPores[i]->info().saturation * listOfPores[i]->info().mergedVolume * listOfPores[i]->info().p();
				volumeAirBC += listOfPores[i]->info().mergedVolume;
				Ytop += solver->cellBarycenter(listOfPores[i])[1] * listOfPores[i]->info().mergedVolume;
			}
		}
		Real Stop    = volumeWaterAirBC / volumeAirBC; //air BC
		Real Sbottom = volumeWaterVBC / volumeWaterBC; //Water BC.
		Real Ptop    = volumeWaterPressureAirBC / volumeWaterAirBC;
		Real Pbottom = volumeWaterPressureBC / volumeWaterVBC;
		Real z       = (((Ytop / volumeAirBC) - (Ybottom / volumeWaterBC)) / 2.0) + (Ybottom / volumeWaterBC);
		Real gradP   = -1.0 * (Stop - Sbottom) + (Stop * Ptop - Sbottom * Pbottom);
		Real gradZ   = -1.0 * (YDimension / waterVolume) * (Stop - Sbottom) + ((Stop * Ytop / volumeAirBC) - (Sbottom * Ybottom / volumeWaterBC));
		centroidAverageWaterPressure = waterPressure + (1.0 / gradZ) * (z - (YDimension / waterVolume)) * gradP;
	}
}

void TwoPhaseFlowEngine::imposeDeformationFluxTPF()
{
	RTriangulation&     tri     = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		cell->info().dv() = cell->info().dvTPF; //Only relevant for imposed deformation from python-shell
	}
	imposeDeformationFluxTPFSwitch = true;
}

void TwoPhaseFlowEngine::updateDeformationFluxTPF()
{
	RTriangulation&     tri     = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	Real                dv = 0.0, summ = 0.0, summBC = 0.0, summMC = 0.0, SolidVolume = 0.0, dvSwelling = 0.0;

	if (!imposeDeformationFluxTPFSwitch) {
		setPositionsBuffer(true);
		updateVolumes(*solver);


		if (swelling) {
			Real volume = 0.0, invTime = (1.0 / scene->dt);
			if (scene->dt == 0.0) {
				std::cerr << " No dt found!";
			}
			for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
				cell->info().dv() = 0.0;
				if (!cell->info().isFictious) {
					Real solidVol = getSolidVolumeInCell(cell);
					if (solidVol < 0.0) {
						std::cerr << "Error! negative pore body volume! " << solidVol;
						solidVol = 0.0;
					}
					volume = cell->info().volume() * cell->info().volumeSign - solidVol;
					if (volume < 0.0) {
						volume                                              = cell->info().poreBodyVolume;
						listOfPores[cell->info().poreId]->info().isNWRes    = true;
						listOfPores[cell->info().poreId]->info().saturation = truncationPrecision;
					}
					if (cell->info().apparentSolidVolume <= 0.0) {
						cell->info().apparentSolidVolume = solidVol;
					}
					cell->info().dvSwelling = (volume - cell->info().poreBodyVolume + solidVol - cell->info().apparentSolidVolume) * invTime
					        - cell->info().dv();
					if (cell->info().isNWRes || listOfPores[cell->info().poreId]->info().isNWRes) {
						cell->info().dvSwelling = 0.0;
					}
					cell->info().dv() = (volume - cell->info().poreBodyVolume) * invTime;
					SolidVolume += solidVol;
					summ += cell->info().dv();
					summBC += cell->info().dv();
				}
			}
		}
	}


	for (unsigned int i = 0; i < numberOfPores; i++) {
		dv = 0.0, dvSwelling = 0.0;
		for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (cell->info().poreId == (int)i) {
				dv += cell->info().dv();
				dvSwelling += cell->info().dvSwelling;
			}
		}

		listOfPores[i]->info().accumulativeDV         = dv;
		listOfPores[i]->info().accumulativeDVSwelling = dvSwelling;
		summMC += dv;
	}

	if (swelling) {
		//Account for swelling of particles into a non-existing pore (i.e. boundary pores).
		for (unsigned int i = 0; i < numberOfPores; i++) {
			if (listOfPores[i]->info().isNWRes) {
				Real count = 0.0;
				for (unsigned int j = 0; j < listOfPores[i]->info().poreNeighbors.size(); j++) {
					if (!listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().isNWRes) {
						count += 1.0;
					}
				}
				for (unsigned int j = 0; j < listOfPores[i]->info().poreNeighbors.size(); j++) {
					if (!listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().isNWRes) {
						if (count != 0.0) {
							listOfPores[listOfPores[i]->info().poreNeighbors[j]]->info().accumulativeDVSwelling
							        += listOfPores[i]->info().accumulativeDVSwelling / count;
						}
					}
				}
				listOfPores[i]->info().accumulativeDVSwelling = 0.0;
			}
		}
	}
}

Real TwoPhaseFlowEngine::getSolidVolumeInCell(CellHandle cell)
{
	//Dublicate function that depends on position buffer of particles
	//FIXME this function be replaced if function of void volume can be made dependent on updated location of particles
	Real Vsolid                      = 0;
	cell->info().apparentSolidVolume = 0.0;
	for (int i = 0; i < 4; i++) {
		const Vector3r& p0v = positionBufferCurrent[cell->vertex(solver->permut4[i][0])->info().id()].pos;
		const Vector3r& p1v = positionBufferCurrent[cell->vertex(solver->permut4[i][1])->info().id()].pos;
		const Vector3r& p2v = positionBufferCurrent[cell->vertex(solver->permut4[i][2])->info().id()].pos;
		const Vector3r& p3v = positionBufferCurrent[cell->vertex(solver->permut4[i][3])->info().id()].pos;
		Point           p0(p0v[0], p0v[1], p0v[2]);
		Point           p1(p1v[0], p1v[1], p1v[2]);
		Point           p2(p2v[0], p2v[1], p2v[2]);
		Point           p3(p3v[0], p3v[1], p3v[2]);
		Real            rad = positionBufferCurrent[cell->vertex(solver->permut4[i][0])->info().id()].radius;

		Real angle                          = solver->fastSolidAngle(p0, p1, p2, p3);
		cell->info().particleSurfaceArea[i] = rad * rad * angle;


		if (setFractionParticles[cell->vertex(i)->info().id()] > 0) { //should be moved
			cell->info().apparentSolidVolume += rad * rad * angle
			        / (setFractionParticles[cell->vertex(i)->info().id()]
			           * setFractionParticles[cell->vertex(i)->info().id()]); //Coupling to account for swelling in dry pores
		}
		Vsolid += (1. / 3.) * math::pow(rad, 3) * math::abs(angle);
	}
	return Vsolid;
}


void TwoPhaseFlowEngine::updatePoreUnitProperties()
{
	//FIXME clean-up this function (computePoreThroatRadiusMethod2() does not include update of particle location, thus this is a quick fix


	RTriangulation&     tri     = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (!cell->info().isFictious) {
			for (unsigned int j = 0; j < 4; j++) {
				if (cell->info().poreId != cell->neighbor(j)->info().poreId && cell->info().id > cell->neighbor(j)->info().id) {
					Real    rA = positionBufferCurrent[cell->vertex(facetVertices[j][0])->info().id()].radius;
					Real    rB = positionBufferCurrent[cell->vertex(facetVertices[j][1])->info().id()].radius;
					Real    rC = positionBufferCurrent[cell->vertex(facetVertices[j][2])->info().id()].radius;
					CVector posA(
					        positionBufferCurrent[cell->vertex(facetVertices[j][0])->info().id()].pos[0],
					        positionBufferCurrent[cell->vertex(facetVertices[j][0])->info().id()].pos[1],
					        positionBufferCurrent[cell->vertex(facetVertices[j][0])->info().id()].pos[2]);
					CVector posB(
					        positionBufferCurrent[cell->vertex(facetVertices[j][1])->info().id()].pos[0],
					        positionBufferCurrent[cell->vertex(facetVertices[j][1])->info().id()].pos[1],
					        positionBufferCurrent[cell->vertex(facetVertices[j][1])->info().id()].pos[2]);
					CVector posC(
					        positionBufferCurrent[cell->vertex(facetVertices[j][2])->info().id()].pos[0],
					        positionBufferCurrent[cell->vertex(facetVertices[j][2])->info().id()].pos[1],
					        positionBufferCurrent[cell->vertex(facetVertices[j][2])->info().id()].pos[2]);
					CVector B = posB
					        - posA; //positionBufferCurrent[cell->vertex(facetVertices[j][1])->info().id()].pos - positionBufferCurrent[cell->vertex(facetVertices[j][0])->info().id()].pos;
					CVector x = B / sqrt(B.squared_length());
					CVector C = posC
					        - posA; //positionBufferCurrent[cell->vertex(facetVertices[j][2])->info().id()].pos - positionBufferCurrent[cell->vertex(facetVertices[j][0])->info().id()].pos;
					CVector z = CGAL::cross_product(x, C);
					CVector y = CGAL::cross_product(x, z);
					y         = y / math::sqrt(y.squared_length());

					Real b1[2];
					b1[0] = B * x;
					b1[1] = B * y;
					Real c1[2];
					c1[0] = C * x;
					c1[1] = C * y;

					Real A = ((math::pow(rA, 2)) * (1 - c1[0] / b1[0]) + ((math::pow(rB, 2) * c1[0]) / b1[0]) - math::pow(rC, 2)
					          + pow(c1[0], 2) + math::pow(c1[1], 2) - ((math::pow(b1[0], 2) + math::pow(b1[1], 2)) * c1[0] / b1[0]))
					        / (2 * c1[1] - 2 * b1[1] * c1[0] / b1[0]);
					Real BB = (rA - rC - ((rA - rB) * c1[0] / b1[0])) / (c1[1] - b1[1] * c1[0] / b1[0]);
					Real CC = (math::pow(rA, 2) - math::pow(rB, 2) + math::pow(b1[0], 2) + math::pow(b1[1], 2)) / (2 * b1[0]);
					Real D  = (rA - rB) / b1[0];
					Real E  = b1[1] / b1[0];
					Real F  = math::pow(CC, 2) + math::pow(E, 2) * math::pow(A, 2) - 2 * CC * E * A;

					Real c = -F - math::pow(A, 2) + pow(rA, 2);
					Real b = 2 * rA - 2 * (D - BB * E) * (CC - E * A) - 2 * A * BB;
					Real a = 1 - math::pow((D - BB * E), 2) - math::pow(BB, 2);

					if ((math::pow(b, 2) - 4 * a * c) < 0) {
						std::cout << "NEGATIVE DETERMINANT" << endl;
					}
					Real reff = (-b + math::sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);


					if (cell->vertex(facetVertices[j][2])->info().isFictious || cell->vertex(facetVertices[j][1])->info().isFictious
					    || cell->vertex(facetVertices[j][2])->info().isFictious) {
						reff = -1 * reff;
					}

					cell->info().poreThroatRadius[j]                                      = reff;
					cell->neighbor(j)->info().poreThroatRadius[tri.mirror_index(cell, j)] = reff;
				}
			}
		}
	}


	makeListOfPoresInCells(true);
}

void TwoPhaseFlowEngine::makeListOfPoresInCells(bool fast)
{
	RTriangulation&     tri     = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	bool                cancel = false, firstCheck = true;
	;
	for (unsigned int j = 0; j < numberOfPores; j++) {
		firstCheck = true;
		std::vector<int>  poreNeighbors;
		std::vector<Real> listOfkNorm;
		std::vector<Real> listOfEntrySaturation;
		std::vector<Real> listOfEntryPressure;
		std::vector<Real> listOfThroatArea;
		for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (cell->info().poreId == (int)j) {
				for (unsigned int ngb = 0; ngb < 4; ngb++) {
					if (cell->neighbor(ngb)->info().poreId != (int)j && cell->neighbor(ngb)->info().poreId != -1) {
						cancel = false;
						for (unsigned int checkID = 0; checkID < poreNeighbors.size(); checkID++) {
							if (poreNeighbors[checkID] == cell->neighbor(ngb)->info().poreId) {
								cancel = true;
								// 		    std::cerr<<"skipCell";
							}
						}
						if ((firstCheck || !cancel) || poreNeighbors.size() == 0) {
							if (!fast) {
								poreNeighbors.push_back(cell->neighbor(ngb)->info().poreId);
							}
							if (!fast) {
								listOfkNorm.push_back(cell->info().kNorm()[ngb]);
							}
							listOfEntryPressure.push_back(
							        entryMethodCorrection * surfaceTension / cell->info().poreThroatRadius[ngb]);
							Real saturation = poreSaturationFromPcS(cell, -1.0 * cell->info().entryPressure[ngb]);
							listOfEntrySaturation.push_back(saturation);
							if (saturation > 1.0 || saturation < 0.0 || saturation != saturation) {
								std::cerr << endl
								          << "Time to update triangulation, entry saturation not correct: " << saturation;
							}


							if (!fast) {
								const CVector& Surfk = cell->info().facetSurfaces[ngb];
								Real           area  = sqrt(Surfk.squared_length());
								listOfThroatArea.push_back(area * cell->info().facetFluidSurfacesRatio[ngb]);
							}
							if (firstCheck) {
								firstCheck = false;
							}
						}
					}
				}
			}
		}
		if (fast) {
			//                 listOfPores[j]->info().poreNeighbors = poreNeighbors;
			listOfPores[j]->info().listOfEntrySaturation = listOfEntrySaturation;
			listOfPores[j]->info().listOfEntryPressure   = listOfEntryPressure;
			// 		   listOfPores[j]->info().listOfThroatArea = listOfThroatArea;
			// 		   listOfPores[j]->info().listOfkNorm = listOfkNorm;
			// 		   listOfPores[j]->info().listOfkNorm2 = listOfkNorm;
		}

		if (!fast) {
			for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
				if (cell->info().poreId == (int)j) {
					cell->info().poreNeighbors         = poreNeighbors;
					cell->info().listOfEntrySaturation = listOfEntrySaturation;
					cell->info().listOfEntryPressure   = listOfEntryPressure;
					cell->info().listOfThroatArea      = listOfThroatArea;
					cell->info().listOfkNorm           = listOfkNorm;
					cell->info().listOfkNorm2          = listOfkNorm;
				}
			}
		}
	}
}

void TwoPhaseFlowEngine::copyPoreDataToCells()
{
	//NOTE: Don't apply this function via python directly after applying reTriangulation() via Python, this will give segment fault.
	RTriangulation&     tri     = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (!cell->info().isFictious) {
			cell->info().saturation   = listOfPores[cell->info().poreId]->info().saturation;
			cell->info().p()          = listOfPores[cell->info().poreId]->info().p();
			cell->info().hasInterface = bool(hasInterfaceList[cell->info().poreId]);
			cell->info().flux         = listOfFlux[cell->info().poreId];
			cell->info().isNWRes      = listOfPores[cell->info().poreId]->info().isNWRes;
			cell->info().airWaterArea = listOfPores[cell->info().poreId]->info().airWaterArea;
			if (deformation) {
				cell->info().mergedVolume = listOfPores[cell->info().poreId]->info().mergedVolume; //NOTE ADDED AFTER TRUNK UPDATE
				cell->info().poreBodyRadius
				        = getChi(cell->info().numberFacets) * math::pow(listOfPores[cell->info().poreId]->info().mergedVolume, (1. / 3.));
			} //NOTE ADDED AFTER TRUNK UPDATE
			  //NOTE ADDED AFTER TRUNK UPDATE
		}
	}
}

void TwoPhaseFlowEngine::actionTPF()
{
	iterationTPF += 1;
	if (firstDynTPF) {
		std::cout << endl
		          << "Welcome to the two-phase flow Engine" << endl
		          << "by T.Sweijen, B.Chareyre and S.M.Hassanizadeh" << endl
		          << "For contact: T.Sweijen@uu.nl";
		solver->computePermeability();
		scene->time = 0.0;
		initialization();
		actionMergingAlgorithm();
		calculateResidualSaturation();
		setInitialConditions();
		setBoundaryConditions();
		verifyCompatibilityBC();
		setPoreNetwork();
		scene->dt = 1e-20;
		setListOfPores();
		solvePressure();
		getQuantities();
		firstDynTPF = false;
	}
	if (!firstDynTPF && !stopSimulation) {
		//    bool remesh = false;
		//Time steps + deformation, but no remeshing
		scene->time = scene->time + scene->dt;
		if (deformation && !remesh) {
			updateDeformationFluxTPF();
			if (int(float(iterationTPF) / 10.0) == float(iterationTPF) / 10.0) {
				updatePoreUnitProperties();
			}
		}
		//Update pore throat radii etc.

		if (deformation && remesh) {
			reTriangulate(); //retriangulation + merging
			calculateResidualSaturation();
			transferConditions(); //get saturation, hasInterface from previous network
			setBoundaryConditions();
			setPoreNetwork();
		}
		setListOfPores();
		if (solvePressureSwitch) {
			solvePressure();
		}
		if (deformation) {
			if (int(float(iterationTPF) / 50.0) == float(iterationTPF) / 50.0) {
				getQuantities();
			}
		} //FIXME update of quantities has to be made more appropiate


		//    getQuantities();//NOTE FIX

		if (!deformation) {
			if (!getQuantitiesUpdateCont) {
				if (int(float(iterationTPF) / 100.0) == float(iterationTPF) / 100.0) {
					getQuantities();
				}
			}
			if (getQuantitiesUpdateCont) {
				getQuantities();
			}
		}


		if (remesh) {
			remesh = false;
		} //Remesh bool is also used in solvePressure();
	}
}


void TwoPhaseFlowEngine::updateReservoirLabel()
{
	RTriangulation& tri = solver->T[solver->currentTes].Triangulation();
	if (clusters.size() < 2) {
		clusters.resize(2);
		clusters[0] = shared_ptr<PhaseCluster>(new PhaseCluster(solver->tesselation()));
		clusters[1] = shared_ptr<PhaseCluster>(new PhaseCluster(solver->tesselation()));
	}
	clusters[0]->reset();
	clusters[0]->label = 0;
	clusters[1]->reset();
	clusters[1]->label          = 1;
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (cell->info().isNWRes)
			clusterGetPore(clusters[0].get(), cell);
		else if (cell->info().isWRes) {
			clusterGetPore(clusters[1].get(), cell);
			for (int facet = 0; facet < 4; facet++)
				if ((not tri.is_infinite(cell->neighbor(facet))) and !cell->neighbor(facet)->info().isWRes)
					clusterGetFacet(clusters[1].get(), cell, facet);
		} else if (cell->info().label > 1)
			continue;
		else
			cell->info().label = -1;
	}
}

void TwoPhaseFlowEngine::clusterGetFacet(PhaseCluster* cluster, CellHandle cell, int facet)
{
	cell->info().hasInterface = true;
	Real interfArea           = sqrt((cell->info().facetSurfaces[facet] * cell->info().facetFluidSurfacesRatio[facet]).squared_length());
	cluster->interfaces.push_back(PhaseCluster::Interface(std::pair<std::pair<unsigned int, unsigned int>, Real>(
	        std::pair<unsigned int, unsigned int>(cell->info().id, cell->neighbor(facet)->info().id), interfArea)));
	cluster->interfaces.back().outerIndex = facet;
	cluster->interfaces.back().innerCell  = cell;
	cluster->interfacialArea += interfArea;
	if (cluster->entryRadius < cell->info().poreThroatRadius[facet]) {
		cluster->entryRadius = cell->info().poreThroatRadius[facet];
		cluster->entryPore   = cell->info().id;
	}
}

void TwoPhaseFlowEngine::clusterGetPore(PhaseCluster* cluster, CellHandle cell)
{
	cell->info().label      = cluster->label;
	cell->info().saturation = cluster->label == 0 ? 0 : 1;
	cell->info().isNWRes    = cluster->label == 0 ? true : false;
	cell->info().isWRes     = cluster->label == 0 ? false : true;
	cluster->volume += cell->info().poreBodyVolume;
	cluster->pores.push_back(cell);
}

vector<int> TwoPhaseFlowEngine::clusterOutvadePore(unsigned startingId, unsigned imbibedId, int /*index*/)
{
	CellHandle&   origin  = solver->tesselation().cellHandles[startingId];
	CellHandle&   newPore = solver->tesselation().cellHandles[imbibedId];
	PhaseCluster* cluster = clusters[origin->info().label].get();
	cluster->resetSolver(); //reset the linear system
	clusterGetPore(cluster, newPore);
	//NOTE: the code below could be a starting point for more efficient removal, it's currently useless (and parameter index as well)
	// Further, removing from lists should be faster than from vectors, OTOH we probably also need access by index.
	/*unsigned facetIdx;
	if (	index>=0 and unsigned(index)<cluster->interfaces.size() and
		cluster->interfaces[index].first.first == startingId and
		cluster->interfaces[index].first.second == imbibedId)  {
		  facetIdx=index;
	} else {
	  if (index>=0) LOG_WARN("index mismatch wrt. cell ids");
	  for (facetIdx=0; cluster->interfaces[facetIdx].first.first != startingId or cluster->interfaces[facetIdx].first.second!=imbibedId; facetIdx++)
		{if ((facetIdx+1)>=cluster->interfaces.size()) LOG_WARN("interface not found");}
	}*/
	bool        updateIntfs = false; //if turned true later we will have to clean interfaces
	vector<int> merged      = { cluster->label };

	for (int k = 0; k < 4; k++) {
		if (INFT(newPore->neighbor(k)) or (newPore->neighbor(k) == origin))
			continue;
		if (newPore->neighbor(k)->info().label == 0)
			clusterGetFacet(cluster, newPore, k);
		else {
			// 			updateIntfs=true;
			if (newPore->neighbor(k)->info().label != cluster->label) {
				merged.push_back(newPore->neighbor(k)->info().label);
				cluster->mergeCluster(*clusters[newPore->neighbor(k)->info().label], newPore);
			} else
				updateIntfs = true; //one more interface needs to be removed, FIXME: this may lead to long copy operations
		}
	}
	if (updateIntfs) {
		for (int k = cluster->interfaces.size() - 1; k >= 0; k--)
			if (solver->tesselation().cellHandles[cluster->interfaces[k].first.second]->info().label == cluster->label)
				cluster->interfaces.erase(cluster->interfaces.begin() + k);
	}
	PhaseCluster* cluster0 = clusters[0].get();
	for (auto p = cluster0->pores.begin();;) { //slow search...
		if (p == cluster0->pores.end()) {
			LOG_WARN("pore " << newPore->info().id << "not found in cluster" << cluster0->label << " of size " << cluster0->pores.size());
			break;
		} else {
			if ((*p) != newPore)
				p++;    // warning: this is not equivalent to p.id==cell.id for some reason, some wrong positive it seems
				        // 			if ((*p)->info().id!=cell->info().id) p++;
			else {
				cluster0->pores.erase(p);
				break;
			}
		}
	}

	for (int k = cluster->interfaces.size() - 1; k >= 0; k--)
		if (solver->tesselation().cellHandles[cluster->interfaces[k].first.second]->info().label == cluster->label) {
			//TODO: what happens to the other interfaces on the same pore? they will be removed but should they give bridges or something?
			// 			cluster->interfacialArea-=cluster->interfaces[k].second;
			cluster->interfaces.erase(cluster->interfaces.begin() + k);
		}
	return merged;
}

vector<int> TwoPhaseFlowEngine::clusterInvadePore(PhaseCluster* cluster, CellHandle cell)
{
	//invade the pore and attach to NW reservoir, label is assigned after reset
	int label               = cell->info().label;
	cell->info().saturation = 0;
	cell->info().isNWRes    = true;
	cell->info().isWRes     = false;
	clusterGetPore(clusters[0].get(), cell);

	//update the cluster(s)
	unsigned    nPores = cluster->pores.size();
	vector<int> newClusters; //for returning the list of possible sub-clusters, empty if we are removing the last pore of the base cluster
	if (nPores == 0) {
		LOG_WARN("Invading the empty cluster id=" << label);
	}
	if (nPores == 1) {
		cluster->reset();
		cluster->label = label;
		return newClusters;
	}
	FOREACH(CellHandle & cell2, cluster->pores) { cell2->info().label = -1; } //mark all pores, and get them back in again below
	cell->info().label = 0;                                                   //mark the invaded one

	//find a remaining pore
	unsigned neighborStart = 0;
	while ((cell->neighbor(neighborStart)->info().label != -1 or solver->T[solver->currentTes].Triangulation().is_infinite(cell->neighbor(neighborStart)))
	       and neighborStart < 3)
		++neighborStart;
	if (neighborStart == 3 and cell->neighbor(neighborStart)->info().label != -1)
		cerr << "This is not supposed to happen (line " << __LINE__ << ")" << endl;

	auto nCell          = cell->neighbor(neighborStart); //use the remaining pore to start reconstruction of the cluster
	nCell->info().label = label;                         //assign the label of the original cluster
	cluster->reset();                                    //reset pores, volume, entryRadius, area... but restore label again after that
	cluster->label = label;
	updateSingleCellLabelRecursion(nCell, cluster); //rebuild
	newClusters.push_back(cluster->label);          //we will return the original cluster itself if not empty

	// gen new clusters on the fly from the other neighbors of the invaded pore (for disconnected subclusters)
	for (int neighborId = neighborStart + 1; neighborId <= 3; neighborId++) { //should be =1 if the cluster remain the same -1 removed pore
		const CellHandle& nCell2 = cell->neighbor(neighborId);
		if (nCell2->info().label != -1 or solver->T[solver->currentTes].Triangulation().is_infinite(nCell2))
			continue; //already reached from another neighbour (connected domain): skip, else this is a new cluster
		shared_ptr<PhaseCluster> clst(new PhaseCluster(solver->tesselation()));
		clst->label = clusters.size();
		newClusters.push_back(clst->label);
		clusters.push_back(clst);
		updateSingleCellLabelRecursion(nCell2, clusters.back().get());
	}
	return newClusters; // return list of created clusters
}

bool TwoPhaseFlowEngine::connectedAroundEdge(const RTriangulation& Tri, CellHandle& cell, unsigned facet1, unsigned facet2)
{
	revertEdge(facet1, facet2);
	RTriangulation::Cell_circulator cell1 = Tri.incident_cells(cell, facet1, facet2, cell);
	RTriangulation::Cell_circulator cell0 = cell1++;
	const int&                      label = cell1->info().label;
	while (cell1 != cell0 and !Tri.is_infinite(cell1) and (cell1->info().label == label))
		cell1++; //around edge looking for contiguous labels
	return (cell1 == cell0);
}

vector<int> TwoPhaseFlowEngine::clusterInvadePoreFast(PhaseCluster* cluster, CellHandle cell)
{
	//invade the pore and attach to NW reservoir, label is assigned after reset
	int label = cell->info().label;
	if (label != cluster->label)
		LOG_WARN("wrong label");
	if (cell->info().Pcondition) {
		if (solver->debugOut)
			LOG_WARN("invading a Pcondition pore (ignored)");
		return vector<int>(1, label);
	}
	const RTriangulation& Tri = solver->T[solver->currentTes].Triangulation();
#ifdef LINSOLV
	cluster->resetSolver();
#endif
	unsigned id             = cell->info().id;
	cell->info().saturation = 0;
	cell->info().isNWRes    = true;
	cell->info().isWRes     = false;
	// 	cell->info().Pcondition=true;
	clusterGetPore(clusters[0].get(), cell); //this will update cell label as well
	//update the cluster(s)
	unsigned    nPores = cluster->pores.size();
	vector<int> newClusters; //for returning the list of possible sub-clusters, empty if we are removing the last pore of the base cluster
	if (nPores == 0) {
		LOG_WARN("Invading an empty cluster id=" << label);
		return newClusters;
	}
	if (nPores == 1) {
		LOG_WARN("Invading last pore of cluster id=" << label);
		cluster->reset();
		cluster->label = label;
		return newClusters;
	}
	//count neighbors from the same cluster
	vector<CellHandle> clustNeighbors;
	vector<unsigned>   clustNIdx;
	for (int k = 0; k < 4; k++)
		if (not INFT(cell->neighbor(k)) and cell->neighbor(k)->info().label == label) {
			clustNeighbors.push_back(cell->neighbor(k));
			clustNIdx.push_back(k);
		}
	unsigned nN     = clustNeighbors.size();
	unsigned nFaces = 4 - nN;
	//update interfaces, first remove old ones
	for (int k = cluster->interfaces.size() - 1; (k >= 0 and nFaces > 0); k--)
		if (cluster->interfaces[k].first.first == id) {
			//TODO: what happens to the other interfaces on the same pore? they will be removed but should they give bridges or something?
			cluster->interfacialArea -= cluster->interfaces[k].second;
			cluster->interfaces.erase(cluster->interfaces.begin() + k);
			nFaces--;
		}
	//then add new ones (TODO: set capillary pressure and volume?)
	for (auto cn = clustNeighbors.begin(); cn != clustNeighbors.end(); cn++)
		clusterGetFacet(cluster, *cn, (*cn)->index(cell));
	//now remove the invaded pore
	for (auto p = cluster->pores.begin();;) { //slow search...
		if (p == cluster->pores.end()) {
			LOG_WARN("pore " << cell->info().id << "not found in cluster" << cluster->label << " of size " << cluster->pores.size());
			break;
		} else {
			if ((*p) != cell)
				p++;    // warning: this is not equivalent to p.id==cell.id for some reason, some wrong positive it seems
				        // 			if ((*p)->info().id!=cell->info().id) p++;
			else {
				cluster->pores.erase(p);
				break;
			}
		}
	}
	//it could be that the cluster has been splitted in smaller clusters, before going to complex rebuilding method we try to exit the trivial cases below
	newClusters.push_back(cluster->label); //we will return at least the original cluster itself

	//1. case of only one neighbor from cluster connected to the one being erased
	if (nN == 1) { /*LOG_WARN("nN==1 ?!");*/
		return newClusters;
	}
	if (nN == 2) { /*LOG_WARN("nN==2 ?!");*/
		// check if pores connected by the removed one are still directly connected locally around one edge
		unsigned i = clustNIdx[0];
		unsigned j = clustNIdx[1];
		if (connectedAroundEdge(Tri, cell, i, j))
			return newClusters;
	}
	if (nN == 3) { /*LOG_WARN("nN==3 ?!");*/
		// check if pores connected by the removed one are still connected locally around at least two edges
		unsigned i  = clustNIdx[0];
		unsigned j  = clustNIdx[1];
		unsigned k  = clustNIdx[2];
		unsigned ij = unsigned(connectedAroundEdge(Tri, cell, i, j));
		unsigned jk = unsigned(connectedAroundEdge(Tri, cell, j, k));
		unsigned ki = unsigned(connectedAroundEdge(Tri, cell, k, i));
		if ((ij + jk + ki) >= 2) //the cluster is not splitted by the invasion (2 connexions between three cells)
			return newClusters;
	}
	if (nN == 4)
		LOG_WARN(
		        "nN==4 ?! for cell" << cell->info().id << " " << cell->neighbor(0)->info().id << " " << cell->neighbor(1)->info().id << " "
		                            << cell->neighbor(2)->info().id << " " << cell->neighbor(3)->info().id);
	//not a trivial case, go for a split (possibly still resuting in one single cluster)
	CellHandle startingCell = cluster->pores[0]; //will be changed below in the case label=1, to keep cluster 1 as the W-reservoir cluster
	if (cluster->label == 1) {
		bool foundImpP1 = false;
		int  k          = cluster->pores.size() - 1;
		while (not foundImpP1) {
			if (cluster->pores[k]->info().Pcondition) {
				startingCell = cluster->pores[k];
				foundImpP1   = true;
			} else if (k > 0)
				k--;
			else {
				LOG_WARN("no Pcondition pore in cluster 1");
				break;
			}
		}
	}
	return splitCluster(cluster, startingCell); // return list of created clusters
}

vector<int> TwoPhaseFlowEngine::splitCluster(PhaseCluster* cluster, CellHandle cellInit)
{
	unsigned oldSize = cluster->pores.size();
	if (oldSize == 0) {
		LOG_WARN("empty call ");
		return vector<int>();
	}
	unsigned nextLabel = clusters.size();
	FOREACH(CellHandle & cell, cluster->pores) { cell->info().label = nextLabel; } //mark all pores, and get them back in again below
	unsigned nPoresOld = markRecursively(cellInit, cluster->label);
	if (nPoresOld == oldSize) { /*LOG_WARN("no split, return (label="<<cluster->label <<")");*/
		return vector<int>(1, cluster->label);
	}
	clusters.push_back(shared_ptr<PhaseCluster>(new PhaseCluster(*cluster->tes)));
	auto clst   = clusters.back();
	clst->label = nextLabel;

	unsigned countNew = 0;
	for (int k = cluster->pores.size() - 1; k >= 0; k--) {
		const CellHandle& c = cluster->pores[k];
		if (c->info().label == (int)nextLabel) {
			cluster->volume -= c->info().poreBodyVolume;
			clusterGetPore(clst.get(), c);
			cluster->pores.erase(
			        cluster->pores.begin() + k); //FIXME: definitely needs a rebuild of two lists instead of such 'erase', which is very slow
			countNew++;
		}
	}
	for (int k = cluster->interfaces.size() - 1; k >= 0; k--) {
		const CellHandle& c = solver->tesselation().cellHandles[cluster->interfaces[k].first.first];
		if (c->info().label == (int)nextLabel) {
			clst->interfaces.push_back(PhaseCluster::Interface(cluster->interfaces[k]));
			cluster->interfaces.erase(cluster->interfaces.begin() + k);
		}
	}
	if (countNew > 1) {
		vector<int> clusterList = splitCluster(clst.get(), clst->pores[0]);
		clusterList.push_back(cluster->label);
		return clusterList;
	} else {
		vector<int> clusterList = { cluster->label, (int)nextLabel };
		return clusterList;
	}
}

unsigned TwoPhaseFlowEngine::markRecursively(const CellHandle& cell, int newLabel)
{
	// 	LOG_WARN("markRecursively "<<cell->info().id<<" "<<cell->info().label<<" "<<newLabel <<" "<<solver->tesselation().Triangulation().is_infinite(cell)<<" "<<clusters[1]->pores[0]->info().id);
	if (solver->tesselation().Triangulation().is_infinite(cell) or cell->info().label == newLabel)
		return 0;
	int originalLabel  = cell->info().label;
	cell->info().label = newLabel;
	unsigned count     = 1;
	for (int facet = 0; facet < 4; facet++)
		if (cell->neighbor(facet)->info().label == originalLabel)
			count += markRecursively(cell->neighbor(facet), newLabel);
	return count;
}

// int TwoPhaseFlowEngine:: getMaxCellLabel()
// {
//     int maxLabel=-1;
//     RTriangulation& tri = solver->T[solver->currentTes].Triangulation();
//     FiniteCellsIterator cellEnd = tri.finite_cells_end();
//     for ( FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++ ) {
//         if (cell->info().label>maxLabel) maxLabel=cell->info().label;
//     }
//     return maxLabel;
// }

void TwoPhaseFlowEngine::updateCellLabel()
{
	//     int currentLabel = getMaxCellLabel();//FIXME: A loop on cells for each new label?? is it serious??
	updateReservoirLabel();
	int                 currentLabel = clusters.size();
	RTriangulation&     tri          = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd      = tri.finite_cells_end();
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (cell->info().label == -1) {
			shared_ptr<PhaseCluster> clst(new PhaseCluster(solver->tesselation()));
			clst->label = currentLabel;
			clusters.push_back(clst);
			updateSingleCellLabelRecursion(cell, clusters.back().get());
			currentLabel++;
		}
	}
}

void TwoPhaseFlowEngine::updateSingleCellLabelRecursion(CellHandle cell, PhaseCluster* cluster)
{
	clusterGetPore(cluster, cell);
	//     cell->info().label=label;
	//     cluster->volume+=cell->info().
	//     cluster->pores.push_back(cell);
	for (int facet = 0; facet < 4; facet++) {
		CellHandle nCell = cell->neighbor(facet);
		if (solver->T[solver->currentTes].Triangulation().is_infinite(nCell))
			continue;
		//         if (nCell->info().Pcondition) continue;
		//         if ( (nCell->info().isFictious) && (!isInvadeBoundary) ) continue;
		//TODO:the following condition may relax to relate to nCell->info().hasInterface
		if ((nCell->info().saturation == cell->info().saturation) && (nCell->info().label != cell->info().label))
			updateSingleCellLabelRecursion(nCell, cluster);
		else if (nCell->info().isNWRes)
			clusterGetFacet(cluster, cell, facet);
	}
}


boost::python::list TwoPhaseFlowEngine::pyClusters()
{
	boost::python::list ret;
	for (vector<shared_ptr<PhaseCluster>>::iterator it = clusters.begin(); it != clusters.end(); ++it)
		ret.append(*it);
	return ret;
}

void TwoPhaseFlowEngine::updatePressure()
{
	boundaryConditions(*solver);
	solver->pressureChanged = true;
	solver->reApplyBoundaryConditions();
	RTriangulation&     tri     = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (cell->info().isWRes == true) {
			cell->info().p() = bndCondValue[2];
		}
		if (cell->info().isNWRes == true) {
			cell->info().p() = bndCondValue[3];
		}
		if (isPhaseTrapped) {
			if (cell->info().isTrapW) {
				cell->info().p() = bndCondValue[3] - cell->info().trapCapP;
			}
			if (cell->info().isTrapNW) {
				cell->info().p() = bndCondValue[2] + cell->info().trapCapP;
			}
			//check cell reservoir info.
			if (!cell->info().isWRes && !cell->info().isNWRes && !cell->info().isTrapW && !cell->info().isTrapNW) {
				cerr << "ERROR! NOT FIND Cell Info!";
			}
			// 	{cell->info().p()=bndCondValue[2]; if (isInvadeBoundary) cerr<<"Something wrong in updatePressure.(isInvadeBoundary)";}
		}
	}
}

void TwoPhaseFlowEngine::invasion()
{
	if (isPhaseTrapped)
		invasion1();
	else
		invasion2();
}

///mode1 and mode2 can share the same invasionSingleCell(), invasionSingleCell() ONLY change neighbor pressure and neighbor saturation, independent of reservoirInfo.
void TwoPhaseFlowEngine::invasionSingleCell(CellHandle cell)
{
	Real localPressure   = cell->info().p();
	Real localSaturation = cell->info().saturation;
	if (useFastInvasion and cell->info().label > 0)
		clusterInvadePoreFast(clusters[cell->info().label].get(), cell);
	for (int facet = 0; facet < 4; facet++) {
		CellHandle nCell = cell->neighbor(facet);
		if (solver->T[solver->currentTes].Triangulation().is_infinite(nCell))
			continue;
		if (nCell->info().Pcondition)
			continue; //FIXME:defensive
			          //         if ( (nCell->info().isFictious) && (!isInvadeBoundary) ) continue;
		if (cell->info().poreThroatRadius[facet] < 0)
			continue;

		if ((nCell->info().saturation == localSaturation) && (nCell->info().p() != localPressure)
		    && ((nCell->info().isTrapNW) || (nCell->info().isTrapW))) {
			nCell->info().p() = localPressure;
			if (solver->debugOut) {
				cerr << "merge trapped phase" << endl;
			}
			invasionSingleCell(nCell);
		} ///here we merge trapped phase back to reservoir
		else if ((nCell->info().saturation > localSaturation)) {
			Real nPcThroat = surfaceTension / cell->info().poreThroatRadius[facet];
			Real nPcBody   = surfaceTension / nCell->info().poreBodyRadius;
			if ((localPressure - nCell->info().p() > nPcThroat) && (localPressure - nCell->info().p() > nPcBody)) {
				nCell->info().p()          = localPressure;
				nCell->info().saturation   = localSaturation;
				nCell->info().hasInterface = false;
				if (solver->debugOut) {
					cerr << "drainage" << endl;
				}
				if (recursiveInvasion)
					invasionSingleCell(nCell);
			}
			////FIXME:Introduce cell.hasInterface
			// 	  else if( (localPressure-nCell->info().p()>nPcThroat) && (localPressure-nCell->info().p()<nPcBody) && (cell->info().hasInterface==false) && (nCell->info().hasInterface==false) ) {
			// 	    if(solver->debugOut) {cerr<<"invasion paused into pore interface "<<endl;}
			// 	    nCell->info().hasInterface=true;
			// 	  }
			// 	  else continue;
		} else if ((nCell->info().saturation < localSaturation)) {
			Real nPcThroat = surfaceTension / cell->info().poreThroatRadius[facet];
			Real nPcBody   = surfaceTension / nCell->info().poreBodyRadius;
			if ((nCell->info().p() - localPressure < nPcBody) && (nCell->info().p() - localPressure < nPcThroat)) {
				nCell->info().p()        = localPressure;
				nCell->info().saturation = localSaturation;
				if (solver->debugOut) {
					cerr << "imbibition" << endl;
				}
				if (recursiveInvasion)
					invasionSingleCell(nCell);
			}
			//// FIXME:Introduce cell.hasInterface
			// 	  else if ( (nCell->info().p()-localPressure<nPcBody) && (nCell->info().p()-localPressure>nPcThroat) /*&& (cell->info().hasInterface==false) && (nCell->info().hasInterface==false)*/ ) {
			// 	    nCell->info().p() = localPressure;
			// 	    nCell->info().saturation=localSaturation;
			// 	    if(solver->debugOut) {cerr<<"imbibition paused pore interface"<<endl;}
			// 	    nCell->info().hasInterface=true;
			// 	  }
			// 	  else continue;
		} else
			continue;
	}
}
///invasion mode 1: withTrap
void TwoPhaseFlowEngine::invasion1()
{
	if (solver->debugOut) {
		cout << "----start invasion1----" << endl;
	}

	///update Pw, Pn according to reservoirInfo.
	updatePressure();
	if (solver->debugOut) {
		cout << "----invasion1.updatePressure----" << endl;
	}

	///invasionSingleCell by Pressure difference, change Pressure and Saturation.
	RTriangulation&     tri     = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	if (isDrainageActivated) {
		for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (cell->info().isNWRes)
				invasionSingleCell(cell);
		}
	}
	if (isImbibitionActivated) {
		for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (cell->info().isWRes)
				invasionSingleCell(cell);
		}
	}
	if (solver->debugOut) {
		cout << "----invasion1.invasionSingleCell----" << endl;
	}

	///update W, NW reservoirInfo according to cell->info().saturation
	updateReservoirs1();
	if (solver->debugOut) {
		cout << "----invasion1.update W, NW reservoirInfo----" << endl;
	}

	///search new trapped W-phase/NW-phase, assign trapCapP, isTrapW/isTrapNW flag for new trapped phases. But at this moment, the new trapped W/NW cells.P= W/NW-Res.P. They will be updated in next updatePressure() func.
	checkTrap(bndCondValue[3] - bndCondValue[2]);
	if (solver->debugOut) {
		cout << "----invasion1.checkWTrap----" << endl;
	}

	///update trapped W-phase/NW-phase Pressure //FIXME: is this necessary?
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (cell->info().isTrapW) {
			cell->info().p() = bndCondValue[3] - cell->info().trapCapP;
		}
		if (cell->info().isTrapNW) {
			cell->info().p() = bndCondValue[2] + cell->info().trapCapP;
		}
	}
	if (solver->debugOut) {
		cout << "----invasion1.update trapped W-phase/NW-phase Pressure----" << endl;
	}

	if (isCellLabelActivated and !useFastInvasion)
		updateCellLabel();
	if (solver->debugOut) {
		cout << "----update cell labels----" << endl;
	}
}

///search trapped W-phase or NW-phase, define trapCapP=Pn-Pw. assign isTrapW/isTrapNW info.
void TwoPhaseFlowEngine::checkTrap(Real pressure)
{
	RTriangulation&     tri     = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		//       if( (cell->info().isFictious) && (!cell->info().Pcondition) && (!isInvadeBoundary) ) continue;
		if ((cell->info().isWRes) || (cell->info().isNWRes) || (cell->info().isTrapW) || (cell->info().isTrapNW))
			continue;
		cell->info().trapCapP = pressure;
		if (cell->info().saturation == 1.0)
			cell->info().isTrapW = true;
		if (cell->info().saturation == 0.0)
			cell->info().isTrapNW = true;
	}
}

void TwoPhaseFlowEngine::updateReservoirs1()
{
	RTriangulation&     tri     = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (cell->info().Pcondition)
			continue;
		cell->info().isWRes  = false;
		cell->info().isNWRes = false;
	}

	for (FlowSolver::VCellIterator it = solver->boundingCells[2].begin(); it != solver->boundingCells[2].end(); it++) {
		if ((*it) == NULL)
			continue;
		WResRecursion(*it);
	}

	for (FlowSolver::VCellIterator it = solver->boundingCells[3].begin(); it != solver->boundingCells[3].end(); it++) {
		if ((*it) == NULL)
			continue;
		NWResRecursion(*it);
	}
}

void TwoPhaseFlowEngine::WResRecursion(CellHandle cell)
{
	for (int facet = 0; facet < 4; facet++) {
		CellHandle nCell = cell->neighbor(facet);
		if (solver->T[solver->currentTes].Triangulation().is_infinite(nCell))
			continue;
		if (nCell->info().Pcondition)
			continue;
		//         if ( (nCell->info().isFictious) && (!isInvadeBoundary) ) continue;
		if (nCell->info().saturation != 1.0)
			continue;
		if (nCell->info().isWRes == true)
			continue;
		nCell->info().isWRes   = true;
		nCell->info().isNWRes  = false;
		nCell->info().isTrapW  = false;
		nCell->info().trapCapP = 0.0;
		WResRecursion(nCell);
	}
}

void TwoPhaseFlowEngine::NWResRecursion(CellHandle cell)
{
	for (int facet = 0; facet < 4; facet++) {
		CellHandle nCell = cell->neighbor(facet);
		if (solver->T[solver->currentTes].Triangulation().is_infinite(nCell))
			continue;
		if (nCell->info().Pcondition)
			continue;
		//         if ( (nCell->info().isFictious) && (!isInvadeBoundary) ) continue;
		if (nCell->info().saturation != 0.0)
			continue;
		if (nCell->info().isNWRes == true)
			continue;
		nCell->info().isNWRes  = true;
		nCell->info().isWRes   = false;
		nCell->info().isTrapNW = false;
		nCell->info().trapCapP = 0.0;
		NWResRecursion(nCell);
	}
}

///invasion mode 2: withoutTrap
void TwoPhaseFlowEngine::invasion2()
{
	if (solver->debugOut) {
		cout << "----start invasion2----" << endl;
	}

	///update Pw, Pn according to reservoirInfo.
	updatePressure();
	if (solver->debugOut) {
		cout << "----invasion2.updatePressure----" << endl;
	}

	///drainageSingleCell by Pressure difference, change Pressure and Saturation.
	RTriangulation&     tri     = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	if (isDrainageActivated) {
		for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (cell->info().isNWRes)
				invasionSingleCell(cell);
		}
	}
	///drainageSingleCell by Pressure difference, change Pressure and Saturation.
	if (isImbibitionActivated) {
		for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (cell->info().isWRes)
				invasionSingleCell(cell);
		}
	}

	if (solver->debugOut) {
		cout << "----invasion2.invasionSingleCell----" << endl;
	}

	///update W, NW reservoirInfo according to Pressure
	updateReservoirs2();
	if (solver->debugOut) {
		cout << "----drainage2.update W, NW reservoirInfo----" << endl;
	}
}

void TwoPhaseFlowEngine::updateReservoirs2()
{
	RTriangulation&     tri     = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (cell->info().p() == bndCondValue[2]) {
			cell->info().isWRes  = true;
			cell->info().isNWRes = false;
		} else if (cell->info().p() == bndCondValue[3]) {
			cell->info().isNWRes = true;
			cell->info().isWRes  = false;
		} else {
			cerr << "drainage mode2: updateReservoir Error!" << endl;
		}
	}
}

Real TwoPhaseFlowEngine::getMinDrainagePc()
{
	Real                nextEntry = 1e50;
	RTriangulation&     tri       = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd   = tri.finite_cells_end();
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (cell->info().isNWRes == true) {
			for (int facet = 0; facet < 4; facet++) {
				CellHandle nCell = cell->neighbor(facet);
				if (tri.is_infinite(nCell))
					continue;
				if (nCell->info().Pcondition)
					continue;
				//                 if ( (nCell->info().isFictious) && (!isInvadeBoundary) ) continue;
				if (nCell->info().isWRes == true && cell->info().poreThroatRadius[facet] > 0) {
					Real nCellP = math::max(
					        (surfaceTension / cell->info().poreThroatRadius[facet]), (surfaceTension / nCell->info().poreBodyRadius));
					//                     Real nCellP = surfaceTension/cell->info().poreThroatRadius[facet];
					nextEntry = math::min(nextEntry, nCellP);
				}
			}
		}
	}

	if (nextEntry == 1e50) {
		cout << "End drainage !" << endl;
		return nextEntry = 0;
	} else
		return nextEntry;
}

Real TwoPhaseFlowEngine::getMaxImbibitionPc()
{
	Real                nextEntry = -1e50;
	RTriangulation&     tri       = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd   = tri.finite_cells_end();
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (cell->info().isWRes == true) {
			for (int facet = 0; facet < 4; facet++) {
				CellHandle nCell = cell->neighbor(facet);
				if (tri.is_infinite(nCell))
					continue;
				if (nCell->info().Pcondition)
					continue;
				//                 if ( (nCell->info().isFictious) && (!isInvadeBoundary) ) continue;
				if (nCell->info().isNWRes == true && cell->info().poreThroatRadius[facet] > 0) {
					Real nCellP = math::min(
					        (surfaceTension / nCell->info().poreBodyRadius), (surfaceTension / cell->info().poreThroatRadius[facet]));
					nextEntry = math::max(nextEntry, nCellP);
				}
			}
		}
	}

	if (nextEntry == -1e50) {
		cout << "End imbibition !" << endl;
		return nextEntry = 0;
	} else
		return nextEntry;
}

Real TwoPhaseFlowEngine::getSaturation(bool isSideBoundaryIncluded)
{
	if ((!isInvadeBoundary) && (isSideBoundaryIncluded))
		cerr << "In isInvadeBoundary=false drainage, isSideBoundaryIncluded can't set true." << endl;
	RTriangulation&     tri         = solver->T[solver->currentTes].Triangulation();
	Real                poresVolume = 0.0; //total pores volume
	Real                wVolume     = 0.0; //NW-phase volume
	FiniteCellsIterator cellEnd     = tri.finite_cells_end();

	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (cell->info().Pcondition)
			continue;
		if ((cell->info().isFictious) && (!isSideBoundaryIncluded))
			continue;
		poresVolume = poresVolume + cell->info().poreBodyVolume;
		if (cell->info().saturation > 0.0) {
			wVolume = wVolume + cell->info().poreBodyVolume * cell->info().saturation;
		}
	}
	return wVolume / poresVolume;
}

///compute forces
void TwoPhaseFlowEngine::computeFacetPoreForcesWithCache(bool onlyCache)
{
	RTriangulation& Tri = solver->T[solver->currentTes].Triangulation();
	CVector         nullVect(0, 0, 0);
	//reset forces
	if (!onlyCache)
		for (FiniteVerticesIterator v = Tri.finite_vertices_begin(); v != Tri.finite_vertices_end(); ++v)
			v->info().forces = nullVect;
	// 	#ifdef parallel_forces
	// 	if (solver->noCache) {
	// 		solver->perVertexUnitForce.clear(); solver->perVertexPressure.clear();
	// 		solver->perVertexUnitForce.resize(solver->T[solver->currentTes].maxId+1);
	// 		solver->perVertexPressure.resize(solver->T[solver->currentTes].maxId+1);}
	// 	#endif
	// 	CellHandle neighbourCell;
	// 	VertexHandle mirrorVertex;
	CVector tempVect;
	//FIXME : Ema, be carefull with this (noCache), it needs to be turned true after retriangulation
	if (solver->noCache) { //WARNING:all currentTes must be solver->T[solver->currentTes], should NOT be solver->T[currentTes]
		for (FlowSolver::VCellIterator cellIt = solver->T[solver->currentTes].cellHandles.begin();
		     cellIt != solver->T[solver->currentTes].cellHandles.end();
		     cellIt++) {
			CellHandle& cell = *cellIt;
			//reset cache
			for (int k = 0; k < 4; k++)
				cell->info().unitForceVectors[k] = nullVect;

			for (int j = 0; j < 4; j++)
				if (!Tri.is_infinite(cell->neighbor(j))) {
					const CVector& Surfk = cell->info().facetSurfaces[j];
					//FIXME : later compute that fluidSurf only once in hydraulicRadius, for now keep full surface not modified in cell->info for comparison with other forces schemes
					//The ratio void surface / facet surface
					Real area = sqrt(Surfk.squared_length());
					if (area <= 0)
						cerr << "AREA <= 0!! AREA=" << area << endl;
					CVector                     facetNormal   = Surfk / area;
					const std::vector<CVector>& crossSections = cell->info().facetSphereCrossSections;
					CVector                     fluidSurfk    = cell->info().facetSurfaces[j] * cell->info().facetFluidSurfacesRatio[j];
					/// handle fictious vertex since we can get the projected surface easily here
					if (cell->vertex(j)->info().isFictious) {
						Real projSurf                  = math::abs(Surfk[solver->boundary(cell->vertex(j)->info().id()).coordinate]);
						tempVect                       = -projSurf * solver->boundary(cell->vertex(j)->info().id()).normal;
						cell->vertex(j)->info().forces = cell->vertex(j)->info().forces + tempVect * cell->info().p();
						//define the cached value for later use with cache*p
						cell->info().unitForceVectors[j] = cell->info().unitForceVectors[j] + tempVect;
					}
					/// Apply weighted forces f_k=sqRad_k/sumSqRad*f
					CVector facetUnitForce = -fluidSurfk * cell->info().solidLine[j][3];
					CVector facetForce     = cell->info().p() * facetUnitForce;

					for (int y = 0; y < 3; y++) {
						cell->vertex(facetVertices[j][y])->info().forces
						        = cell->vertex(facetVertices[j][y])->info().forces + facetForce * cell->info().solidLine[j][y];
						//add to cached value
						cell->info().unitForceVectors[facetVertices[j][y]]
						        = cell->info().unitForceVectors[facetVertices[j][y]] + facetUnitForce * cell->info().solidLine[j][y];
						//uncomment to get total force / comment to get only pore tension forces
						if (!cell->vertex(facetVertices[j][y])->info().isFictious) {
							cell->vertex(facetVertices[j][y])->info().forces = cell->vertex(facetVertices[j][y])->info().forces
							        - facetNormal * cell->info().p() * crossSections[j][y];
							//add to cached value
							cell->info().unitForceVectors[facetVertices[j][y]]
							        = cell->info().unitForceVectors[facetVertices[j][y]] - facetNormal * crossSections[j][y];
						}
					}
					// 	#ifdef parallel_forces
					// 	solver->perVertexUnitForce[cell->vertex(j)->info().id()].push_back(&(cell->info().unitForceVectors[j]));
					// 	solver->perVertexPressure[cell->vertex(j)->info().id()].push_back(&(cell->info().p()));
					// 	#endif
				}
		}
		solver->noCache = false; //cache should always be defined after execution of this function
		if (onlyCache)
			return;
	} else { //use cached values when triangulation doesn't change
		 // 		#ifndef parallel_forces
		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); cell++) {
			for (int yy = 0; yy < 4; yy++)
				cell->vertex(yy)->info().forces = cell->vertex(yy)->info().forces + cell->info().unitForceVectors[yy] * cell->info().p();
		}
		/*		#else
		#pragma omp parallel for num_threads(ompThreads)
		for (int vn=0; vn<= solver->T[solver->currentTes].maxId; vn++) {
			VertexHandle& v = solver->T[solver->currentTes].vertexHandles[vn];
			const int& id =  v->info().id();
			CVector tf (0,0,0);
			int k=0;
			for (vector<const Real*>::iterator c = solver->perVertexPressure[id].begin(); c != solver->perVertexPressure[id].end(); c++)
				tf = tf + (*(solver->perVertexUnitForce[id][k++]))*(**c);
			v->info().forces = tf;
		}
		#endif*/
	}
	if (solver->debugOut) {
		CVector totalForce = nullVect;
		for (FiniteVerticesIterator v = Tri.finite_vertices_begin(); v != Tri.finite_vertices_end(); ++v) {
			if (!v->info().isFictious)
				totalForce = totalForce + v->info().forces;
			else if (solver->boundary(v->info().id()).flowCondition == 1)
				totalForce = totalForce + v->info().forces;
		}
		cout << "totalForce = " << totalForce << endl;
	}
}

bool TwoPhaseFlowEngine::detectBridge(RTriangulation::Finite_edges_iterator& edge)
{
	bool                            dryBridgeExist = true;
	const RTriangulation&           Tri            = solver->T[solver->currentTes].Triangulation();
	RTriangulation::Cell_circulator cell1          = Tri.incident_cells(*edge);
	RTriangulation::Cell_circulator cell0          = cell1++;
	if (cell0->info().saturation == 1) {
		dryBridgeExist = false;
		return dryBridgeExist;
	} else {
		while (cell1 != cell0) {
			if (cell1->info().saturation == 1) {
				dryBridgeExist = false;
				break;
			} else
				cell1++;
		}
		return dryBridgeExist;
	}
}

bool TwoPhaseFlowEngine::isCellNeighbor(unsigned int cell1, unsigned int cell2)
{
	bool neighbor = false;
	for (unsigned int i = 0; i < 4; i++) {
		if (solver->T[solver->currentTes].cellHandles[cell1]->neighbor(i)->info().id == cell2) {
			neighbor = true;
			break;
		}
	}
	return neighbor;
}

void TwoPhaseFlowEngine::setPoreThroatRadius(unsigned int cell1, unsigned int cell2, Real radius)
{
	if (isCellNeighbor(cell1, cell2) == false) {
		cout << "cell1 and cell2 are not neighbors." << endl;
	} else {
		for (unsigned int i = 0; i < 4; i++) {
			if (solver->T[solver->currentTes].cellHandles[cell1]->neighbor(i)->info().id == cell2)
				solver->T[solver->currentTes].cellHandles[cell1]->info().poreThroatRadius[i] = radius;
			if (solver->T[solver->currentTes].cellHandles[cell2]->neighbor(i)->info().id == cell1)
				solver->T[solver->currentTes].cellHandles[cell2]->info().poreThroatRadius[i] = radius;
		}
	}
}
Real TwoPhaseFlowEngine::getPoreThroatRadius(unsigned int cell1, unsigned int cell2)
{
	Real r = -1.;
	if (isCellNeighbor(cell1, cell2) == false) {
		cerr << "cell1 and cell2 are not neighbors." << endl;
	} else {
		for (unsigned int i = 0; i < 4; i++) {
			if (solver->T[solver->currentTes].cellHandles[cell1]->neighbor(i)->info().id == cell2) {
				r = solver->T[solver->currentTes].cellHandles[cell1]->info().poreThroatRadius[i];
				break;
			}
		}
	}
	return r;
}

} // namespace yade

#endif //TwoPhaseFLOW
