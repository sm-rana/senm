/* SenM - Sensor-driven hydraulic Model of water distribution systems
Copyright (C) 2013  Jinduan Chen <jinduan.uc@gmail.com>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License v2
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License v2 (http://www.gnu.org/licenses/gpl-2.0.html)
for more details.   */

#pragma once
#include "Network.h"

/** EPANET hydraulic solver 
Thread-safe EPANET-style hydraulic solver. Instead of using global variables,
the class maintains its own set of intermediate variables and result variables

*/
class Solver
{
public:  //type definitions
	typedef enum Err {
		OK, 
		NETWORK_NOT_EXIST,
		MALLOC_ERROR
	} ErrorCode;

	typedef enum _status {
		UN_ALLOCATED,
		ALLOCATED,
		SOLVED,
		ERR_SOLVER
	} status;

private:   //members
	int N, M;  //Node counts, link counts
	double Htol; /* head tolerance */ 
	double Hacc; /*hydraulic accuracy*/
	int MaxIter;  /*  max. hydraulic trials  */
	int ExtraIter;              /* Stop if network unbalanced     */
	double RQtol;  /* hydraulics parameters  */
	int CheckFreq, MaxCheck, DampLimit;

	static const double MISSING, QZERO;

	// the network
	Network* _net;

	double Hexp;    /* Exponent in headloss formula */

	/*
	** NOTE: Hydraulic analysis of the pipe network at a given point in time
	**       is done by repeatedly solving a linearized version of the 
	**       equations for conservation of flow & energy:
	**
	**           A*H = F
	**
	**       where H = vector of heads (unknowns) at each node,
	**             F = vector of right-hand side coeffs.
	**             A = square matrix of coeffs.
	**       and both A and F are updated at each iteration until there is
	**       negligible change in pipe flows.
	**
	**       Each row (or column) of A corresponds to a junction in the pipe
	**       network. Each link (pipe, pump or valve) in the network has a
	**       non-zero entry in the row-column of A that corresponds to its
	**       end points. This results in A being symmetric and very sparse.
	**       The following arrays are used to efficiently manage this sparsity:
	*/

	double   *D,                    /* Node actual demand           */
			*E,                    /* Emitter flows                */
			*K,                    /* Link settings                */
			*Q,                    /* Link flows                   */
			*X;                    /* General purpose array        */


	double   *H;                    /* Node heads                   */
	char     *S;                    /* Link status                  */



	double   *Aii,        /* Diagonal coeffs. of A               */
		*Aij,        /* Non-zero, off-diagonal coeffs. of A */
		*F;          /* Right hand side coeffs.             */
	double   *P,          /* Inverse headloss derivatives        */
		*Y;          /* Flow correction factors             */

public:    //api
	// Create a solver to compute against a network, allocate memory if needed
	static ErrorCode createSolver(Network* net, Solver** outsolver);

	// (Re)set control info collected by a snapshot
	// If no snapshots are provided, use 
	// default control (in the solver's network)
	ErrorCode setSnapshot(double* snapshot);

	// Solve hydraulic equations use a set of demands
	ErrorCode solveHyd(double* demands);

	// compute log likelihood of hydraulic measurements given 
	// the simulation results. If Snapshot is not provided, use the previously
	// set snapshot, if snapshot is not previously set, use default control
	// in the solver's network. network must be solved already
	double logLHyd(double* snapshot);

	// Get solver status
	status getStatus();

	// Query results
	double getPressure(const char* node_id);
	double getFlowrate(const char* link_id);

	// free memory
	~Solver();


protected:
	Solver(); //no public constructor

	void initlinkflow(int link, char status, double setting);

};

