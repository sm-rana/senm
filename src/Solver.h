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
#include "DataSource.h"

/** EPANET hydraulic solver 
Thread-safe EPANET-style hydraulic solver. Instead of using global variables,
the class maintains its own set of intermediate variables and result variables
so that a network with multiple demand (loading) scenarios can be solved in many
threads at the same time.
*/
/* Instances of the Solver class is able to read the Network instance, 
however, to ensure other Solver's hydraulic simulation is not 
interfered, don't change the Network's internal data*/
struct Solver
{
// Error and warning facilities
	/// Error Codes
	enum EWICode {
		OK, 

		NETWORK_NOT_EXIST,
		MALLOC_ERROR,
		NETWORK_NOT_READY,
		DS_NOT_EXIST,

		NO_L_TANK,
		NO_F_PUMP_GPV,
		NO_B_PUMP_GPV,
		NO_V_TCV_FCV,
		CHAN_MINDEX_ERR,

        NO_CHANNELS,
		DUP_CHANNELS,
		D_CHANNEL_AT_TANK,
		FCV_REPLACED_BY_TCV,

		NOT_ENOUGH_XD,
		TOO_MANY_XD,
        B_EMITTER_CONFLICT,
        EQN_ILL_COND,

        VALVE_CAUSE_ILL_COND,
        UNABLE_TO_SOLVE,
        MAX_ITER_REACHED,
        CV_PSV_PRV_PROB,

		CHAN_DATA_NOT_MATCH,
		LAST_DUMMY_EWI
	};

    /// > reporting error, warning, and information messages
    static void reportEWI(EWICode err, unsigned thdN = 0);
	static void reportEWI(EWICode err, Network::ComponentType ctype, int id1, int id2, unsigned thdN = 0);
	static void reportEWI(EWICode err, Network::ComponentType ctype, int id1, unsigned thdN = 0);

// creation and destruction
	///> Factory method.
	/** Create a hydraulic solver for a given network-data source pair,
	\param[in] net    Network instance representing infrastructure
	\param[in] ds    Data source representing SCADA database for control and monitoring data
    \param[in] thdN   thread id 
	\param[out] outsolver    Pointer to where the pointer to the new Solver instance is stored.
	\return   Error information in creating the solver
	*/
	static EWICode createSolver(Network* net, DataSource* ds, unsigned thdN, Solver** outsolver);

	// free memory
	~Solver();

//solver functionalities

	/// Solve the hydraulic equations 
	/** \param [in]  xd  pointer to an array of stochastic demands. 
	           The array must be allocated and filled 
	   \param [in]   nXd   number of water demands. favorably ==_nXd;
	   \return   Number of run-time warnings issued 
	   \sa report()
   */
	void run(double* xd, int nXd, double* chd, int n_chd);

	/// solve the hydraulics using log-demand
	void runlogd(double *lxd, int nxd, double* chd, int n_chd);

	// compute log likelihood of hydraulic measurements given 
	// the simulation results. network must be solved already
	EWICode logL(double* ch_data, int n_ch, double* ll_out, bool disp_res=false);

	//> Get solver status
	//status getStatus();

	// Query results
	//double getPressure(const char* node_id);
	//double getFlowrate(const char* link_id);



	/// the network
	Network* _net;

	/// the datasource registered
	DataSource* _ds;


    unsigned threadN;  ///> which thread the solver is running in?
	int N, M;  //Node counts, link counts


	// Channel availability table (CAT), 1 ~ N: node; N+1 ~ N+M: link
	// using bit-wise operation to determine what channels are 
	// availabel for an asset (LFVBAQCPD)
	Channel::Type  *_tabCAT; 

	// number of channels
	int _nChan;

	// list of channels
	Channel* _lsChan; 

	//double* _ss; //Channel data snapshot pointer(data stored in DataSource class)

	// stochastic/variable water demand (xd)

	// number of stochastic demands
	int  _nXd; 

	//pointer to the water demands array (data stored in Population class)
	double * _xd;  

	// stochastic (i.e. variable) water demand - node index lookup table
	// size = MaxJuncs (number of junctions) - number of D channels
	// however, only the junctions with baseline demand > 0 will have an
	// entry here.
	int * _tabXd;  

	// junction index to xd index. If no xd for this junction, index = -1
	//size = MaxJuncs+1
	int * _tabIdx2Xd;  

	



	// hydarulic engine configuration
	double Htol; /* head tolerance */ 
	double Qtol; /* flow rate tolerance */ 
	double Hacc; /*hydraulic accuracy*/
	int MaxIter;  /*  max. hydraulic trials  */
	int ExtraIter;              /* Stop if network unbalanced     */
	double RQtol;  /* hydraulics parameters  */
//    double RelaxFactor;  // Relaxation factor for the iterative solver (flow)
//	double DampLimit; // if residual error is less than it, will decrease Ralaxation

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

	// Hydraulic engine internals
	double   *D,                    /* Node actual demand           */
		*E,                    /* Emitter flows                */
		*K,                    /* Link settings                */
		*Q,                    /* Link flows                   */
		*X;                    /* General purpose array        */

	double   *H;                    /* Node heads                   */

public:
    /// if a link is enabled, status of it
    enum LinkStatus {
        FULL_OPEN, /// opened no valving headloss
        PARTIAL, /// partially opened
        CLOSED, ///> closed, but may be opened in subsequent iterations
        DISABLED, ///> keep closed
	};
	LinkStatus     *S;                    /* Link status                  */

	// Sparse matrix solver, XLNZ, NZSUB, and LNZ are fixed and thus stored in Network class
	double   *Aii,        /* Diagonal coeffs. of A               */
		*Aij,        /* Non-zero, off-diagonal coeffs. of A */
		*F;          /* Right hand side coeffs.             */
	double   *P,          /* Inverse headloss derivatives        */
		*Y;          /* Flow correction factors             */

    double  *B;  ///> B channel data (discharge-side pressure)


	void initlinkflow(int link, char status, double setting);

	/// set channel data for Channels [L][F][B][V][C][D]
	EWICode setLFBVCD(double* chd, int n_ch);

    /// renew coefficients at the beginning of each iteration
    void newcoeffs();

    /// @{
    /// renew coefficients for each type of component
    void linkcoeffs();
    void pipecoeff(int k);

    void pumpcoeffB(int k); // SenM specific

    void tcvcoeff(int k);
    void prvcoeff(int k, int n1, int n2);
    void psvcoeff(int k, int n1, int n2);

    void valvecoeff(int k);

    double Solver::DWcoeff(int k, double *dfdq);
    void nodecoeffs();
    void emittercoeffs();

    void valvecoeffs();

    /// @}
//	double SGN(double x) {return (x<0?-1:1);};
//	double SQR(double x) {return x*x;};
    
    ///> iterative solver
    int linsolve(int n ,double* Aii, double* Aij, double *B);

    int badvalve(int n);
    double newflows();

    double emitflowchange(int);

    int valvestatus();
    int linkstatus();
    LinkStatus cvstatus(LinkStatus s, double dh, double q);
       

    LinkStatus  psvstatus(int k, LinkStatus  s, double hset, double h1, double h2);
    LinkStatus  prvstatus(int k, LinkStatus  s, double hset, double h1, double h2);




protected:
	Solver(); //no public constructor

};

