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
/* Instances of the Solver class is able to read private members of the singleton 
Network (friend class), however, to ensure other Solver's hydraulic simulation is not 
interfered, don't change the Network's internal data*/
class Solver
{
public:  // Error and warning facilities
	/// Error Codes
	enum ECode {
		OK, 
		NETWORK_NOT_EXIST,
		MALLOC_ERROR,
		NETWORK_NOT_READY,
		DS_NOT_EXIST,
		NO_L_TANK,
		NO_F_PUMP_GPV,
		NO_B_PUMP_GPV,
		NO_V_TCV_FCV,
		JUNC_CHAN_MINDEX_ERR,
		TANK_CHAN_MINDEX_ERR,
		VALVE_CHAN_MINDEX_ERR,
		PUMP_CHAN_MINDEX_ERR,
		PIPE_CHAN_MINDEX_ERR,

		LAST_DUMMY_E
	}; 
	/// Warning Codes
	enum WCode {
		NONE,
		NO_CHANNELS,
		DUP_CHANNELS,
		D_CHANNEL_AT_TANK,
		FCV_REPLACED_BY_TCV,
		NOT_ENOUGH_XD,
		TOO_MANY_XD,
        B_EMITTER_CONFLICT,
        EQN_ILL_COND,
        VALVE_CAUSE_ILL_COND,
        MAX_ITER_REACHED,
        CV_PSV_PRV_PROB,

		LAST_DUMMY_W
	};
	/// Warning record
	struct Warn {
		WCode wc; ///> warn code
		Network::FieldType ctype; ///> component type
		int cindex; ///> component index (epanet)
		Warn(): wc(NONE), ctype(Network::UNKNOWN), cindex(0) {}; ///>default contructor
		///> trivial constructor
		Warn(WCode wc_in, Network::FieldType ctype_in, int cindex_in=0):
			wc(wc_in), ctype(ctype_in), cindex(cindex_in) {};
		///> constructor with only WCode specified
		Warn(WCode wc_in) : wc(wc_in), ctype(Network::UNKNOWN), cindex(0) {};
	};
	/// Error record
	struct Error {
		ECode ec; ///>Creation time error
		Network::FieldType ctype; ///>Creation time error component type
		int cindex;  ///Creation time error index (epanet)
		Error(): ec(OK), ctype(Network::UNKNOWN), cindex(0) {}; ///>default contructor
	};



	static Error errorC; ///>creation-time error

	/// @{
	/// Look up the text explanation of an error or a warning
	/** \param[in]   err/warn   Error or warning code.
	\param[out]   txt      A tchar string pre-allocated to be filled with info
	*/
	static void fillEWinfo(Error err, TCHAR* txt);
	static void fillEWinfo(Warn warn, TCHAR* txt);
	/// @}

protected:
    ///Creation-time warning list, the list is built during solver creation
    /// and never emptied
	std::list<Warn> _warnListC; 
    ///>Run-time warning list, the list is re-built with each hydraulic simulation
	std::list<Warn> _warnListR; 

public: // creation and destruction
	///> Factory method.
	/** Create a hydraulic solver for a given network-data source pair,
	\param[in] net    Network instance representing infrastructure
	\param[in] ds    Data source representing SCADA database for control and monitoring data
	\param[out] outsolver    Pointer to where the point to the new Solver instance is stored.
	\return   Error information in creating the solver
	*/
	static ECode createSolver(Network* net, DataSource* ds, Solver** outsolver);

	/// Report solver creation error
	static void reportCreationError(); 

	// free memory
	~Solver();

public:    //solver functionalities

	/// Solve the hydraulic equations 
	/** \param [in]  xd  pointer to an array of stochastic demands. 
	           The array must be allocated and filled 
	   \param [in]   nXd   number of water demands. favorably ==_nXd;
	   \return   Number of run-time warnings issued 
	   \sa report()
   */
	int run(double* xd, int nXd);

	/// report all warnings in the both warning list to stderr. 
	void report();


	// compute log likelihood of hydraulic measurements given 
	// the simulation results. If Snapshot is not provided, use the previously
	// set snapshot, if snapshot is not previously set, use default control
	// in the solver's network. network must be solved already
	//double logLHyd(double* snapshot);

	//> Get solver status
	//status getStatus();

	// Query results
	//double getPressure(const char* node_id);
	//double getFlowrate(const char* link_id);



protected:   //members
	// the network
	Network* _net;

	// the datasource registered
	DataSource* _ds;


protected:
	int N, M;  //Node counts, link counts

	// Channel availability table, 1 ~ N: node; N+1 ~ N+M: link
	Channel::Type  *_tabCAT; 
	int _nChan;
	Channel* _lsChan; //channel list
	double* _ss; //Channel data snapshot pointer(data stored in DataSource class)

	// stochastic/variable water demand (xd)
	int  _nXd; // number of stochastic demands
	int * _tabXd;  // stochastic (i.e. variable) water demand - node index lookup table
	int * _tabIdx2Xd;  // junction index to xd index
	double * _xd;  //pointer to the water demands array (data stored in Population class)

	// hydarulic engine configuration
	double Htol; /* head tolerance */ 
	double Qtol; /* flow rate tolerance */ 
	double Hacc; /*hydraulic accuracy*/
	int MaxIter;  /*  max. hydraulic trials  */
	int ExtraIter;              /* Stop if network unbalanced     */
	double RQtol;  /* hydraulics parameters  */
//    double RelaxFactor;  // Relaxation factor for the iterative solver (flow)
//	double DampLimit; // if residual error is less than it, will decrease Ralaxation

	// constants
	//static const double MISSING, QZERO;
    static const double CBIG, CSMALL;
    // constants for DW hydraulic equation
    static const double A1, A2, A3, A4, A8, A9, AA, AB, AC;

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
protected:
	LinkStatus     *S;                    /* Link status                  */

	// Sparse matrix solver, XLNZ, NZSUB, and LNZ are fixed and thus stored in Network class
	double   *Aii,        /* Diagonal coeffs. of A               */
		*Aij,        /* Non-zero, off-diagonal coeffs. of A */
		*F;          /* Right hand side coeffs.             */
	double   *P,          /* Inverse headloss derivatives        */
		*Y;          /* Flow correction factors             */

    double  *B;  ///> B channel data (discharge-side pressure)

protected:
	Solver(); //no public constructor

	void initlinkflow(int link, char status, double setting);

	/// set channel data for Channels [L][F][B][V][C][D]
	void setLFBVCD();

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
	double SGN(double x) {return (x<0?-1:1);};
	double SQR(double x) {return x*x;};
    
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





};

