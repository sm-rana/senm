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
/* The SenM class Network 
Singleton, read-only data container of an EPANET network's static data, including

[for hydraulic solver]
network topology; 
pipe list with roughness, diameter, and length info; 
junction list with elevation; 
pump, reservoir and tank list;
time and unit settings

[for initialization of demand model]
demand patterns and baseline demands

*/
#include <winnt.h>
#include <windows.h>
#include <map>
#include <tchar.h>
#include "vtkincludes.h"

class Network {
public:
	enum ErrorCode {
		OK, 
		NET_NOT_CREATED,
		CANT_OPEN_FILE,
		PATTERN_ID_EXISTS, 
		CURVE_ID_EXISTS,
		NODE_ID_EXISTS,
		LINK_ID_EXISTS,
		MALLOC_ERROR,
		NOT_ENOUGH_NODES,
		NO_TANKS,
		INPUT_LINE_TOO_LONG,
		SYNTAX_ERR_OF_SECT_NAME,
		SYNTAX_ERR_OF_JUNC,
		SYNTAX_ERR_OF_TANK_OR_RESERV,
		SYNTAX_ERR_OF_PIPE,
		SYNTAX_ERR_OF_PUMP,
		SYNTAX_ERR_OF_VALVE,
		SYNTAX_ERR_OF_PATTERN,
		SYNTAX_ERR_OF_CURVE,
		SYNTAX_ERR_OF_DEMAND,
		SYNTAX_ERR_OF_STATUS,
		SYNTAX_ERR_OF_OPTIONS,
		SYNTAX_ERR_OF_COORDINATES,
		ILLEGAL_NUMBER,
		DEMAND_PATTERN_UNDEFINED,
		GRADE_PATTERN_UNDEFINED,
		PATTERN_UNDEFINED,
		PUMP_CURVE_UNDEFINED,
		PUMP_PATTERN_UNDEFINED,
		VOLUME_CURVE_ERROR,
		VOLUME_CURVE_UNDEFINED,
		CURVE_UNDEFINED,
		NODE_UNDEFINED,
		NODE_COORDINATES_UNDEFINED,
		LINK_UNDEFINED,
		LINK_IS_A_LOOP,
		NODE_NOT_LINKED,
		CV_STATUS_ERROR,
		CURVE_HAS_NO_DATA,
		CURVE_DATA_NOT_ASC,
		PUMP_CURVE_INVALID,
		TANK_LEVEL_INVALID
	} ;

	  enum FieldType { /* Network variables:                  */
	//this enum is also used by Channel definition (see datasource.h),  
	//the database schema (see .sql files), and graphic objects (see sensor.h)
	//to denote the type of the physical varaible of a Scada system sensor
		  
		  //-----NODE VARIABLES -----//
		  //------------- sensor type 'L' -------------//
		  ELEV,         /*   nodal elevation                   */
						// also for tank levels, sensor type 'L', dim: Length
		  
						
		  DEMAND,       /*   nodal demand flow                 */
		  HEAD,         /*   nodal hydraulic head              */
						// total head, dim: Length
		  
		  //------------- sensor type 'P' -------------//			
		  PRESSURE,     /*   nodal pressure                    */
						//  free head, sensor type 'P', 
						// dim: Mass/(length*time^2) for hydraulic solver, 
						//      Length for sensor / scada db
		 
		  //-----LINK VARIABLES -----//
		  LENGTH,       /*   link length                       */
		  DIAM,         /*   link diameter                     */

		  //------------- sensor type 'Q' -------------//
		  FLOW,         /*   link flow rate                    */
						//  pipe flow rate, sensor type 'Q',
						//  dim: Length^3/time


		  VELOCITY,     /*   link flow velocity                */
		  HEADLOSS,     /*   link head loss                    */

		  //------------- sensor type 'C' -------------//
		  LSTATUS,      /*   link status                       */
						// on/off status for links, sensor type 'C', no dim

		  SETTING,      /*   pump/valve setting                */
		  FRICTION,     /*   link friction factor              */

		  //-----OTHER VARIABLES -----//
		  POWER,        /*   pump power output                 */
		  TIME,         /*   simulation time                   */
		  VOLUME,       /*   tank volume                       */
		  CLOCKTIME,    /*   simulation time of day            */
		  FILLTIME,     /*   time to fill a tank               */
		  DRAINTIME};   /*   time to drain a tank              */

	// get the Network instance, if no one exists,  created it.
	// If a network has been created, return it, in this case, inpfilename can be NULL.
	static ErrorCode getNetwork(LPCTSTR inpfilename, Network **net);

	// build a 2d and a 3d vtkNetwork polydata, 
	// net2d and net3d must have been allocated,
	// the scalars of the Polydatasets (attribute) at nodes (vertex-typed cells) 
	// are the baseline water demands
	ErrorCode get2d3dNet(vtkPolyData* net2d, vtkPolyData* net3d);

	// lookup the net2d/3d graphic cell id from network index (1-order)
	// the cellid < number of vertex cells (=MaxNodes) denotes nodes
	// the cellid > number of vertex cells (=MaxNodes) denotes links
	vtkIdType index2CellId(int index, FieldType type) ;

	// look up the graphic cell id from network component id (string)
	vtkIdType netId2CellId(char* id, FieldType type) ;

	// look up the node id (index but 0-order) from a water user id (0-order)
	// the 0-order node index is equivalent to the graphic cell id
	int userId2NodeId(vtkIdType userid) {
		if (userid >=0 && userid < Nusers)
			return _tabUser[userid];
		else return -1;
	};

	// number of nodes
	int getNnodes() {return MaxNodes;};
	int getNusers() {return Nusers;};

	// look up the node id (string) from a node index (1-order)
	void nodeId2NodeName(int nodeid, char** name) {
		if (nodeid>0 && nodeid<=MaxNodes)
			*name = Node[nodeid].ID;
	};

	// get the max baseline demand
	double getMaxBaselineDemand() {
		_maxBD = -DBL_MAX;
		for (int i=1;i<MaxNodes;++i)
			if (Node[i].D->Base > _maxBD)
				_maxBD = Node[i].D->Base;
		return _maxBD;
	}
	
	// get the value of baseline demand value at a certain percentile (0 to 1)
	double getBaseDemandPercentile(double pc);

	LPTSTR report() {//report network information in a string
		if (singleton) {
			_stprintf(_info, TEXT(
			"Network Info:\n\tTotal Nodes: %d, Tanks+Reservoirs: %d, Junctions: %d\n\tTotal Links: %d, Pipes: %d, Pumps: %d, Valves: %d\n"),
			MaxNodes, MaxTanks, MaxJuncs, MaxLinks, MaxPipes, MaxPumps, MaxValves); 
		} else {
			_stprintf(_info, TEXT("Network not loaded.\n"));
		}
		return _info;
	};

protected:
	TCHAR _info[512];
	int* _tabUser; //lookup table, user id to node id (both start with 0)
	vtkIdType *_tabNodeIndex, *_tabLinkIndex; //look up table for net2d/3d cellid
	double _maxBD; //max baseline demand

public:


	// Free memory
	~Network();


	// Maximum characters allowed in one line in the input file name
	static const int MAX_LINE = 255;

	// Maximum character allowed in the ID name
	static const int MAX_ID = 31;

	// Maximum tokens allowed in a line
	static const int MAX_TOKS = 40;

	// CONSTANTS
	static const double _PI, _TINY, _BIG;  //init in constructor

	// definitions of data structions for the network

	/* Element of list of floats */
	typedef struct  Floatlist { 
		double  value;
		struct  Floatlist *next;
	} SFloatlist;

	/* Element of temp list for Pattern & Curve data */
	typedef struct  Tmplist  {
		int        i;
		char       ID[MAX_ID+1];
		SFloatlist *x;
		SFloatlist *y;
		struct     Tmplist  *next;
	} STmplist;

	typedef struct        /* TIME PATTERN OBJECT */
	{
		char   ID[MAX_ID+1]; /* Pattern ID       */
		int    Length;      /* Pattern length   */
		double *F;          /* Pattern factors  */
	}  Spattern;

	typedef struct        /* CURVE OBJECT */
	{
		char   ID[MAX_ID+1]; /* Curve ID         */
		int    Type;        /* Curve type       */
		int    Npts;        /* Number of points */
		double *X;          /* X-values         */
		double *Y;          /* Y-values         */
	}  Scurve;

	struct Sdemand            /* DEMAND CATEGORY OBJECT */
	{
		double Base;            /* Baseline demand  */
		int    Pat;             /* Pattern index    */
		struct Sdemand *next;   /* Next record      */
	};
	typedef struct Sdemand *Pdemand; /* Pointer to demand object */

	typedef struct            /* NODE OBJECT */
	{
		char    ID[MAX_ID+1];    /* Node ID          */
		double  El;             /* Elevation        */
		double  x;				/* x-y coordinates of the node*/
		double  y;
		char    coordFlag;		/* the coordinates of the node is specified*/
		Pdemand D;              /* Demand pointer   */
	}  Snode;

	typedef struct            /* LINK OBJECT */
	{
		char    ID[MAX_ID+1];    /* Link ID           */
		int     N1;             /* Start node index  */
		int     N2;             /* End node index    */
		double  Diam;           /* Diameter          */
		double  Len;            /* Length            */
		double  Kc;             /* Roughness         */
		double  Km;             /* Minor loss coeff. */
		double  Kb;             /* Bulk react. coeff */
		double  Kw;             /* Wall react. coeff */
		double  R;              /* Flow resistance   */
		char    Type;           /* Link type         */
		char    Stat;           /* Initial status    */

	}  Slink;

	typedef struct     /* TANK OBJECT */
	{
		int    Node;     /* Node index of tank       */
		double A;        /* Tank area                */
   double Hmin;     /* Minimum water elev       */
   double Hmax;     /* Maximum water elev       */
   double H0;       /* Initial water elev       */
   double Vmin;     /* Minimum volume           */
   double Vmax;     /* Maximum volume           */
   double V0;       /* Initial volume           */
   double Kb;       /* Reaction coeff. (1/days) */
   double V;        /* Tank volume              */
   double C;        /* Concentration            */
   int    Pat;      /* Fixed grade time pattern */
   int    Vcurve;   /* Vol.- elev. curve index  */

	}  Stank;

	typedef struct     /* PUMP OBJECT */
	{
		int    Link;     /* Link index of pump          */
		int    Ptype;    /* Pump curve type             */
		/* (see PumpType below)        */
		double Q0;       /* Initial flow                */
		double Qmax;     /* Maximum flow                */
		double Hmax;     /* Maximum head                */
		double H0;       /* Shutoff head                */
		double R;        /* Flow coeffic.               */
		double N;        /* Flow exponent               */
		int    Hcurve;   /* Head v. flow curve index    */
		int    Ecurve;   /* Effic. v. flow curve index  */
		int    Upat;     /* Utilization pattern index   */
		int    Epat;     /* Energy cost pattern index   */
		double Ecost;    /* Unit energy cost            */
		double Energy[6];  /* Energy usage statistics:  */
		/* 0 = pump utilization      */
		/* 1 = avg. efficiency       */
		/* 2 = avg. kW/flow          */
		/* 3 = avg. kwatts           */
		/* 4 = peak kwatts           */
		/* 5 = cost/day              */
	}  Spump;

	typedef struct     /* VALVE OBJECT */
	{
		int   Link;     /* Link index of valve */
	}  Svalve;

	
	struct   Sadjlist         /* NODE ADJACENCY LIST ITEM */
	{
		int    node;            /* Index of connecting node */
		int    link;            /* Index of connecting link */
		struct Sadjlist *next;  /* Next item in list        */
	};
	/* Pointer to adjacency list item */
	typedef struct Sadjlist *Padjlist; 

	typedef enum _SectNum  
	
	{TITLE,JUNCTIONS,RESERVOIRS,TANKS,PIPES,PUMPS,
                  VALVES,CONTROLS,RULES,DEMANDS,SOURCES,EMITTERS,
                  PATTERNS,CURVES,QUALITY,STATUS,ROUGHNESS,ENERGY,
                  REACTIONS,MIXING,REPORT,TIMES,OPTIONS,
                  COORDS,VERTICES,LABELS,BACKDROP,TAGS,END

	} SectNum;
	static char* SectTxt[];  //defined in network.cpp, any changes made for Section
	//num much also change sect texts.

	typedef enum _Units {US, SI} UnitsType;

	enum FlowUnitsType             /* Flow units:                         */
	{CFS,          /*   cubic feet per second             */
	GPM,          /*   gallons per minute                */
	MGD,          /*   million gallons per day           */
	IMGD,         /*   imperial million gal. per day     */
	AFD,          /*   acre-feet per day                 */
	LPS,          /*   liters per second                 */
	LPM,          /*   liters per minute                 */
	MLD,          /*   megaliters per day                */
	CMH,          /*   cubic meters per hour             */
	CMD};         /*   cubic meters per day              */

	enum PressUnitsType            /* Pressure units:                     */
	{PSI,          /*   pounds per square inch            */
	KPA,          /*   kiloPascals                       */
	METERS};      /*   meters                            */

	enum FormType                  /* Head loss formula:                  */
	{HW,           /*   Hazen-Williams                    */
	DW,           /*   Darcy-Weisbach                    */
	CM};          /*   Chezy-Manning                     */

	enum NodeType                  /* Type of node:                       */
                {JUNC,          /*    junction                         */
                 RESERV,        /*    reservoir                        */
                 TANK};         /*    tank                             */

	enum LinkType                  /* Type of link:                       */
                 {CV,           /*    pipe with check valve            */
                  PIPE,         /*    regular pipe                     */
				  PUMP,			/*    pump  */
				  VALVE};                                    
	enum StatType                  /* Link/Tank status:                   */
                 {CLOSED,       /*   closed                            */
				  OPEN};         /*   open                              */
	enum PumpType                  /* Type of pump curve:                 */
                {CONST_HP,      /*    constant horsepower              */
                 POWER_FUNC,    /*    power function                   */
                 CUSTOM,        /*    user-defined custom curve        */
                 NOCURVE};



private:
	// The only instance allowed for this class
	static Network* singleton; 

	// Make Solver class a friend class, so all private members are accessible
	friend class Solver;

	// Options/configurations
	UnitsType Unitsflag;
	FormType	Formflag;
	FlowUnitsType Flowflag;
	PressUnitsType Pressflag;
	
    

    int DefPat;               /* Default demand pattern index   */
	char DefPatID[MAX_ID];	 /*Default pattern id */
    int Dmult;             /* Demand multiplier              */ 

	// Component counts
	int MaxJuncs, MaxTanks, MaxPipes, MaxPumps, MaxValves;
	int MaxNodes, MaxLinks;

	// Component iterators (for input processing only)
	int Njuncs, Nnodes, Ntanks, Npipes, Npumps, Nvalves, Nlinks;
	int Nusers;  //water users (base demand >0)
	int Ntokens;  /* Number of tokens in input line    */
	STmplist *PrevPat, *PrevCurve;

	// Components
	Snode    *Node;                 /* Node data                    */
	Slink    *Link;                 /* Link data                    */
	Stank    *Tank;                 /* Tank data                    */
	Spump    *Pump;                 /* Pump data                    */
	Svalve   *Valve;                /* Valve data                   */

	Padjlist *Adjlist;              /* Node adjacency lists         */	
	int *Degree;    /* Number of links adjacent to each node, used by the solver  */
	int Ncoeffs;               /* Number of non-0 matrix coeffs*/
	int     *Ndx;        /* Index of link's coeff. in Aij       */
	int     *Order;      /* Node-to-row of A                    */
	int 	*Row;        /* Row-to-node of A                    */
	
	/** The following arrays store the positions of the non-zero coeffs.    
	** of the lower triangular portion of A whose values are stored in Aij:
	*/
	int      *XLNZ,       /* Start position of each column in NZSUB  */
		*NZSUB,      /* Row index of each coeff. in each column */
		*LNZ;        /* Position of each coeff. in Aij array    */

	// Pump curves
	int MaxCurves;			// number
	STmplist *Curvelist;    //temporary list 
	Scurve *Curve; 

	// (demand) patterns
	int MaxPats;
	STmplist *Patlist;		//temporary list
	Spattern *Pattern;

	
	// Maps for id-index lookup (nodes and links), original EPANET hashtable is deprecated
	
	struct cmp_str 	{
		bool operator()(char const *a, char const *b)
		{
			return std::strcmp(a, b) < 0;
		}
	};
	std::map<char*, int, cmp_str> Nht, Lht;
	typedef std::pair<char*, int> HTPair;
	typedef std::map<char*, int>::iterator HTIt;

private:
	// default constructor, not accessible publicly
	Network(); 
	//members
	TCHAR _inpfilename[MAX_PATH];

	// Load data from INP file
	ErrorCode _loadInp(FILE*);

	// process a line and store characters to Tok[]
	ErrorCode coordata();
	ErrorCode juncdata();
	ErrorCode tankdata();
		ErrorCode inittanks();

	ErrorCode pipedata();

	ErrorCode pumpdata();  
		ErrorCode getpumpcurve(int, double*);
		int powercurve(double h0, double h1, double h2, double q1, double q2, 
			double *a, double* b, double* c);
		ErrorCode getpumpparams();

	ErrorCode valvedata();
	ErrorCode unlinked();

	// hydraulic solver preparation
	ErrorCode buildlists(int); // build adjacent lists
		int paralink(int node1, int node2, int link); //parallel link check
		void xparalinks();   // remove parallel links
		void countdegree();
		ErrorCode reordernodes();
		int mindegree(int, int);
		int growlist(int);
		int newlink(Padjlist);
		int linked(int, int);
		int addlink(int, int, int);
	ErrorCode storesparse(int);
	void freelists();
	ErrorCode ordersparse(int);
		void transpose(int n, int* il, int *jl, int *xl, 
			int* ilt, int* jlt, int* xlt, int* nzt);

	ErrorCode patterndata();  ErrorCode getpatterns();
	ErrorCode curvedata();  ErrorCode getcurves();
	ErrorCode demanddata();
	ErrorCode statusdata();
		void changestatus(int index, char status, double setting);
	ErrorCode optiondata();
		ErrorCode optionchoice(int);

	// unit conversion
	double Ucf[18]; //see FieldType
	void initunits();
	void convertunits();

	// compute link major headloss
	void resistence(int);
	double Hexp; /* Exponent in headloss formula */
	

	// Create a new pattern
	ErrorCode addpattern(char*);

	ErrorCode addcurve(char*);

	// misc post loading stuff
	void adjustdata();

	//helper functions:
	//interpolation using a curve
	double interp(int n, double* x, double* y, double);

	// find an item in a templist by id
	STmplist* findID(char*, STmplist* list);

	// determines which keyword appears on input line
	int  findmatch(char *line, char *keyword[]);

	// sees if substr matches any part of str, (case insensitive)
	int  match(char *str, char *substr);

	//scans string for tokens, saving pointers to them
    // in module global variable Tok[]
	int  gettokens(char*);
	char* Tok[MAX_TOKS];

	//get a double
	int getfloat(char*, double*);

	// free temp list
	void freeTmplist(STmplist *t);
	void freeFloatlist(SFloatlist *f);


};

