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

#include <map>
#include <list>
#include "SenmCoreIncs.h"

/// Max string size of a component type string
#define MAX_NET_CTYPE_STR_SIZE 64

/// EPANET network infrastructure
/**Singleton data container of distribution network infrastructure based on an EPANET inp file. 
Data managed by the Network class including:

[for hydraulic solver]
- network topology; 
- pipe list with roughness, diameter, and length info; 
- junction list with elevation; 
- pump, reservoir and tank list;
- time and unit settings

[for initialization of the demand model]
- demand patterns 
- baseline demands

Any component in the network have 2 different types of IDs. 
- (EPANET) index:  Index defined by the EPANET engine. Starts from 1. 
Node index: 1 <= index <= MaxJuncs: Junctions
MaxJuncs < index <= MaxNodes :  tanks(including reservoirs).
Link index: 1 <= index <= MaxPipes:  Pipes, 
MaxPipes < index <= MaxPumps: pumps, 
MaxPumps < index <= MaxLinks: valves
- (EPANET) ID:  ID string (i.e., name) defined by the EPANET engine. One-on-one map to index. 

The PBV cannot be modeled in SenM.
The FCV is treated as a TCV, (its setting is the minor head loss coeff.)
*/
struct Network {
	//error and warning system
	/** Errors produced when creating the network instance.*/
	enum ErrorCode  {
		OK, 
		UNLOADED,
		NET_NOT_CREATED,
		CANT_OPEN_FILE,
		PATTERN_ID_EXISTS, 
		CURVE_ID_EXISTS,
		NODE_ID_EXISTS,
		LINK_ID_EXISTS,
		MALLOC_ERROR,
		NO_JUNCS,
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
		SYNTAX_ERR_OF_EMITTERS,
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
		SET_CV_STATUS,
		CURVE_HAS_NO_DATA,
		CURVE_DATA_NOT_ASC,
		PUMP_CURVE_INVALID,
		TANK_LEVEL_INVALID,
		UNKNOWN_UNIT,
	} ;

	/// Warnings in creating the SenM network
	enum WarnCode {
		NONE,
		PBV_SPECIFIED_IN_INP,
		TCV_REPLACE_FCV,

	};     

	/// Return the pointer to the text of an error code stored in a member char array
	LPCTSTR problemText(ErrorCode ec); 
	/// Return the text of a warning message
	LPCTSTR problemText(WarnCode ec); 

	/// Report the list of warnings to stdout
	void reportWarns(); 

	ErrorCode _ecCur; // Current error code.
	int _lineEc; // Line in which an error is detected
	TCHAR _strProblem[1024]; // Text of error or warning messages.

	int _lineno;  // current line number
	typedef std::pair<int, WarnCode> Warn; //a warning instance = line no. + warncode
	std::list<Warn> _lsWC; // list of warning codes.   


	// Creation and destruction
	///> Factory method. 
	/** obtain a Network instance based on the EPANET inpfile, 
	The Network instance is a singleton. If a network has not been created, 
	create one from the input file. Otherwise, return the existing one. 
	In the latter case, inpfilename can be NULL.
	\param [in] inpfilename        Name of the input file. May include path
	\param [out] net               Pointer to where a pointer to the network is stored
	\return        Error code typed \sa Network::ErrorCode
	*/
	static ErrorCode getNetwork(LPCTSTR inpfilename, Network **net);

	///> Destroy the network singleton and free memory
	~Network();

	///> report network information to stdout
	static void report();

	// Publiclly accessible types and constants
	/** Variable types*/
	/** EPANET has its own typing system for the attributes of components. 
	This enum is also used by Channel definition (\sa datasource.h),  Network class
	the database schema (see .sql files), and graphic objects (\sa sensor.h)
	to specify the type of the physical varaible of a Scada system sensor */
	enum FieldType { 

		//-----NODE VARIABLES -----//
		///   nodal elevation, 
		/** The enum is also for monitored tank levels , 
		and channel type 'L', The dimension of an ELEV variable is Length
		\sa SensorTank
		*/
		ELEV,             

		DEMAND,       ///>   nodal demand flow                 
		HEAD,         ///> nodal hydraulic (total) head, dim: Length

		///   nodal pressure, 
		/**  In the Solver the item is for free pressure, and has the dimension
		of Mass/(length*time^2). 
		the item is also used for monitored pressure head, 
		there it has the dimension of Length
		\sa SensorHead
		*/
		PRESSURE,

		//-----LINK VARIABLES -----//
		LENGTH,       ///>   link length                       
		DIAM,         ///>   link diameter                     

		/// Pipe/pump/valve flow rate
		/** In both Solver and SensorFlow the item is for flow rate,
		with dimension of Length^3/time
		\sa SensorFlow
		*/
		FLOW,         /*   link flow rate                    */

		VELOCITY,     ///>   link flow velocity                
		HEADLOSS,     ///>   link head loss                    

		///   link status                       
		/// on/off status for links,  no dimension.
		///\sa SensorControl
		LSTATUS,      

		SETTING,      ///>   pump/valve setting                
		FRICTION,     ///>   link friction factor              

		//-----OTHER VARIABLES -----//
		POWER,        ///>   pump power output                 
		TIME,         ///>   simulation time                  
		VOLUME,       ///>   tank volume                     
		CLOCKTIME,    ///>   simulation time of day         
		FILLTIME,     ///>   time to fill a tank           
		DRAINTIME,   ///>   time to drain a tank         
		UNKNOWN,  ///> unkown type of variable
	};

	/** Head loss formula:                  */
	enum FormType
	{HW,           /**>   Hazen-Williams                    */
	DW,           /**>   Darcy-Weisbach                    */
	CM};          /**>   Chezy-Manning                     */

	/** Type of node:                       */
	enum ComponentType                  {
		JUNC,          /**>    junction                         */
		RESERV,        /**>    reservoir                        */
		TANK,       ///> Water storage facility
		CV,           /**>    pipe with check valve            */
		PIPE,         /**>    regular pipe                     */
		PUMP,			/**>    pump  */
		PRV,          /**>    pressure reducing valve          */
		PSV,          /**>    pressure sustaining valve        */
		PBV,          /**>    pressure breaker valve           */
		FCV,          /**>    flow control valve               */
		TCV,          /**>    throttle control valve           */
		GPV,         /**>    general purpose valve            */
        UNDETERMINED,  ///> not determined component or no component
	};

    /// fetch text description of a component
    static TCHAR* typeStr(ComponentType, TCHAR*);

	/** Link/Tank status:                   */
	enum StatType                  
	{CLOSED,       /**>   closed                          */
	OPEN,         /**>   open                              */
	//    ACTIVE,  ///>  
	};
	/** Type of pump curve:                 */
	enum PumpType                  
	{CONST_HP,      /**>    constant horsepower              */
	POWER_FUNC,    /**>    power function                   */
	CUSTOM,        /**>    user-defined custom curve        */
	NOCURVE};
	///> Maximum characters allowed in one line in the input file name
	static const int MAX_LINE = 255;

	///> Maximum character allowed in the ID name
	static const int MAX_ID = 31;

	///> Maximum tokens allowed in a line
	static const int MAX_TOKS = 40;

	//Unit system
	///> General units type
	/** Specifies type of units for variables other than flow rate 
	and pressure */
	enum UnitsType {US, SI}; 

	/** Flow units:                         */
	enum FlowUnitsType             
	{CFS,          /**>   cubic feet per second             */
	GPM,          /**>   gallons per minute                */
	MGD,          /**>   million gallons per day           */
	IMGD,         /**>   imperial million gal. per day     */
	AFD,          /**>   acre-feet per day                 */
	LPS,          /**>   liters per second                 */
	LPM,          /**>   liters per minute                 */
	MLD,          /**>   megaliters per day                */
	CMH,          /**>   cubic meters per hour             */
	CMD};         /**>   cubic meters per day              */

	/** Pressure units:                     */
	enum PressUnitsType 
	{PSI,          /**>   pounds per square inch            */
	KPA,          /**>   kiloPascals                       */
	METERS};      /**>   meters                            */

	/** Text representation abbreviation level*/
	enum AbbrLev
	{UNIT_SHORT, UNIT_FULL};

	/// @{
	///> Get string representation of a unit
	/** \param [in] type     Type of the variable
	\param [in] type     Unit of the variable
	\param [out] str     String/text representation of the unit.
	Storage of the char string should be been allocated
	\return  Errors
	*/
	static ErrorCode unitText(FieldType type, UnitsType unit, TCHAR* str, AbbrLev al=UNIT_SHORT); 
	static ErrorCode unitText(FieldType type, FlowUnitsType unit, TCHAR* str, AbbrLev al=UNIT_SHORT); 
	static ErrorCode unitText(FieldType type, PressUnitsType unit, TCHAR* str, AbbrLev al=UNIT_SHORT); 

	/// @}

	//Utilities

	/// Get number of nodes
	int getNnodes() {return MaxNodes;};
	/// Get number of water users.
	int getNusers() {return Nusers;};
	/// Get reference total system average demand
	float* getRefD() {return ref_demand;};

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

	///> Find the EPANET ID string for an EPANET index
	/** \param [in] nodeid        EPANET index
	\param [out] name         Pointer to where the EPANET id string is stored
	\param [in] name_size     Size of the name char array allocated
	if the index does not exist or name_size>MAX_ID, return empty string
	*/
	void nodeId2NodeName(int nodeid, char* name, int name_size) {
		if (nodeid>0 && nodeid<=MaxNodes)
			strcpy(name, Node[nodeid].ID);
		else name[0] = '\0'; // empty string
	};

    /// Find index for the EPANET ID
    int name2idx(char* name) {
        HTIt it;
		it = Nht.find(name);
		if (it!=Nht.end()) return it->second;

        it = Lht.find(name);
		if (it!=Lht.end()) return it->second;

        return 0; //not found
	}

	/// Get the maxiumu baseline demand
	double getMaxBaselineDemand() {
		_maxBD = -DBL_MAX;
		for (int i=1;i<MaxNodes;++i)
			if (Node[i].D->Base > _maxBD)
				_maxBD = Node[i].D->Base;
		return _maxBD;
	}

	/// Get a percentile of the baseline demand
	/** \param [in] pc      percentile, must between 0 and 1
	*/
	double getBaseDemandPercentile(double pc);



	// definitions of data structures for the network
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
		double  Ke;             /* Emitter coeff */
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

	//defined in network.cpp, any changes made to _SectNum
	//must also change SectTxt[].
	static char* SectTxt[];  


	TCHAR _info[1024]; ///> information string about the loaded network

	double _maxBD; //max baseline demand



	// The only instance allowed for this class
	static Network* singleton; 

	/// Make Solver class a friend class, so all private members are accessible there
	/// \sa Solver
	friend class Solver;

	/// Make VisNetwork class a friend class, 
	/// \sa VisNetwork
	friend class VisNetwork;

	// Options/configurations
	UnitsType Unitsflag;
	FormType	Formflag;
	FlowUnitsType Flowflag;
	PressUnitsType Pressflag;



	int DefPat;               /* Default demand pattern index   */
	char DefPatID[MAX_ID];	 /*Default pattern id */
	double Dmult;             /* Demand multiplier              */ 
	float Qexp;  /* Emitter flow exponent */

	// Component counts
	int MaxJuncs, MaxTanks, MaxPipes, MaxPumps, MaxValves;
	int MaxNodes, MaxLinks;
	float* ref_demand; // base hourly water demand for users

	// Component iterators (for input processing only)
	int Njuncs, Nnodes, Ntanks, Npipes, Npumps, Nvalves, Nlinks;
	int Nusers;  //water users (base demand >0)
	int Nemitters;  //count of junctions with emitters
	int Ntokens;  /* Number of tokens in input line    */
	STmplist *PrevPat, *PrevCurve;

	// Components
	Snode    *Node;                 /* Node data                    */
	Slink    *Link;                 /* Link data                    */
	Stank    *Tank;                 /* Tank data                    */
	Spump    *Pump;                 /* Pump data                    */
	Svalve   *Valve;                /* Valve data                   */

	// Quick indexing arrays, the listed items have different treatment
	// in EPANET and SenM
	int* tPumpGPV;    // pumps and GPVs
	int nPumpGPV;  // number of pumps and GPVs
	int* tCRSV;   // Check valves, PRV, PSV
	int nCRSV;     // number of check valves, PRV, PSV.
	int* tTFV;  // throttle control valves, flow control vavles
	int nTFV;   //number of TCV, FCV

	Padjlist *Adjlist;              /* Node adjacency lists         */	
	int *Degree;    /* Number of links adjacent to each node, used by the solver  */
	int Ncoeffs;               /* Number of non-0 matrix coeffs*/
	int     *Ndx;        /* Index of link's coeff. in Aij       */
	int     *Order;      /* which node corresponds to the [i]th row in A?           */
	int 	*Row;        /* which row in A corresponds to the [i]th node?           */

	/** Sparse matrix representation 
	The following arrays store the positions of the non-zero coeffs.    
	of the lower triangular portion of A whose values are stored in Aij:
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


	// disable default constructor
	Network(); 
	//members
	TCHAR _inpfilename[MAX_PATH];

	// Load data from INP file
	ErrorCode loadInp(FILE*);

	// process a line and store characters to Tok[]
	ErrorCode coordata();
	ErrorCode juncdata();
	ErrorCode tankdata();
	ErrorCode emitterdata();
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

	// unit conversion related
	double Ucf[18]; //see FieldType
	void initunits();
	void convertunits();

	// compute link's major headloss
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

