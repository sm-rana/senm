#include <cstring>
#include <cstdio>
#include <cmath>
#include <tchar.h>
#include <Windows.h>
#include "Percentile.h"
#include "Network.h"
#include "senmCoreIncs.h"

Network* Network::singleton = NULL;

//strings of input file's section names, must be consistent with
//the enum type SectNum
char* Network::SectTxt[] = {
	"[TITLE]","[JUNCTIONS]","[RESERVOIRS]","[TANKS]","[PIPES]","[PUMPS]",
	"[VALVES]","[CONTROLS]","[RULES]","[DEMANDS]","[SOURCES]","[EMITTERS]",
	"[PATTERNS]","[CURVES]","[QUALITY]","[STATUS]","[ROUGHNESS]","[ENERGY]",
	"[REACTIONS]","[MIXING]","[REPORT]","[TIMES]","[OPTIONS]",
	"[COORDINATES]","[VERTICES]","[LABELS]","[BACKDROP]","[TAGS]","[END]"
};


Network::Network() { //setdefaults() in input1.c

	Formflag = HW;              /* Use Hazen-Williams formula     */
	Unitsflag= US;				/* US unit system                 */
	Flowflag = GPM;             /* Flow units are gpm             */
	Pressflag = PSI;             /* Pressure units are psi         */

	DefPat    = 0;               /* Default demand pattern index   */
	strncpy(DefPatID,"1",MAX_ID); /* Default pattern id*/
	Dmult     = 1.0;             /* Demand multiplier              */ 
	Qexp = 2.0;  /* Emitter flow exponent */

	MaxJuncs = Njuncs = 0;
	Nusers = 0; Nemitters = 0;
	MaxTanks = Ntanks = 0;
	MaxNodes = Nnodes = 0;
	MaxPipes = Npipes = 0;
	MaxPumps = Npumps = 0;
	MaxValves = Nvalves = 0;
	MaxLinks = Nlinks = 0;
	MaxCurves= 0;

	//sizes of quick indexing arrays
	nPumpGPV = nCRSV = nTFV = 0;

	PrevPat = PrevCurve = NULL;

	Patlist = NULL;
	Curvelist = NULL;

	_ecCur = UNLOADED;  // current error code
	_lineEc = 0;
	_lineno = 0;

}


LPCTSTR Network::problemText(ErrorCode ec) { 

	switch (ec) {
	case OK:
		_tcscpy(_strProblem, TEXT("OK."));
		break;
	case NET_NOT_CREATED:
		_tcscpy(_strProblem, 
			TEXT("No network has been loaded. Please specify a network file."));
		break;
	case CANT_OPEN_FILE:
		_tcscpy(_strProblem,
			TEXT("Cannot open the network file."));
		break;
	case PATTERN_ID_EXISTS: case CURVE_ID_EXISTS: case NODE_ID_EXISTS:
	case LINK_ID_EXISTS:
		_tcscpy(_strProblem, TEXT("INP file refers to an ID not existed."));
		break;
	case MALLOC_ERROR:
		_tcscpy(_strProblem, TEXT("Memory allocation error when building components."));
		break;
	case NO_JUNCS:
		_tcscpy(_strProblem, TEXT("There are no junctions specified in the network."));
		break;
	case NO_TANKS:
		_tcscpy(_strProblem, TEXT("There are no tanks or reservoirs in the network."));
		break;
	case 		INPUT_LINE_TOO_LONG :
	case 		SYNTAX_ERR_OF_SECT_NAME :
	case 		SYNTAX_ERR_OF_JUNC :
	case 		SYNTAX_ERR_OF_TANK_OR_RESERV :
	case 		SYNTAX_ERR_OF_PIPE :
	case 		SYNTAX_ERR_OF_PUMP :
	case 		SYNTAX_ERR_OF_VALVE :
	case 		SYNTAX_ERR_OF_PATTERN :
	case 		SYNTAX_ERR_OF_CURVE :
	case 		SYNTAX_ERR_OF_DEMAND :
	case 		SYNTAX_ERR_OF_STATUS :
	case 		SYNTAX_ERR_OF_OPTIONS :
	case 		SYNTAX_ERR_OF_COORDINATES :
		_tcscpy(_strProblem, TEXT("INP file syntax error."));
		break;
	case  ILLEGAL_NUMBER:
		_tcscpy(_strProblem, TEXT("Number in the INP file is not in the right format."));
		break;
	case 		DEMAND_PATTERN_UNDEFINED: 
	case 		GRADE_PATTERN_UNDEFINED: 
	case 		PATTERN_UNDEFINED: 
	case 		PUMP_CURVE_UNDEFINED: 
	case 		PUMP_PATTERN_UNDEFINED: 
	case 		VOLUME_CURVE_ERROR: 
	case 		VOLUME_CURVE_UNDEFINED: 
	case 		CURVE_UNDEFINED: 
	case 		NODE_UNDEFINED: 
	case 		LINK_UNDEFINED: 
		_tcscpy(_strProblem, TEXT("Component not defined."));
		break;
	case 		NODE_COORDINATES_UNDEFINED: 
		_tcscpy(_strProblem, TEXT("Coordinates of the nodes are not defined. ")
			TEXT("This will not impact hydraulic simulation but 3d visuals can't be built.") );
		break;
	case  LINK_IS_A_LOOP:
		_tcscpy(_strProblem, TEXT("Some link is a loop - cannot build network."));
		break;
	case NODE_NOT_LINKED:
		_tcscpy(_strProblem, TEXT("Some node is not linked to the network."));
		break;
	case SET_CV_STATUS:
		_tcscpy(_strProblem, TEXT("Cannot set the on/off status of a check valve."));
		break;
	case CURVE_HAS_NO_DATA: case CURVE_DATA_NOT_ASC: 
	case PUMP_CURVE_INVALID:
		_tcscpy(_strProblem, TEXT("Curve data error."));
		break;
	case TANK_LEVEL_INVALID:
		_tcscpy(_strProblem, TEXT("Tank level is not in the min/max range."));
		break;
	default:
		_tcscpy(_strProblem, TEXT("Unknown error code supplied."));
	}
	return _strProblem;
}

LPCTSTR Network::problemText(WarnCode wc) {
	switch (wc) {

	case PBV_SPECIFIED_IN_INP:
		_tcscpy(_strProblem, TEXT("PBV cannot be modeled in SenM. it will")
			TEXT("be replaced with a CV"));
		break;
	default:
		_tcscpy(_strProblem, TEXT(""));
		break;
	}
	return _strProblem;
}



Network::ErrorCode Network::getNetwork(LPCTSTR inpfilename, Network **net) {

	// if no file name provided, return current singleton
	if (inpfilename == NULL || _tcscmp(inpfilename, TEXT("")) == 0)  {
		if (singleton == NULL) {
			*net = NULL;
			//		_ftprintf(stderr, TEXT("No network has been loaded. Please specify a network file.\n"));
			return NET_NOT_CREATED;
		} else {  // use existing network
			*net = singleton;
			return singleton->_ecCur;
		}
	} 

	// if file name provided and equal to the old singleton, return old one
	if (singleton != NULL && 
		_tcscmp(inpfilename, singleton->_inpfilename) == 0) {
			// use existing network
			*net = singleton;
			return singleton->_ecCur;
	} 

	//otherwise (file name provided, different from old one)
	//, create a new network
	FILE* infile;
	infile = _tfopen(inpfilename, TEXT("rt"));
	if (infile == NULL) {
		*net = NULL;
		//	_ftprintf(stderr, TEXT("Cannot open the network file.\n"));
		return CANT_OPEN_FILE;
	}

	delete singleton; // delete previous one
	singleton = new Network();
	_tcscpy_s(singleton->_inpfilename, inpfilename);

	singleton->loadInp(infile);
	fclose(infile);

	TCHAR tmptxt[MAX_ERR_STRING_SIZE];
    _stprintf(tmptxt, TEXT("Try loading the network in file (%s).\n%s \n"), 
		inpfilename, singleton->problemText(singleton->_ecCur));
    ewi(tmptxt);
	if (singleton->_lineEc >0)
		_ftprintf(stderr, TEXT("Possible error in Line %d\n"), 
		singleton->_lineEc);

	if (singleton->_ecCur != OK) {
		delete(singleton);
		*net = NULL;
		return singleton->_ecCur;
	} else {
		*net = singleton;
		return OK;
	}
}

Network::~Network()
{
	int j;
	Pdemand demand, nextdemand;

	/* Free memory for node data */
	if (Node != NULL)
	{
		for (j=0; j<=MaxNodes; j++)
		{
			/* Free memory used for demand category list */
			demand = Node[j].D;
			while (demand != NULL)
			{
				nextdemand = demand->next;
				free(demand);
				demand = nextdemand;
			}
		}
		free(Node);
	}

	/* Free memory for other network objects */
	free(Link);
	free(Tank);
	free(Pump);
	free(Valve);

	/* Free memory for time patterns */
	if (Pattern != NULL)
	{
		for (j=0; j<=MaxPats; j++) free(Pattern[j].F);
		free(Pattern);
	}

	/* Free memory for curves */
	if (Curve != NULL)
	{
		for (j=0; j<=MaxCurves; j++)
		{
			free(Curve[j].X);
			free(Curve[j].Y);
		}
		free(Curve);
	}

	free(Degree);
    free(ref_demand);

	/* free memory for node and link hash table*/
	for (HTIt it = Nht.begin(); it!=Nht.end(); ++it) {
		delete[] it->first;
	}

	for (HTIt it = Lht.begin(); it!=Lht.end(); ++it) {
		delete[] it->first;
	}

	// delete quick indexing arrays
	free(tPumpGPV);
	free(tCRSV);
	free(tTFV);
}

Network::ErrorCode Network::loadInp(FILE* InFile) {

	// use the code in EPANET toolkit
	// netsize() in input2.c, first scan, determine various sizes

	char  line[MAX_LINE+1];     /* Line from input data file    */
	char  *tok;                /* First token of line          */
	int   sect,newsect;        /* Input data sections          */

	_ecCur = UNLOADED;

	/* Initialize network component counts */
	sect        = -1;

	/* Add a default pattern 0 */
	MaxPats = -1;
	_ecCur = addpattern("");
	if (_ecCur != OK) return _ecCur;

	/* Make the first pass through the data file,
	counting numbers of components */
	while (fgets(line,MAX_LINE,InFile) != NULL)    {
		/* Skip blank lines & those beginning with a comment */
		tok = strtok(line, " \t\n\r");
		if (tok == NULL) continue;
		if (*tok == ';') continue;

		/* Check if line begins with a new section heading */
		if (*tok == '[')  {
			newsect = findmatch(tok, SectTxt);
			if (newsect >= 0) {
				sect = newsect;
				if (strcmp(SectTxt[sect], "[END]") == 0) break;
				continue;
			}
			else continue;
		}

		/* Add to count of current component */
		switch(sect)  {
		case JUNCTIONS:  MaxJuncs++;    break;
		case RESERVOIRS:
		case TANKS:      MaxTanks++;    break;
		case PIPES:      MaxPipes++;    break;
		case PUMPS:      MaxPumps++;    break;
		case VALVES:     MaxValves++;   break;
		case PATTERNS:   addpattern(tok);
			break;
		case CURVES:     addcurve(tok);
			break;
		}
	}

	MaxNodes = MaxJuncs + MaxTanks;
	MaxLinks = MaxPipes + MaxPumps + MaxValves;
	if (MaxPats < 1) MaxPats = 1;
	if (_ecCur == OK)    {
		if (MaxJuncs < 1) {
			//		_ftprintf(stderr, TEXT("No nodes in the network.\n"));
			_ecCur = NO_JUNCS;   
		} else if (MaxTanks == 0) {
			//		_ftprintf(stderr, TEXT("No tank/reservoir in the network.\n"));
			_ecCur = NO_TANKS;   
		}
	}
	if (_ecCur!= OK) return _ecCur;



	// allocdata() in epanet.c, allocate memory
	/*************************************************************
	NOTE: Because network components of a given type are indexed
	starting from 1, their arrays must be sized 1
	element larger than the number of components.
	*************************************************************/

	Node = (Snode *)  calloc(MaxNodes + 1, sizeof(Snode));
	Link = (Slink *) calloc(MaxLinks + 1, sizeof(Slink));
	Adjlist = (Padjlist *) calloc(MaxNodes + 1,  sizeof(Padjlist));
	Ndx    = (int *)   calloc(MaxLinks + 1,  sizeof(int));
	Degree = (int *) calloc(MaxNodes + 1, sizeof(int));


	Order  = (int *)   calloc(MaxNodes + 1,  sizeof(int));
	Row    = (int *)   calloc(MaxNodes + 1,  sizeof(int));

	Tank    = (Stank *)    calloc(MaxTanks+1,   sizeof(Stank));
	Pump    = (Spump *)    calloc(MaxPumps+1,   sizeof(Spump));
	Valve   = (Svalve *)   calloc(MaxValves+1,  sizeof(Svalve));

	Pattern = (Spattern *) calloc(MaxPats+1,    sizeof(Spattern));
	Curve   = (Scurve *)   calloc(MaxCurves+1,  sizeof(Scurve));


	if (Node == NULL || Link == NULL || Pump == NULL || Valve == NULL
		|| Pattern == NULL || Curve == NULL) {
			//			_ftprintf(stderr, TEXT("Cannot allocate memory.\n"));
			_ecCur = MALLOC_ERROR;
			return (_ecCur);
	}

	/* Initialize pointers used in patterns, curves, and demand category lists */

	int n;
	for (n=0; n<=MaxPats; n++)
	{
		Pattern[n].Length = 0;
		Pattern[n].F = NULL;
	}
	for (n=0; n<=MaxCurves; n++)
	{
		Curve[n].Npts = 0;
		Curve[n].Type = -1;
		Curve[n].X = NULL;
		Curve[n].Y = NULL;
	}
	for (n=0; n<=MaxNodes; n++) Node[n].D = NULL;


	// second scan, load data into memory
	rewind(InFile);
	// readdata() in input2.c
	///////////////////

	char wline[MAX_LINE+1];  /* Working copy of input line      */
	int inperr; 

	/* Read each line from input file. */
	while (fgets(line,MAX_LINE,InFile) != NULL)    {
		++_lineno;

		/* Make copy of line and scan for tokens */
		strcpy(wline,line);
		Ntokens = gettokens(wline);

		/* Skip blank lines and comments */
		if (Ntokens == 0) continue;
		if (*Tok[0] == ';') continue;

		/* Check if max. length exceeded */
		if (strlen(line) >= MAX_LINE)
			return INPUT_LINE_TOO_LONG;

		/* Check if at the beginning of a new input section */
		if (*Tok[0] == '[')
		{
			newsect = findmatch(Tok[0],SectTxt);
			if (newsect >= 0)
			{
				sect = newsect;
				if (strcmp(SectTxt[sect], "[END]")  == NULL) break;
				continue;
			}
			else  { // do nothing
				//return SYNTAX_ERR_OF_SECT_NAME;
			}
		} else { // in a section's body

			/* Otherwise process next line of input in current section */
			/* newline() in input2.c */
			switch(sect)  {
			case JUNCTIONS:  _ecCur=juncdata();    break;
			case RESERVOIRS:
			case TANKS:      _ecCur=tankdata();    break;

			case EMITTERS:   _ecCur=emitterdata(); break;
			case PIPES:      _ecCur=pipedata();    break;
			case PUMPS:      _ecCur=pumpdata();   break;
			case VALVES:     _ecCur=valvedata();   break;

			case PATTERNS:   _ecCur=patterndata();   break;
			case CURVES:     _ecCur=curvedata();  break;
			case DEMANDS:	_ecCur=demanddata(); break;
			case STATUS:  _ecCur=statusdata(); break;
			case OPTIONS:	_ecCur=optiondata(); break;
				//load x-y coordinates for later building vtk graphics
			case COORDS: _ecCur=coordata(); break; 
				// ignore other sections
			}  
			if (_ecCur != OK)  {
				_lineEc = _lineno;	
				return _ecCur;
			}
		}
	}

	if ((_ecCur = unlinked()) || 
		(_ecCur = getpatterns()) || 
		(_ecCur = getcurves()) || 
		(_ecCur = getpumpparams())
		) return _ecCur;

	for (int ii = 1; ii<=MaxNodes; ++ii) 
		if (Node[ii].coordFlag == 0) {

			_ecCur = NODE_COORDINATES_UNDEFINED;
			return _ecCur;
		}

		// post-processing
		adjustdata();
		initunits(); //set up unit conversion factors

		_ecCur = inittanks();
		convertunits();

		//pre-compute link resistence
		for (int i=1; i<=Nlinks; ++i) resistence(i);

		freeTmplist(Patlist);
		freeTmplist(Curvelist);

		//file processing done.

		// compute Nusers;
		Nusers = 0;
		Nemitters = 0;
		for (int ii = 1; ii<=MaxJuncs; ++ii) {
			if (Node[ii].D && Node[ii].D->Base >0) {
				Nusers++; 
			}
			if (Node[ii].Ke  >0) Nemitters++; 
		}
        ref_demand = (float*)calloc(Nusers, sizeof(float));
        int iuser=0;
		for (int ii = 1; ii<=MaxJuncs; ++ii) {
			if (Node[ii].D && Node[ii].D->Base >0) {
                ref_demand[iuser++] = Node[ii].D->Base;
			}
		}

		// build quick indexing arrays 
		if (nPumpGPV) 
			tPumpGPV = (int*)calloc(nPumpGPV, sizeof(int));
		else
			tPumpGPV = NULL;

		if (nCRSV) 
			tCRSV = (int*)calloc(nCRSV, sizeof(int));
		else
			tCRSV = NULL;

		if (nTFV) 
			tTFV = (int*)calloc(nTFV, sizeof(int));
		else
			tTFV = NULL;

		int iPG = 0; 
		int iCRSV = 0;
		int iTFV = 0;
		for (int ii = 1; ii<=MaxLinks; ++ii) {
			switch (Link[ii].Type) {
			case CV: case PRV: case PSV:
				tCRSV[iCRSV++] = ii; break;
			case PUMP: case GPV:
				tPumpGPV[iPG++] = ii; break;
			case TCV: case FCV:
				tTFV[iTFV++] = ii; break;
			default: //do nothing
				continue;
			}
		}

		//Now ENopenH() in epanet.c and openhyd() in hydraul.c

		//createsparse() in smatrix.c,  prepare

		/* Build node-link adjacency lists with parallel links removed. */
		_ecCur = buildlists(1);  
		xparalinks();    /* Remove parallel links */
		countdegree();   /* Find degree of each junction */
		//Adjlist and Degree are prepared

		Ncoeffs = Nlinks;

		/* Re-order nodes to minimize number of non-zero coeffs.    */
		/* in factorized solution matrix. At same time, adjacency   */
		/* list is updated with links representing non-zero coeffs. */
		_ecCur = reordernodes();

		/* Allocate memory for sparse storage of positions of non-zero */
		/* coeffs. and store these positions in vector NZSUB. */
		_ecCur = storesparse(Njuncs);

		/* Free memory used for adjacency lists and sort */
		/* row indexes in NZSUB to optimize linsolve().  */
		freelists();
		_ecCur = ordersparse(Njuncs);

		/* Re-build adjacency lists without removing parallel */
		/* links for use in future connectivity checking.     */
		buildlists(FALSE);

		free(Degree);

		//we can go as far as here
		//allocmatrix() must be done in class Solver since we need parallelism
		//see Solver.cpp


		return _ecCur;
}

double Network::getBaseDemandPercentile(double pc) {
	if (pc<0 || pc>1) return -1;
	IterativePercentile<double> ip(pc);
	for (int i=1;i<MaxNodes;++i) {
		if (Node[i].D != NULL && Node[i].D->Base > 0)
			ip.add(Node[i].D->Base);
	}
	return ip.get();

}


Network::ErrorCode  Network::ordersparse(int n)
	/*
	**--------------------------------------------------------------
	** Input:   n = number of rows in solution matrix               
	** Output:  returns eror code                                   
	** Purpose: puts row indexes in ascending order in NZSUB        
	**--------------------------------------------------------------
	*/
{
	int  i, k;
	int  *xlnzt, *nzsubt, *lnzt, *nzt;
	ErrorCode  errcode = OK;

	xlnzt  = (int *) calloc(n+2, sizeof(int));
	nzsubt = (int *) calloc(Ncoeffs+2, sizeof(int));
	lnzt   = (int *) calloc(Ncoeffs+2, sizeof(int));
	nzt    = (int *) calloc(n+2, sizeof(int));

	if (xlnzt && nzsubt && lnzt && nzt)
	{

		/* Count # non-zeros in each row */
		for (i=1; i<=n; i++) nzt[i] = 0;
		for (i=1; i<=n; i++)
		{
			for (k=XLNZ[i]; k<XLNZ[i+1]; k++) nzt[NZSUB[k]]++;
		}
		xlnzt[1] = 1;
		for (i=1; i<=n; i++) xlnzt[i+1] = xlnzt[i] + nzt[i];

		/* Transpose matrix twice to order column indexes */
		transpose(n,XLNZ,NZSUB,LNZ,xlnzt,nzsubt,lnzt,nzt);
		transpose(n,xlnzt,nzsubt,lnzt,XLNZ,NZSUB,LNZ,nzt);
	} else errcode = MALLOC_ERROR;

	/* Reclaim memory */
	free(xlnzt);
	free(nzsubt);
	free(lnzt);
	free(nzt);
	return(errcode);
}                        /* End of ordersparse */

void  Network::transpose(int n, int *il, int *jl, int *xl, int *ilt, int *jlt,
						 int *xlt, int *nzt)
						 /*
						 **---------------------------------------------------------------------
						 ** Input:   n = matrix order                                           
						 **          il,jl,xl = sparse storage scheme for original matrix       
						 **          nzt = work array                                           
						 ** Output:  ilt,jlt,xlt = sparse storage scheme for transposed matrix  
						 ** Purpose: Determines sparse storage scheme for transpose of a matrix 
						 **---------------------------------------------------------------------
						 */
{
	int  i, j, k, kk;

	for (i=1; i<=n; i++) nzt[i] = 0;
	for (i=1; i<=n; i++)
	{
		for (k=il[i]; k<il[i+1]; k++)
		{
			j = jl[k];
			kk = ilt[j] + nzt[j];
			jlt[kk] = i;
			xlt[kk] = xl[k];
			nzt[j]++;
		}
	}
}                        /* End of transpose */



void  Network::freelists()
	/*
	**--------------------------------------------------------------
	** Input:   none                                                
	** Output:  none                                                
	** Purpose: frees memory used for nodal adjacency lists         
	**--------------------------------------------------------------
	*/
{
	int   i;
	Padjlist alink;

	for (i=0; i<=Nnodes; i++)
	{
		for (alink = Adjlist[i]; alink != NULL; alink = Adjlist[i])
		{
			Adjlist[i] = alink->next;
			free(alink);
		}
	}
}                        /* End of freelists */



Network::ErrorCode  Network::storesparse(int n)
	/*
	**--------------------------------------------------------------
	** Input:   n = number of rows in solution matrix               
	** Output:  returns error code                                  
	** Purpose: stores row indexes of non-zeros of each column of   
	**          lower triangular portion of factorized matrix       
	**--------------------------------------------------------------
	*/
{
	Padjlist alink;
	int   i, ii, j, k, l, m;
	ErrorCode   errcode = OK;

	/* Allocate sparse matrix storage */
	XLNZ  = (int *) calloc(n+2, sizeof(int));
	NZSUB = (int *) calloc(Ncoeffs+2, sizeof(int));
	LNZ   = (int *) calloc(Ncoeffs+2, sizeof(int));
	if (XLNZ == NULL || NZSUB == NULL || LNZ == NULL) return(MALLOC_ERROR);

	/* Generate row index pointers for each column of matrix */
	k = 0;
	XLNZ[1] = 1;
	for (i=1; i<=n; i++)             /* column */
	{
		m = 0;
		ii = Order[i];
		for (alink = Adjlist[ii]; alink != NULL; alink = alink->next)
		{
			j = Row[alink->node];    /* row */
			l = alink->link;
			if (j > i && j <= n)
			{
				m++;
				k++;
				NZSUB[k] = j;
				LNZ[k] = l;
			}
		}
		XLNZ[i+1] = XLNZ[i] + m;
	}
	return(errcode);
}                        /* End of storesparse */


void  Network::xparalinks()
	/*
	**--------------------------------------------------------------
	** Input:   none                                                
	** Output:  none                                                
	** Purpose: removes parallel links from nodal adjacency lists   
	**--------------------------------------------------------------
	*/
{
	int    i;
	Padjlist  alink,       /* Current item in adjacency list */
		blink;       /* Previous item in adjacency list */

	/* Scan adjacency list of each node */
	for (i=1; i<=Nnodes; i++)
	{
		alink = Adjlist[i];              /* First item in list */
		blink = NULL;
		while (alink != NULL)
		{
			if (alink->node == 0)      /* Parallel link marker found */
			{
				if (blink == NULL)      /* This holds at start of list */
				{
					Adjlist[i] = alink->next;
					free(alink);             /* Remove item from list */
					alink = Adjlist[i];
				}
				else                    /* This holds for interior of list */
				{
					blink->next = alink->next;
					free(alink);             /* Remove item from list */
					alink = blink->next;
				}
			}
			else
			{
				blink = alink;          /* Move to next item in list */
				alink = alink->next;
			}
		}
	}
}                        /* End of xparalinks */

void  Network::countdegree()
	/*
	**----------------------------------------------------------------
	** Input:   none                                                
	** Output:  none                                                
	** Purpose: counts number of nodes directly connected to each node            
	**----------------------------------------------------------------
	*/
{
	int   i;
	Padjlist alink;
	memset(Degree,0,(Nnodes+1)*sizeof(int));

	/* NOTE: For purposes of node re-ordering, Tanks (nodes with  */
	/*       indexes above Njuncs) have zero degree of adjacency. */

	for (i=1; i<=Njuncs; i++)
		for (alink = Adjlist[i]; alink != NULL; alink = alink->next)
			if (alink->node > 0) Degree[i]++;
}

Network::ErrorCode  Network::reordernodes()
	/*
	**--------------------------------------------------------------
	** Input:   none                                                
	** Output:  returns 1 if successful, 0 if not                   
	** Purpose: re-orders nodes to minimize # of non-zeros that     
	**          will appear in factorized solution matrix           
	**--------------------------------------------------------------
	*/
{
	int k, knode, m, n;
	for (k=1; k<=Nnodes; k++)
	{
		Row[k] = k;
		Order[k] = k;
	}
	n = Njuncs;
	for (k=1; k<=n; k++)                   /* Examine each junction    */
	{
		m = mindegree(k,n);                 /* Node with lowest degree  */
		knode = Order[m];                   /* Node's index             */
		if (!growlist(knode)) return(MALLOC_ERROR);  /* Augment adjacency list   */
		Order[m] = Order[k];                /* Switch order of nodes    */
		Order[k] = knode;
		Degree[knode] = 0;                  /* In-activate node         */
	}
	for (k=1; k<=n; k++)                   /* Assign nodes to rows of  */
		Row[Order[k]] = k;                   /*   coeff. matrix          */
	return(OK);
} 


int  Network::mindegree(int k, int n)
	/*
	**--------------------------------------------------------------
	** Input:   k = first node in list of active nodes              
	**          n = total number of junction nodes                  
	** Output:  returns node index with fewest direct connections     
	** Purpose: finds active node with fewest direct connections
	**--------------------------------------------------------------
	*/
{
	int i, m;
	int min = n,
		imin = n;

	for (i=k; i<=n; i++)
	{
		m = Degree[Order[i]];
		if (m < min)
		{
			min = m;
			imin = i;
		}
	}
	return(imin);
}                        /* End of mindegree */


int  Network::growlist(int knode)
	/*
	**--------------------------------------------------------------
	** Input:   knode = node index                                  
	** Output:  returns 1 if successful, 0 if not                   
	** Purpose: creates new entries in knode's adjacency list for   
	**          all unlinked pairs of active nodes that are         
	**          adjacent to knode                                   
	**--------------------------------------------------------------
	*/
{
	int   node;
	Padjlist alink;

	/* Iterate through all nodes connected to knode */
	for (alink = Adjlist[knode]; alink != NULL; alink = alink -> next)
	{
		node = alink->node;       /* End node of connecting link  */
		if (Degree[node] > 0)     /* End node is active           */
		{
			Degree[node]--;        /* Reduce degree of adjacency   */
			if (!newlink(alink))   /* Add to adjacency list        */
				return(0);
		}
	}
	return(1);
}                        /* End of growlist */


int  Network::newlink(Padjlist alink)
	/*
	**--------------------------------------------------------------
	** Input:   alink = element of node's adjacency list            
	** Output:  returns 1 if successful, 0 if not                   
	** Purpose: links end of current adjacent link to end nodes of  
	**          all links that follow it on adjacency list          
	**--------------------------------------------------------------
	*/
{
	int   inode, jnode;
	Padjlist blink;

	/* Scan all entries in adjacency list that follow anode. */
	inode = alink->node;             /* End node of connection to anode */
	for (blink = alink->next; blink != NULL; blink = blink->next)
	{
		jnode = blink->node;          /* End node of next connection */

		/* If jnode still active, and inode not connected to jnode, */
		/* then add a new connection between inode and jnode.       */
		if (Degree[jnode] > 0)        /* jnode still active */
		{
			if (!linked(inode,jnode))  /* inode not linked to jnode */
			{

				/* Since new connection represents a non-zero coeff. */
				/* in the solution matrix, update the coeff. count.  */
				Ncoeffs++;

				/* Update adjacency lists for inode & jnode to */
				/* reflect the new connection.                 */
				if (!addlink(inode,jnode,Ncoeffs)) return(0);
				if (!addlink(jnode,inode,Ncoeffs)) return(0);
				Degree[inode]++;
				Degree[jnode]++;
			}
		}
	}
	return(1);
}                        /* End of newlink */


int  Network::linked(int i, int j)
	/*
	**--------------------------------------------------------------
	** Input:   i = node index                                      
	**          j = node index                                      
	** Output:  returns 1 if nodes i and j are linked, 0 if not     
	** Purpose: checks if nodes i and j are already linked.         
	**--------------------------------------------------------------
	*/
{
	Padjlist alink;
	for (alink = Adjlist[i]; alink != NULL; alink = alink->next)
		if (alink->node == j) return(1);
	return(0);
}                        /* End of linked */


int  Network::addlink(int i, int j, int n)
	/*
	**--------------------------------------------------------------
	** Input:   i = node index                                      
	**          j = node index                                      
	**          n = link index                                      
	** Output:  returns 1 if successful, 0 if not                   
	** Purpose: augments node i's adjacency list with node j        
	**--------------------------------------------------------------
	*/
{
	Padjlist alink;
	alink = (struct Sadjlist *) malloc(sizeof(struct Sadjlist));
	if (alink == NULL) return(0);
	alink->node = j;
	alink->link = n;
	alink->next = Adjlist[i];
	Adjlist[i] = alink;
	return(1);
}                        /* End of addlink */

void Network::resistence(int k /* link index */) {
	double e,d,L;
	Link[k].R = _TINY;
	switch (Link[k].Type)
	{

		/* Link is a pipe. Compute resistance based on headloss formula. */
		/* Friction factor for D-W formula gets included during solution */
		/* process in pipecoeff() function.                              */
	case CV:
	case PIPE: 
		e = Link[k].Kc;                 /* Roughness coeff. */
		d = Link[k].Diam;               /* Diameter */
		L = Link[k].Len;                /* Length */
		switch(Formflag)
		{
		case HW: Link[k].R = 4.727*L/pow(e,Hexp)/pow(d,4.871);
			break;
		case DW: Link[k].R = L/2.0/32.2/d/((_PI*d*d/4.0)*(_PI*d*d/4.0));
			break;
		case CM: Link[k].R = (4.0*e/(1.49*_PI*d*d))*(4.0*e/(1.49*_PI*d*d))
					 *pow((d/4.0),-1.333)*L;
		}
		break;

		/* Link is a pump. Use negligible resistance. */
	case PUMP:
		//??   Link[k].R = CBIG;  //CSMALL;
		break;
	}
}

void Network::initunits()
	/*
	**--------------------------------------------------------------
	**  Input:   none
	**  Output:  none
	**  Purpose: determines unit conversion factors
	**			Ucfs convert everything to 
	cfs for flow rates
	psi for (free) pressure
	feet for elevation, total head, headloss, tank levels, pipe length
	inch for pipe diameter
	feet per second for velocity
	**--------------------------------------------------------------
	*/
{
	double  dcf,  /* distance conversion factor      */
		qcf,  /* flow conversion factor          */
		hcf,  /* head conversion factor          */
		pcf,  /* pressure conversion factor      */
		wcf;  /* energy conversion factor        */

	if (Unitsflag == SI)                            /* SI units */
	{
		dcf = 304.8; //  km/feet
		qcf = 28.317;  // lps/cfs
		if (Flowflag == LPM) qcf = 1699;
		if (Flowflag == MLD) qcf = 2.4466;
		if (Flowflag == CMH) qcf = 101.94;
		if (Flowflag == CMD) qcf = 2446.6;
		hcf = 0.3048;
		if (Pressflag == METERS) pcf = 0.3048;//  m/feet
		else pcf = 6.895*0.4333;
		wcf = 0.7457;
	}
	else                                         /* US units */
	{
		dcf = 12.0;
		qcf = 1.0;
		if (Flowflag == GPM) qcf = 448.831;
		if (Flowflag == MGD) qcf = 0.64632;
		if (Flowflag == IMGD)qcf = 0.5382;
		if (Flowflag == AFD) qcf = 1.9837;
		hcf = 1.0;
		pcf = 0.4333;
		wcf = 1.0;
	}

	Ucf[DEMAND]    = qcf;
	Ucf[ELEV]      = hcf;
	Ucf[HEAD]      = hcf;
	Ucf[PRESSURE]  = pcf;
	Ucf[LENGTH]    = hcf;
	Ucf[DIAM]      = dcf;
	Ucf[FLOW]      = qcf;
	Ucf[VELOCITY]  = hcf;
	Ucf[HEADLOSS]  = hcf;
	Ucf[FRICTION]  = 1.0;
	Ucf[POWER]     = wcf;
	Ucf[VOLUME]    = hcf*hcf*hcf;
}

void  Network::convertunits()
	/*
	**--------------------------------------------------------------
	**  Input:   none
	**  Output:  none
	**  Purpose: converts units of input data
	**--------------------------------------------------------------
	*/
{
	int   i,j,k;
	double ucf;        /* Unit conversion factor */
	Pdemand demand;   /* Pointer to demand record */

	/* Convert nodal elevations & initial WQ */
	/* (WQ source units are converted in QUALITY.C */
	for (i=1; i<=Nnodes; i++)
	{
		Node[i].El /= Ucf[ELEV];
	}

	/* Convert demands */
	for (i=1; i<=Njuncs; i++)
	{
		for (demand = Node[i].D; demand != NULL; demand = demand->next)
			demand->Base /= Ucf[DEMAND];
	}
	/* Convert emitter discharge coeffs. to head loss coeff. */
	ucf = pow(Ucf[FLOW],Qexp)/Ucf[PRESSURE];
	for (i=1; i<=Njuncs; i++)
		if (Node[i].Ke > 0.0) Node[i].Ke = ucf/pow(Node[i].Ke,Qexp);


	/* Initialize tank variables (convert tank levels to elevations) */
	for (j=1; j<=Ntanks; j++)
	{
		i = Tank[j].Node;
		Tank[j].H0 = Node[i].El + Tank[j].H0/Ucf[ELEV];
		Tank[j].Hmin = Node[i].El + Tank[j].Hmin/Ucf[ELEV];
		Tank[j].Hmax = Node[i].El + Tank[j].Hmax/Ucf[ELEV];
		Tank[j].A = _PI*(Tank[j].A/Ucf[ELEV])*(Tank[j].A/Ucf[ELEV])/4.0;
		Tank[j].V0 /= Ucf[VOLUME];
		Tank[j].Vmin /= Ucf[VOLUME];
		Tank[j].Vmax /= Ucf[VOLUME];
		Tank[j].Kb /= 86400;
		Tank[j].V = Tank[j].V0;
	}

	/* Convert units of link parameters */
	for (k=1; k<=Nlinks; k++)
	{
		if (Link[k].Type <= PIPE)
		{
			/* Convert pipe parameter units:                         */
			/*    - for Darcy-Weisbach formula, convert roughness    */
			/*      from millifeet (or mm) to ft (or m)              */
			/*    - for US units, convert diameter from inches to ft */
			if (Formflag  == DW) Link[k].Kc /= (1000.0*Ucf[ELEV]);
			Link[k].Diam /= Ucf[DIAM];
			Link[k].Len /= Ucf[LENGTH];

			/* Convert minor loss coeff. from V^2/2g basis to Q^2 basis */
			Link[k].Km = 0.02517*Link[k].Km / 
				((Link[k].Diam) * (Link[k].Diam))
				/ ((Link[k].Diam) * (Link[k].Diam));

		} else if (Link[k].Type == PUMP )
		{
			/* Convert units for pump curve parameters */
			i = (int)Link[k].Diam;
			if (Pump[i].Ptype == CONST_HP)
			{
				/* For constant hp pump, convert kw to hp */
				if (Unitsflag == SI) Pump[i].R /= Ucf[POWER];
			}
			else
			{
				/* For power curve pumps, convert     */
				/* shutoff head and flow coefficient  */
				if (Pump[i].Ptype == POWER_FUNC)
				{
					Pump[i].H0 /= Ucf[HEAD];
					Pump[i].R  *= (pow(Ucf[FLOW],Pump[i].N)/Ucf[HEAD]);
				}
				/* Convert flow range & max. head units */
				Pump[i].Q0   /= Ucf[FLOW];
				Pump[i].Qmax /= Ucf[FLOW];
				Pump[i].Hmax /= Ucf[HEAD];
			}
		}

		else
		{      /* For flow control valves, convert flow setting    */
			/* while for other valves convert pressure setting  */
			Link[k].Diam /= Ucf[DIAM];
			Link[k].Km = 0.02517*Link[k].Km/pow(Link[k].Diam, 4);      
			if (Link[k].Kc != 0) switch (Link[k].Type)
			{
				//            case FCV: Link[k].Kc /= Ucf[FLOW]; break;
			case PRV:
			case PSV:
			case PBV: Link[k].Kc /= Ucf[PRESSURE]; break;
			}
		}

	}

}                       /*  End of convertunits  */

Network::ErrorCode Network::inittanks() 
	/*
	**---------------------------------------------------------------
	**  Input:   none
	**  Output:  returns error code
	**  Purpose: initializes volumes in non-cylindrical tanks
	**---------------------------------------------------------------
	*/
{
	int   i,j,n = 0;
	double a;
	int   levelerr;
	ErrorCode errcode = OK;

	for (j=1; j<=Ntanks; j++)
	{

		/* Skip reservoirs */
		if (Tank[j].A == 0.0) continue;

		/* Check for valid lower/upper tank levels */
		levelerr = 0;
		if (Tank[j].H0   > Tank[j].Hmax ||
			Tank[j].Hmin > Tank[j].Hmax ||
			Tank[j].H0   < Tank[j].Hmin
			) levelerr = 1;

		/* Check that tank heights are within volume curve */
		i = Tank[j].Vcurve;
		if (i > 0)
		{
			n = Curve[i].Npts - 1;
			if (Tank[j].Hmin < Curve[i].X[0] ||
				Tank[j].Hmax > Curve[i].X[n]
			) levelerr = 1;
		}

		/* Report error in levels if found */
		if (levelerr)
		{
			return(TANK_LEVEL_INVALID);
		}

		/* Else if tank has a volume curve, */
		else if (i > 0)
		{
			/* Find min., max., and initial volumes from curve */
			Tank[j].Vmin = interp(Curve[i].Npts,Curve[i].X,
				Curve[i].Y,Tank[j].Hmin);
			Tank[j].Vmax = interp(Curve[i].Npts,Curve[i].X,
				Curve[i].Y,Tank[j].Hmax);
			Tank[j].V0   = interp(Curve[i].Npts,Curve[i].X,
				Curve[i].Y,Tank[j].H0);

			/* Find a "nominal" diameter for tank */
			a = (Curve[i].Y[n] - Curve[i].Y[0])/
				(Curve[i].X[n] - Curve[i].X[0]);
			Tank[j].A = sqrt(4.0*a/_PI);
		}
	}
	return(errcode);
}                       /* End of inittanks */

double  Network::interp(int n, double x[], double y[], double xx)
	/*----------------------------------------------------------------
	**  Input:   n  = number of data pairs defining a curve
	**           x  = x-data values of curve
	**           y  = y-data values of curve
	**           xx = specified x-value
	**  Output:  none
	**  Returns: y-value on curve at x = xx
	**  Purpose: uses linear interpolation to find y-value on a
	**           data curve corresponding to specified x-value.
	**  NOTE:    does not extrapolate beyond endpoints of curve.
	**----------------------------------------------------------------
	*/
{
	int    k,m;
	double  dx,dy;

	m = n - 1;                          /* Highest data index      */
	if (xx <= x[0]) return(y[0]);       /* xx off low end of curve */
	for (k=1; k<=m; k++)                /* Bracket xx on curve     */
	{
		if (x[k] >= xx)                 /* Interp. over interval   */
		{
			dx = x[k]-x[k-1];
			dy = y[k]-y[k-1];
			if (abs(dx) < _TINY) return(y[k]);
			else return(y[k] - (x[k]-xx)*dy/dx);
		}
	}
	return(y[m]);                       /* xx off high end of curve */
}                       /* End of interp */

void Network::adjustdata()
{
	int   i;
	double ucf;                   /* Unit conversion factor */
	Pdemand demand;              /* Pointer to demand record */

	/* Determine unit system based on flow units */
	switch (Flowflag)
	{
	case LPS:          /* Liters/sec */
	case LPM:          /* Liters/min */
	case MLD:          /* megaliters/day  */
	case CMH:          /* cubic meters/hr */
	case CMD:          /* cubic meters/day */
		Unitsflag = SI;
		break;
	default:
		Unitsflag = US;
	}
	if (Formflag == HW) Hexp = 1.852;
	else                Hexp = 2.0;

	/* Revise pressure units depending on flow units */
	if (Unitsflag != SI) Pressflag = PSI;
	else if (Pressflag == PSI) Pressflag = METERS;

	/* Use default pattern if none assigned to a demand */
	for (i=1; i<=Nnodes; i++)
	{
		for (demand = Node[i].D; demand != NULL; demand = demand->next)
			if (demand->Pat == 0) demand->Pat = DefPat;
	}

}                       /*  End of adjustdata  */


Network::ErrorCode Network::getpumpparams()
	/*
	**-------------------------------------------------------------
	**  Input:   none
	**  Output:  returns error code
	**  Purpose: computes & checks pump curve parameters
	**--------------------------------------------------------------
	*/
{
	int   i, j = 0, k, m, n = 0;
	double a,b,c,
		h0 = 0.0, h1 = 0.0, h2 = 0.0, q1 = 0.0, q2 = 0.0;

	for (i=1; i<=Npumps; i++)
	{
		k = Pump[i].Link;
		if (Pump[i].Ptype == CONST_HP)      /* Constant Hp pump */
		{
			Pump[i].H0 = 0.0;
			Pump[i].R  = -8.814*Link[k].Km;
			Pump[i].N  = -1.0;
			Pump[i].Hmax  = _BIG;             /* No head limit      */
			Pump[i].Qmax  = _BIG;             /* No flow limit      */
			Pump[i].Q0 = 1.0;                /* Init. flow = 1 cfs */
			continue;
		}

		/* Set parameters for pump curves */
		else if (Pump[i].Ptype == NOCURVE)  /* Pump curve specified */
		{
			j = Pump[i].Hcurve;              /* Get index of head curve */
			if (j == 0)
			{                                /* Error: No head curve */
				return(PUMP_CURVE_UNDEFINED);
			}
			n = Curve[j].Npts;
			if (n == 1)                      /* Only a single h-q point */
			{                                /* supplied so use generic */
				Pump[i].Ptype = POWER_FUNC;   /* power function curve.   */
				q1 = Curve[j].X[0];
				h1 = Curve[j].Y[0];
				h0 = 1.33334*h1;
				q2 = 2.0*q1;
				h2 = 0.0;
			}
			else if (n == 3
				&&  Curve[j].X[0] == 0.0)   /* 3 h-q points supplied with */
			{                                /* shutoff head so use fitted */   
				Pump[i].Ptype = POWER_FUNC;   /* power function curve.      */
				h0 = Curve[j].Y[0];
				q1 = Curve[j].X[1];
				h1 = Curve[j].Y[1];
				q2 = Curve[j].X[2];
				h2 = Curve[j].Y[2];
			}
			else Pump[i].Ptype = CUSTOM;     /* Else use custom pump curve.*/

			/* Compute shape factors & limits of power function pump curves */
			if (Pump[i].Ptype == POWER_FUNC)
			{
				if (!powercurve(h0,h1,h2,q1,q2,&a,&b,&c))
				{                             /* Error: Invalid curve */ 
					return(PUMP_CURVE_INVALID);
				}
				else
				{
					Pump[i].H0 = -a;
					Pump[i].R  = -b;
					Pump[i].N  = c;
					Pump[i].Q0 = q1;
					Pump[i].Qmax  = pow((-a/b),(1.0/c));
					Pump[i].Hmax  = h0;
				}
			}
		}

		/* Assign limits to custom pump curves */
		if (Pump[i].Ptype == CUSTOM)
		{
			for (m=1; m<n; m++)
			{
				if (Curve[j].Y[m] >= Curve[j].Y[m-1])
				{                             /* Error: Invalid curve */
					return(PUMP_CURVE_INVALID);
				}
			}
			Pump[i].Qmax  = Curve[j].X[n-1];
			Pump[i].Q0    = (Curve[j].X[0] + Pump[i].Qmax)/2.0;
			Pump[i].Hmax  = Curve[j].Y[0];
		}
	}   /* Next pump */
	return(OK);
}

Network::ErrorCode Network::getcurves()
	/*
	**-----------------------------------------------------------
	**  Input:   none
	**  Output:  returns error code
	**  Purpose: retrieves curve data from temporary linked list
	**-----------------------------------------------------------
	*/
{
	int i,j;
	double x;
	SFloatlist *fx, *fy;
	STmplist *c;

	/* Start at head of curve list */
	c = Curvelist;

	/* Traverse list of curves */
	while (c != NULL)
	{
		i = c->i;
		if (i >= 1 && i <= MaxCurves)
		{

			/* Save curve ID */
			strcpy(Curve[i].ID, c->ID);

			/* Check that curve has data points */
			if (Curve[i].Npts <= 0)
			{
				return(CURVE_HAS_NO_DATA);
			}

			/* Allocate memory for curve data */
			Curve[i].X = (double *) calloc(Curve[i].Npts, sizeof(double));
			Curve[i].Y = (double *) calloc(Curve[i].Npts, sizeof(double));
			if (Curve[i].X == NULL || Curve[i].Y == NULL) return(MALLOC_ERROR);

			/* Traverse list of x,y data */
			x = _BIG;
			fx = c->x;
			fy = c->y;
			j = Curve[i].Npts - 1;
			while (fx != NULL && fy != NULL && j >= 0)
			{

				/* Check that x data is in ascending order */
				if (fx->value >= x)
				{
					return(CURVE_DATA_NOT_ASC);
				}
				x = fx->value;

				/* Save x,y data in Curve structure */
				Curve[i].X[j] = fx->value;
				fx = fx->next;
				Curve[i].Y[j] = fy->value;
				fy = fy->next;
				j--;
			}
		}
		c = c->next;
	}
	return(OK);
}


Network::ErrorCode  Network::unlinked()
	/*
	**--------------------------------------------------------------
	** Input:   none                                                
	** Output:  returns error code if any unlinked junctions found  
	** Purpose: checks for unlinked junctions in network            
	**                                                              
	** NOTE: unlinked tanks have no effect on computations.         
	**--------------------------------------------------------------
	*/
{
	char  *marked;
	int   i;
	ErrorCode errcode = OK;

	marked   = (char *) calloc(Nnodes+1,sizeof(char));
	if (marked == NULL) return (MALLOC_ERROR);
	memset(marked,0,(Nnodes+1)*sizeof(char));
	for (i=1; i<=Nlinks; i++)            /* Mark end nodes of each link */
	{
		marked[Link[i].N1]++;
		marked[Link[i].N2]++;
	}
	for (i=1; i<=Njuncs; i++)            /* Check each junction  */
	{
		if (marked[i] == 0)               /* If not marked then error */
		{
			errcode = NODE_NOT_LINKED;
		}
	}

	free(marked);
	return(errcode);
}                        /* End of unlinked */


int  Network::paralink(int i, int j, int k)
	/*
	**--------------------------------------------------------------
	** Input:   i = index of start node of link                    
	**          j = index of end node of link                    
	**          k = link index                                    
	** Output:  returns 1 if link k parallels another link, else 0
	** Purpose: checks for parallel links between nodes i and j   
	**                                                            
	**--------------------------------------------------------------
	*/
{
	Padjlist alink;
	for (alink = Adjlist[i]; alink != NULL; alink = alink->next)
	{
		if (alink->node == j)     /* Link || to k (same end nodes) */
		{
			Ndx[k] = alink->link;  /* Assign Ndx entry to this link */
			return(1);
		}
	}
	Ndx[k] = k;                  /* Ndx entry if link not parallel */
	return(0);
}                        /* End of paralink */


Network::ErrorCode  Network::buildlists(int paraflag)
	/*
	**--------------------------------------------------------------
	** Input:   paraflag = TRUE if list marks parallel links      
	** Output:  returns error code                                
	** Purpose: builds linked list of links adjacent to each node 
	**--------------------------------------------------------------
	*/
{
	int    i,j,k;
	int    pmark = 0;
	ErrorCode    errcode = OK;
	Padjlist  alink;

	/* For each link, update adjacency lists of its end nodes */
	for (k=1; k<=MaxLinks; k++)
	{
		i = Link[k].N1;
		j = Link[k].N2;
		if (paraflag) pmark = paralink(i,j,k);  /* Parallel link check */

		/* Include link in start node i's list */
		alink = (struct Sadjlist *) malloc(sizeof(struct Sadjlist));
		if (alink == NULL) return(MALLOC_ERROR);
		if (!pmark) alink->node = j;
		else        alink->node = 0;           /* Parallel link marker */
		alink->link = k;
		alink->next = Adjlist[i];
		Adjlist[i] = alink;

		/* Include link in end node j's list */
		alink = (struct Sadjlist *) malloc(sizeof(struct Sadjlist));
		if (alink == NULL) return(MALLOC_ERROR);
		if (!pmark) alink->node = i;
		else        alink->node = 0;           /* Parallel link marker */
		alink->link = k;
		alink->next = Adjlist[j];
		Adjlist[j] = alink;
	}
	return(errcode);
}                        /* End of buildlists */


Network::ErrorCode Network::getpatterns()
	/*
	**-----------------------------------------------------------
	**  Input:   none
	**  Output:  returns error code
	**  Purpose: retrieves pattern data from temporary linked list
	**-------------------------------------------------------------
	*/
{
	int i,j;
	SFloatlist *f;
	STmplist *pat;

	/* Start at head of list */
	pat = Patlist;

	/* Traverse list of patterns */
	while (pat != NULL)
	{

		/* Get index of current pattern in Pattern array */
		i = pat->i;

		/* Check if this is the default pattern */
		if (strcmp(pat->ID, DefPatID) == 0) DefPat = i;
		if (i >= 0 && i <= MaxPats)
		{
			/* Save pattern ID */
			strcpy(Pattern[i].ID, pat->ID);

			/* Give pattern a length of at least 1 */
			if (Pattern[i].Length == 0) Pattern[i].Length = 1;
			Pattern[i].F = (double *) calloc(Pattern[i].Length, sizeof(double));
			if (Pattern[i].F == NULL) return(MALLOC_ERROR);

			/* Start at head of pattern multiplier list */
			/* (which holds multipliers in reverse order)*/
			f = pat->x;
			j = Pattern[i].Length - 1;

			/* Use at least one multiplier equal to 1.0 */
			if (f == NULL) Pattern[i].F[0] = 1.0;

			/* Traverse list, storing multipliers in Pattern array */
			else while (f != NULL && j >= 0)
			{
				Pattern[i].F[j] = f->value;
				f = f->next;
				j--;
			}
		}
		pat = pat->next;
	}
	return(OK);
}


Network::ErrorCode  Network::optiondata()
	/*
	**--------------------------------------------------------------
	**  Input:   none                                                
	**  Output:  returns error code                                  
	**  Purpose: processes [OPTIONS] data                            
	**--------------------------------------------------------------
	*/
{
	int n;
	n = Ntokens - 1;
	return(optionchoice(n));         /* Option is a named choice    */
	/* all numerical options are ignored */
}                        /* end of optiondata */


Network::ErrorCode  Network::optionchoice(int n)
	/*
	**--------------------------------------------------------------
	**  Input:   n = index of last input token saved in Tok[]          
	**  Output:  returns error code or 0 if option belongs to        
	**           those listed below, or  otherwise                 
	**  Purpose: processes fixed choice [OPTIONS] data               
	**  Formats:                                                     
	**    UNITS               CFS/GPM/MGD/IMGD/AFD/LPS/LPM/MLD/CMH/CMD/SI
	**    PRESSURE            PSI/KPA/M                          
	**    HEADLOSS            H-W/D-W/C-M                        
	**    PATTERN             id
	**--------------------------------------------------------------
	*/
{
	/* Check if 1st token matches a parameter name and */
	/* process the input for the matched parameter     */
	if (n < 0) return(SYNTAX_ERR_OF_OPTIONS);
	if (match(Tok[0],"UNITS"))
	{
		if (n < 1) return(OK);
		else if (match(Tok[1],"CFS"))  Flowflag = CFS;
		else if (match(Tok[1],"GPM"))  Flowflag = GPM;
		else if (match(Tok[1],"AFD"))  Flowflag = AFD;
		else if (match(Tok[1],"MGD"))  Flowflag = MGD;
		else if (match(Tok[1],"IMGD")) Flowflag = IMGD;
		else if (match(Tok[1],"LPS"))  Flowflag = LPS;
		else if (match(Tok[1],"LPM"))  Flowflag = LPM;
		else if (match(Tok[1],"CMH"))  Flowflag = CMH;
		else if (match(Tok[1],"CMD"))  Flowflag = CMD;
		else if (match(Tok[1],"MLD"))  Flowflag = MLD;
		else if (match(Tok[1],"SI"))   Flowflag = LPS;
		else return(SYNTAX_ERR_OF_OPTIONS);
	}
	else if (match(Tok[0],"PRESSURE"))
	{
		if (n < 1) return(OK);
		else if (match(Tok[1],"PSI"))    Pressflag = PSI;
		else if (match(Tok[1],"KPA"))    Pressflag = KPA;
		else if (match(Tok[1],"METERS")) Pressflag = METERS;
		else return(SYNTAX_ERR_OF_OPTIONS);
	}
	else if (match(Tok[0],"HEADLOSS"))
	{
		if (n < 1) return(OK);
		else if (match(Tok[1],"HW") || match(Tok[1],"H-W")) Formflag = HW;
		else if (match(Tok[1],"DW") || match(Tok[1],"D-W")) Formflag = DW;
		else if (match(Tok[1],"CM")|| match(Tok[1],"C-M")) Formflag = CM;
		else return(SYNTAX_ERR_OF_OPTIONS);
	}

	else if (match(Tok[0],"PATTERN"))            /* Pattern option */
	{
		if (n < 1) return(OK);
		strncpy(DefPatID,Tok[1],MAX_ID); //update default pattern
		// possible mismatch here
	}

	return (OK); //ignore unknown options

}                        /* end of optionchoice */



Network::ErrorCode  Network::juncdata()
	/*
	**--------------------------------------------------------------
	**  Input:   none                                                
	**  Output:  returns error code                                  
	**  Purpose: processes junction data                             
	**  Format:                                                      
	**    [JUNCTIONS]                                              
	**      id  elev.  (demand)  (demand pattern)                  
	**--------------------------------------------------------------
	*/
{
	int      n, p = 0;
	double    el,y = 0.0;
	Pdemand  demand;
	STmplist *pat;


	/* Add new junction to data base */
	n = Ntokens;
	if (Nnodes == MaxNodes) return(SYNTAX_ERR_OF_JUNC);
	Njuncs++;
	Nnodes++;

	char*	junc_name = new char[MAX_ID];
	strncpy_s(junc_name, MAX_ID, Tok[0], MAX_ID-1);
	strncpy(Node[Njuncs].ID, junc_name, MAX_ID-1);
	if (!(Nht.insert(HTPair(junc_name, Njuncs))).second)  //id exists
		return(NODE_ID_EXISTS);

	/* Check for valid data */
	if (n < 2) return(SYNTAX_ERR_OF_JUNC);
	if (!getfloat(Tok[1],&el)) return(ILLEGAL_NUMBER);
	if (n >= 3  && !getfloat(Tok[2],&y)) return(ILLEGAL_NUMBER);
	if (n >= 4)
	{
		pat = findID(Tok[3],Patlist);
		if (pat == NULL) return(DEMAND_PATTERN_UNDEFINED);
		p = pat->i;
	}

	/* Save junction data */
	Node[Njuncs].El  = el;
	Node[Njuncs].Ke  = 0.0;

	/* Create a new demand record */
	/*** Updated 6/24/02 ***/
	if (n >= 3)
	{
		demand = (struct Sdemand *) malloc(sizeof(struct Sdemand));
		if (demand == NULL) return(MALLOC_ERROR);
		demand->Base = y;
		demand->Pat = p;
		demand->next = Node[Njuncs].D;
		Node[Njuncs].D = demand;
	}
	/*** end of update ***/
	return(OK);
}                        /* end of juncdata */

Network::ErrorCode  Network::emitterdata()
	/*
	**--------------------------------------------------------------
	**  Input:   none                                                
	**  Output:  returns error code                                  
	**  Purpose: processes junction emitter data                     
	**  Format:                                                     
	**     [EMITTER]                                               
	**        node   K                                     
	**--------------------------------------------------------------
	*/
{
	int   j,                  /* Node index    */
		n;                  /* # data items  */
	double k;                  /* Flow coeff,   */

	n = Ntokens;
	if (n < 2) return(SYNTAX_ERR_OF_EMITTERS); 

	HTIt inode = Nht.find(Tok[0]);
	if (inode == Nht.end()) 
		return (NODE_UNDEFINED);

	if (inode->second > Njuncs) return(NODE_UNDEFINED);   /* Not a junction.*/
	if (!getfloat(Tok[1],&k)) return(ILLEGAL_NUMBER);
	if (k < 0.0) return(ILLEGAL_NUMBER);
	Node[inode->second].Ke = k;
	return OK;
}


Network::ErrorCode  Network::coordata() {
	int n;
	double x,y;
	n = Ntokens;

	/* Extract data from tokens */
	n = Ntokens;
	if (n < 3) return(SYNTAX_ERR_OF_COORDINATES); 
	if (!getfloat(Tok[1],&x)) return(ILLEGAL_NUMBER);
	if (!getfloat(Tok[2],&y)) return(ILLEGAL_NUMBER);


	/* Otherwise find node (and pattern) being referenced */
	HTIt jnode = Nht.find(Tok[0]);
	if (jnode == Nht.end())
		return(NODE_UNDEFINED);

	Node[jnode->second].x = x;
	Node[jnode->second].y = y;
	Node[jnode->second].coordFlag = 1;

	return OK;

}

Network::ErrorCode  Network::tankdata()
	/*
	**--------------------------------------------------------------
	**  Input:   none                                              
	**  Output:  returns error code                                
	**  Purpose: processes tank & reservoir data                   
	**  Format:                                                    
	**   [RESERVOIRS]                                            
	**     id elev (pattern)                                     
	**   [TANKS]                                                 
	**     id elev (pattern)                                     
	**     id elev initlevel minlevel maxlevel diam (minvol vcurve)                                    

	**--------------------------------------------------------------
	*/
{
	int   i;               /* Node index */
	int    n;              /* # data items */
	int p = 0;  /* Fixed grade time pattern index */
	double el        = 0.0; /* Elevation */


	/* Add new tank to data base */
	n = Ntokens;
	if (Ntanks == MaxTanks
		||  Nnodes == MaxNodes) return(SYNTAX_ERR_OF_TANK_OR_RESERV);
	Ntanks++;
	Nnodes++;
	i = MaxJuncs + Ntanks;                    /* i = node index.     */

	char*	tank_name = new char[MAX_ID];
	strncpy_s(tank_name, MAX_ID, Tok[0], MAX_ID-1);
	strncpy_s(Node[i].ID, MAX_ID, Tok[0], MAX_ID-1);

	if (!(Nht.insert(HTPair(tank_name, i))).second)  //id exists
		return(NODE_ID_EXISTS);

	/* Check for valid data */
	if (n < 2) return(SYNTAX_ERR_OF_TANK_OR_RESERV);                   /* Too few fields.   */
	if (!getfloat(Tok[1],&el)) return(ILLEGAL_NUMBER);   /* Read elevation    */

	STmplist* t;
	// tank level limits
	double initlevel=0;
	double minlevel=0;
	double maxlevel=0; 
	double diam=0;
	double minvol=0;
	int vcurve=0;  //volume curve index

	if (n <= 3) {                         /* Tank is reservoir.*/
		if (n == 3)  {                         /* Pattern supplied  */ 
			t = findID(Tok[2],Patlist);
			if (t == NULL) return(GRADE_PATTERN_UNDEFINED);
			p = t->i;
		}
	}  else { // is a tank
		if (n < 6) {
			return(SYNTAX_ERR_OF_TANK_OR_RESERV);              /* Too few fields for tank.*/
		} else  {
			/* Check for valid input data */
			if (!getfloat(Tok[2],&initlevel)) return(ILLEGAL_NUMBER);
			if (!getfloat(Tok[3],&minlevel))  return(ILLEGAL_NUMBER);
			if (!getfloat(Tok[4],&maxlevel))  return(ILLEGAL_NUMBER);
			if (!getfloat(Tok[5],&diam))      return(ILLEGAL_NUMBER);
			if (diam < 0.0)                   return(ILLEGAL_NUMBER);
			if (n >= 7
				&& !getfloat(Tok[6],&minvol))     return(ILLEGAL_NUMBER);

			/* If volume curve supplied check it exists */
			if (n == 8)  {                           
				t = findID(Tok[7],Curvelist);
				if (t == NULL) return(VOLUME_CURVE_UNDEFINED);
				vcurve = t->i;
			}
		}
	}

	Node[i].El            = el;               /* Elevation.           */
	Node[i].Ke           = 0;
	Tank[Ntanks].Node     = i;                /* Node index.          */
	Tank[Ntanks].H0       = initlevel;        /* Init. level.         */
	Tank[Ntanks].Hmin     = minlevel;         /* Min. level.          */
	Tank[Ntanks].Hmax     = maxlevel;         /* Max level.           */
	Tank[Ntanks].A        = diam;             /* Diameter.            */
	Tank[Ntanks].Pat      = p;                /* Fixed grade pattern. */
	/*
	*******************************************************************
	NOTE: The min, max, & initial volumes set here are based on a     
	nominal tank diameter. They will be modified in INPUT1.C if    
	a volume curve is supplied for this tank.                      
	*******************************************************************
	*/
	double area = _PI*diam*diam/4.0;
	Tank[Ntanks].Vmin = area*minlevel;
	if (minvol > 0.0) Tank[Ntanks].Vmin = minvol;
	Tank[Ntanks].V0 = Tank[Ntanks].Vmin + area*(initlevel - minlevel);
	Tank[Ntanks].Vmax = Tank[Ntanks].Vmin + area*(maxlevel - minlevel);

	Tank[Ntanks].Vcurve   = vcurve;           /* Volume curve         */
	return(OK);
}                        /* end of tankdata */

Network::ErrorCode Network::pipedata() 
	/*
	**--------------------------------------------------------------
	**  Input:   none                                              
	**  Output:  returns error code                                
	**  Purpose: processes pipe data                               
	**  Format:                                                    
	**    [PIPE]                                                
	**    id  node1  node2  length  diam  rcoeff (lcoeff) (status)          
	**--------------------------------------------------------------
	*/
{
	int   j1,                     /* Start-node index  */
		j2,                     /* End-node index    */
		n;                      /* # data items      */
	char  type = PIPE,            /* Link type         */
		status = OPEN;          /* Link status       */
	double length,                 /* Link length       */
		diam,                   /* Link diameter     */
		rcoeff,                 /* Roughness coeff.  */
		lcoeff = 0.0;           /* Minor loss coeff. */

	/* Add new pipe to data base */
	n = Ntokens;
	if (Nlinks == MaxLinks) return(SYNTAX_ERR_OF_PIPE);
	Npipes++;
	Nlinks++;

	//Prepare a hash table entry
	char* pipe_name = new char[MAX_ID];
	strncpy_s(pipe_name, MAX_ID, Tok[0], MAX_ID-1);
	strncpy_s(Link[Nlinks].ID, MAX_ID, Tok[0], MAX_ID-1);

	if (!(Lht.insert(HTPair(pipe_name, Nlinks))).second)  //id exists
		return(LINK_ID_EXISTS);

	/* Check for valid data */
	if (n < 6) return(SYNTAX_ERR_OF_PIPE);
	HTIt j1it = Nht.find(Tok[1]);
	HTIt j2it = Nht.find(Tok[2]);
	if (j1it == Nht.end() || j2it == Nht.end())
		return(NODE_UNDEFINED);
	j1 = j1it->second;  j2 = j2it->second;

	/*** Updated 10/25/00 ***/
	if (j1 == j2) return(LINK_IS_A_LOOP);    

	if (!getfloat(Tok[3],&length) ||
		!getfloat(Tok[4],&diam)   ||
		!getfloat(Tok[5],&rcoeff)
		) return(ILLEGAL_NUMBER);

	if (length <= 0.0 ||
		diam   <= 0.0 ||
		rcoeff <= 0.0
		) return(ILLEGAL_NUMBER);

	/* Case where either loss coeff. or status supplied */
	if (n == 7)
	{
		if      (match(Tok[6],"CV"))        {
			type = CV;
			++nCRSV ;
		}
		else if (match(Tok[6],"CLOSED"))    status = CLOSED;
		else if (match(Tok[6],"OPEN"))      status = OPEN;
		else if (!getfloat(Tok[6],&lcoeff)) return(ILLEGAL_NUMBER);
	}

	/* Case where both loss coeff. and status supplied */
	if (n == 8)
	{
		if (!getfloat(Tok[6],&lcoeff))   return(ILLEGAL_NUMBER);
		if      (match(Tok[7],"CV"))     {
			type = CV;
			++nCRSV;
		}
		else if (match(Tok[7],"CLOSED")) status = CLOSED;
		else if (match(Tok[7],"OPEN"))   status = OPEN;
		else return(ILLEGAL_NUMBER);
	}
	if (lcoeff < 0.0) return(ILLEGAL_NUMBER);

	/* Save pipe data */
	Link[Nlinks].N1    = j1;                  /* Start-node index */
	Link[Nlinks].N2    = j2;                  /* End-node index   */
	Link[Nlinks].Len   = length;              /* Length           */
	Link[Nlinks].Diam  = diam;                /* Diameter         */
	Link[Nlinks].Kc    = rcoeff;              /* Rough. coeff     */
	Link[Nlinks].Km    = lcoeff;              /* Loss coeff       */
	Link[Nlinks].Type  = type;                /* Link type        */
	Link[Nlinks].Stat  = status;              /* Link status      */

	return(OK);
}

Network::ErrorCode  Network::pumpdata() 
	/*
	**--------------------------------------------------------------
	** Input:   none                                                
	** Output:  returns error code                                  
	** Purpose: processes pump data                                 
	** Formats:                                                     
	**  [PUMP]                                                     
	**   (Version 1.x Format):                                              
	**   id  node1  node2  power                                   
	**   id  node1  node2  h1    q1                                
	**   id  node1  node2  h0    h1   q1   h2   q2                 
	**   (Version 2 Format):                                              
	**   id  node1  node2  KEYWORD value {KEYWORD value ...}       
	**   where KEYWORD = [POWER,HEAD,PATTERN,SPEED]                
	**--------------------------------------------------------------
	*/
{
	int   j,
		j1,                    /* Start-node index */
		j2,                    /* End-node index   */
		m, n;                  /* # data items     */
	double y;
	STmplist *t;                 /* Pattern record   */

	/* Add new pump to data base */
	n = Ntokens;
	if (Nlinks == MaxLinks ||
		Npumps == MaxPumps
		) return(SYNTAX_ERR_OF_PUMP);
	Nlinks++;
	Npumps++;
	++nPumpGPV;

	//Prepare a hash table entry
	char* pump_name = new char[MAX_ID];
	strncpy_s(pump_name, MAX_ID, Tok[0], MAX_ID-1);
	strncpy_s(Link[Nlinks].ID, MAX_ID, Tok[0], MAX_ID-1);
	if (!(Lht.insert(HTPair(pump_name, Nlinks))).second)  //id exists
		return(LINK_ID_EXISTS);

	/* Check for valid data */
	if (n < 4) return(SYNTAX_ERR_OF_PUMP);
	HTIt j1it = Nht.find(Tok[1]);
	HTIt j2it = Nht.find(Tok[2]);
	if (j1it == Nht.end() || j2it == Nht.end())
		return(NODE_UNDEFINED);
	j1 = j1it->second;  j2 = j2it->second;

	/*** Updated 10/25/00 ***/
	if (j1 == j2) return(LINK_IS_A_LOOP);    

	/* Save pump data */
	Link[Nlinks].N1    = j1;               /* Start-node index.  */
	Link[Nlinks].N2    = j2;               /* End-node index.    */
	Link[Nlinks].Diam  = Npumps;           /* Pump index.        */
	Link[Nlinks].Len   = 0.0;              /* Link length.       */

	Link[Nlinks].Type  = PUMP;             /* Link type.         */
	Link[Nlinks].Stat  = OPEN;             /* Link status.       */

	Pump[Npumps].Link = Nlinks;            /* Link index.        */
	Pump[Npumps].Ptype = NOCURVE;          /* Type of pump curve */
	Pump[Npumps].Hcurve = 0;               /* Pump curve index   */


	/* If 4-th token is a number then input follows Version 1.x format */
	/* so retrieve pump curve parameters */
	double X[MAX_TOKS];   // input buffer 
	if (getfloat(Tok[3],&X[0]))
	{
		m = 1;
		for (j=4; j<n; j++)
		{
			if (!getfloat(Tok[j],&X[m])) return(ILLEGAL_NUMBER);
			m++;
		}
		return(getpumpcurve(m, X));          /* Get pump curve params */
	}

	/* Otherwise input follows Version 2 format */
	/* so retrieve keyword/value pairs.         */
	m = 4;
	while (m < n)
	{
		if (match(Tok[m-1],"POWE"))          /* Const. HP curve       */
		{
			y = atof(Tok[m]);
			if (y <= 0.0) return(ILLEGAL_NUMBER);
			Pump[Npumps].Ptype = CONST_HP;
			Link[Nlinks].Km = y;
		}
		else if (match(Tok[m-1],"HEAD"))      /* Custom pump curve      */
		{
			t = findID(Tok[m],Curvelist);
			if (t == NULL) return(PUMP_CURVE_UNDEFINED);
			Pump[Npumps].Hcurve = t->i;
		}
		else if (match(Tok[m-1],"PATT"))   /* Speed/status pattern */
		{
			t = findID(Tok[m],Patlist);
			if (t == NULL) return(PUMP_PATTERN_UNDEFINED);
			Pump[Npumps].Upat = t->i;
		}
		else if (match(Tok[m-1],"SPEE"))     /* Speed setting */
		{
			if (!getfloat(Tok[m],&y)) return(ILLEGAL_NUMBER);
			if (y < 0.0) return(ILLEGAL_NUMBER);
			Link[Nlinks].Kc = y;
		}
		else return(SYNTAX_ERR_OF_PUMP);
		m = m + 2;                          /* Skip to next keyword token */
	}
	return(OK);
}

Network::ErrorCode  Network::valvedata() 
	/*
	**--------------------------------------------------------------
	**  Input:   none                                                
	**  Output:  returns error code                                  
	**  Purpose: processes valve data                                
	**  Format:                                                      
	**     [VALVE]                                                 
	**        id  node1  node2  diam  type  setting (lcoeff)       
	**--------------------------------------------------------------
	*/
{
	int   j1,                    /* Start-node index   */
		j2,                    /* End-node index     */
		n;                     /* # data items       */
	char  status = OPEN,       /* Valve status       */
		type;                  /* Valve type         */
	double diam = 0.0,            /* Valve diameter     */
		setting = 0,               /* Valve setting      */
		lcoeff = 0.0;          /* Minor loss coeff.  */
	STmplist *t;                 /* Curve record       */

	/* Add new valve to data base */
	n = Ntokens;
	if (Nlinks == MaxLinks ||
		Nvalves == MaxValves
		) return(SYNTAX_ERR_OF_VALVE);
	Nvalves++;
	Nlinks++;

	//Prepare a hash table entry
	char* valve_name = new char[MAX_ID];
	strncpy_s(valve_name, MAX_ID, Tok[0], MAX_ID-1);
	strncpy_s(Link[Nlinks].ID, MAX_ID, Tok[0], MAX_ID-1);
	if (!(Lht.insert(HTPair(valve_name, Nlinks))).second)  //id exists
		return(LINK_ID_EXISTS);

	/* Check for valid data */
	if (n < 6) return(SYNTAX_ERR_OF_VALVE);
	HTIt j1it = Nht.find(Tok[1]);
	HTIt j2it = Nht.find(Tok[2]);
	if (j1it == Nht.end() || j2it == Nht.end())
		return(NODE_UNDEFINED);
	j1 = j1it->second;  j2 = j2it->second;

	/*** Updated 10/25/00 ***/
	if (j1 == j2) return(LINK_IS_A_LOOP);    
	if      (match(Tok[4],"PRV")) {
		type = PRV;
		++nCRSV;
	}
	else if (match(Tok[4],"PSV")) {
		type = PSV;
		++nCRSV;
	}
	else if (match(Tok[4],"PBV")) {
		type = CV;
		++nCRSV;
		_lsWC.push_back(Warn(_lineno, PBV_SPECIFIED_IN_INP));
	}
	else if (match(Tok[4],"FCV")) {
		// Change all FCV to TCV
		type = TCV;
		_lsWC.push_back(Warn(_lineno, TCV_REPLACE_FCV));
		++nTFV;
	}
	else if (match(Tok[4],"TCV")) {
		type = TCV;
		++nTFV;
	}
	else if (match(Tok[4],"GPV")) {
		type = GPV;
		++nPumpGPV;
	}
	else    return(SYNTAX_ERR_OF_VALVE);    

	if (!getfloat(Tok[5],&setting)) return(SYNTAX_ERR_OF_VALVE);
	if (n >= 7 &&
		!getfloat(Tok[6],&lcoeff)
		) return(SYNTAX_ERR_OF_VALVE);
	/* Save valve data */
	Link[Nlinks].N1     = j1;                 /* Start-node index. */
	Link[Nlinks].N2     = j2;                 /* End-node index.   */
	Link[Nlinks].Diam   = diam;               /* Valve diameter.   */
	Link[Nlinks].Len    = 0.0;                /* Link length.      */
	Link[Nlinks].Kc     = setting;            /* Valve setting.    */
	Link[Nlinks].Km     = lcoeff;             /* Loss coeff        */
	Link[Nlinks].Kb     = 0.0;
	Link[Nlinks].Kw     = 0.0;
	Link[Nlinks].Type   = type;               /* Valve type.       */
	Link[Nlinks].Stat   = status;             /* Valve status.     */

	Valve[Nvalves].Link = Nlinks;             /* Link index.       */
	return(OK);
}
Network::ErrorCode  Network::patterndata()
	/*
	**--------------------------------------------------------------
	**  Input:   none                                                
	**  Output:  returns error code                                  
	**  Purpose: processes time pattern data                         
	**  Format:                                                      
	**     [PATTERNS]                                              
	**        id  mult1  mult2 .....                               
	**--------------------------------------------------------------
	*/
{
	int  i,n;
	double x;
	SFloatlist *f;
	STmplist   *p;
	n = Ntokens - 1;
	if (n < 1) return(SYNTAX_ERR_OF_PATTERN);            /* Too few values        */
	if (                               /* Check for new pattern */
		PrevPat != NULL &&
		strcmp(Tok[0],PrevPat->ID) == 0
		) p = PrevPat;
	else p = findID(Tok[0],Patlist);
	if (p == NULL) return(PATTERN_UNDEFINED);
	for (i=1; i<=n; i++)               /* Add multipliers to list */
	{
		if (!getfloat(Tok[i],&x)) return(ILLEGAL_NUMBER);
		f = (SFloatlist *) malloc(sizeof(SFloatlist));
		if (f == NULL) return(MALLOC_ERROR);
		f->value = x;
		f->next = p->x;
		p->x = f;
	}
	Pattern[p->i].Length += n;         /* Save # multipliers for pattern */
	PrevPat = p;                       /* Set previous pattern pointer */
	return(OK);
}                        /* end of patterndata */

Network::ErrorCode  Network::curvedata()
	/*
	**------------------------------------------------------
	**  Input:   none                                        
	**  Output:  returns error code                          
	**  Purpose: processes curve data                        
	**  Format:                                              
	**     [CURVES]                                        
	**      CurveID   x-value  y-value                    
	**------------------------------------------------------
	*/
{
	double      x,y;
	SFloatlist *fx, *fy;
	STmplist   *c;

	/* Check for valid curve ID */
	if (Ntokens < 3) return(SYNTAX_ERR_OF_CURVE);
	if (
		PrevCurve != NULL &&
		strcmp(Tok[0],PrevCurve->ID) == 0
		) c = PrevCurve;
	else c = findID(Tok[0],Curvelist);
	if (c == NULL) return(CURVE_UNDEFINED);

	/* Check for valid data */
	if (!getfloat(Tok[1],&x)) return(ILLEGAL_NUMBER);
	if (!getfloat(Tok[2],&y)) return(ILLEGAL_NUMBER);

	/* Add new data point to curve's linked list */
	fx = (SFloatlist *) malloc(sizeof(SFloatlist));
	fy = (SFloatlist *) malloc(sizeof(SFloatlist));
	if (fx == NULL || fy == NULL) return(MALLOC_ERROR);
	fx->value = x;
	fx->next = c->x;
	c->x = fx;
	fy->value = y;
	fy->next = c->y;
	c->y = fy;
	Curve[c->i].Npts++;

	/* Save the pointer to this curve */
	PrevCurve = c;
	return(OK);
}


Network::ErrorCode  Network::demanddata()
	/*
	**--------------------------------------------------------------
	**  Input:   none                                                
	**  Output:  returns error code                                  
	**  Purpose: processes node demand data                          
	**  Format:                                                      
	**     [DEMANDS]                                               
	**        MULTIPLY  factor                                     
	**        node  base_demand  (pattern)                         
	**
	**  NOTE: Demands entered in this section replace those 
	**        entered in the [JUNCTIONS] section
	**--------------------------------------------------------------
	*/
{
	int  j,n,p = 0;
	double y;
	Pdemand demand;
	STmplist *pat;

	/* Extract data from tokens */
	n = Ntokens;
	if (n < 2) return(SYNTAX_ERR_OF_DEMAND); 
	if (!getfloat(Tok[1],&y)) return(ILLEGAL_NUMBER);

	/* If MULTIPLY command, save multiplier */
	if (match(Tok[0],"MULTIPLY"))
	{
		if (y <= 0.0) return(ILLEGAL_NUMBER);
		else Dmult = y;
		return(OK);
	}

	/* Otherwise find node (and pattern) being referenced */
	HTIt jit = Nht.find(Tok[0]);
	if (jit == Nht.end() )       return(NODE_UNDEFINED);

	if (jit->second > Njuncs) return(NODE_UNDEFINED);
	j = jit->second;
	if (n >= 3)
	{
		pat = findID(Tok[2],Patlist);
		if (pat == NULL)  return(DEMAND_PATTERN_UNDEFINED);
		p = pat->i;
	}
	/* Replace any demand entered in [JUNCTIONS] section */
	/* (Such demand was temporarily stored in D[]) */

	/*** Updated 6/24/02 ***/
	demand = Node[j].D;
	if (demand)
	{
		demand->Base = y;
		demand->Pat  = p;
	}
	/*** End of update ***/

	/* Otherwise add a new demand to this junction */
	else
	{
		demand = (struct Sdemand *) malloc(sizeof(struct Sdemand));
		if (demand == NULL) return(MALLOC_ERROR);
		demand->Base = y;
		demand->Pat = p;
		demand->next = Node[j].D;
		Node[j].D = demand;
	}
	return(OK);
}                        /* end of demanddata */

Network::ErrorCode  Network::statusdata()
	/*
	**--------------------------------------------------------------
	**  Input:   none                                                
	**  Output:  returns error code                                  
	**  Purpose: processes link initial status data                          
	**  Formats:                                                     
	**    [STATUS]
	**       link   value
	**       link1  (link2)  value                                   
	**--------------------------------------------------------------
	*/
{
	int   j,n;
	long  i,i0,i1;
	double y = 0.0;
	char  status = OPEN;

	if (Nlinks == 0) return(LINK_UNDEFINED);
	n = Ntokens - 1;
	if (n < 1) return(SYNTAX_ERR_OF_STATUS);

	/* Check for legal status setting */
	if      (match(Tok[n],"OPEN"))    status = OPEN;
	else if (match(Tok[n],"CLOSED"))  status = CLOSED;
	else return (OK);  //just ignore it

	HTIt jit;

	/* Single link ID supplied */
	if (n == 1)
	{
		jit = Lht.find(Tok[0]);
		if ( jit == Lht.end()) return(OK); //do nothing
		j = jit->second;
		/* Cannot change status of a Check Valve */
		if (Link[j].Type == CV) return(SET_CV_STATUS);
	}

	/* Range of ID's supplied */
	else
	{
		/* Numerical range supplied */
		if ((i0 = atol(Tok[0])) > 0 && (i1 = atol(Tok[1])) > 0)
		{
			for (j=1; j<=Nlinks; j++)
			{
				i = atol(Link[j].ID);
				if (i >= i0 && i <= i1) Link[j].Stat = status;
			}
		}
		else
			for (j=1; j<=Nlinks; j++)
				if ( (strcmp(Tok[0],Link[j].ID) <= 0) &&
					(strcmp(Tok[1],Link[j].ID) >= 0)
					) Link[j].Stat = status;
	}
	return(OK);
}              /* end of statusdata */


Network::ErrorCode  Network::getpumpcurve(int n, double* X)
	/*
	**--------------------------------------------------------
	**  Input:   n = number of parameters for pump curve
	**  Output:  returns error code
	**  Purpose: processes pump curve data for Version 1.1-
	**           style input data
	**  Notes:
	**    1. Called by pumpdata() in INPUT3.C
	**    2. Current link index & pump index of pump being
	**       processed is found in global variables Nlinks
	**       and Npumps, respectively
	**    3. Curve data read from input line is found in
	**       global variables X[0],...X[n-1]
	**---------------------------------------------------------
	*/
{
	double a,b,c,h0,h1,h2,q1,q2;

	if (n == 1)                /* Const. HP curve       */
	{
		if (X[0] <= 0.0) return(ILLEGAL_NUMBER);
		Pump[Npumps].Ptype = CONST_HP;
		Link[Nlinks].Km = X[0];
	}
	else
	{
		if (n == 2)             /* Generic power curve   */
		{
			q1 = X[1];
			h1 = X[0];
			h0 = 1.33334*h1;
			q2 = 2.0*q1;
			h2 = 0.0;
		}
		else if (n >= 5)        /* 3-pt. power curve     */
		{
			h0 = X[0];
			h1 = X[1];
			q1 = X[2];
			h2 = X[3];
			q2 = X[4];
		}
		else return(ILLEGAL_NUMBER);
		Pump[Npumps].Ptype = POWER_FUNC;
		if (!powercurve(h0,h1,h2,q1,q2,&a,&b,&c)) return(VOLUME_CURVE_ERROR);
		Pump[Npumps].H0 = -a;
		Pump[Npumps].R  = -b;
		Pump[Npumps].N  = c;
		Pump[Npumps].Q0 = q1;
		Pump[Npumps].Qmax  = pow((-a/b),(1.0/c));
		Pump[Npumps].Hmax  = h0;
	}
	return(OK);
}

int  Network::powercurve(double h0, double h1, double h2, double q1,
						 double q2, double *a, double *b, double *c)
						 /*
						 **---------------------------------------------------------
						 **  Input:   h0 = shutoff head
						 **           h1 = design head
						 **           h2 = head at max. flow
						 **           q1 = design flow
						 **           q2 = max. flow
						 **  Output:  *a, *b, *c = pump curve coeffs. (H = a-bQ^c),
						 **           Returns 1 if sucessful, 0 otherwise.
						 **  Purpose: computes coeffs. for pump curve
						 **----------------------------------------------------------
						 */
{
	double h4,h5;
	if (
		h0      < _TINY ||
		h0 - h1 < _TINY ||
		h1 - h2 < _TINY ||
		q1      < _TINY ||
		q2 - q1 < _TINY
		) return(0);
	*a = h0;
	h4 = h0 - h1;
	h5 = h0 - h2;
	*c = log(h5/h4)/log(q2/q1);
	if (*c <= 0.0 || *c > 20.0) return(0);
	*b = -h4/pow(q1,*c);

	/*** Updated 6/24/02 ***/
	if (*b >= 0.0) return(0);

	return(1);
}

Network::ErrorCode  Network::addpattern(char *id) {
	/*
	**-------------------------------------------------------------
	**  Input:   id = pattern ID label
	**  Output:  returns error code 
	**  Purpose: adds a new pattern to the database
	**--------------------------------------------------------------
	*/
	STmplist *p;

	/* Check if ID is same as last one processed */
	if (Patlist != NULL && strcmp(id,Patlist->ID) == 0) 
		return(PATTERN_ID_EXISTS);

	/* Check that pattern was not already created */
	STmplist* itemfound = findID(id, Patlist);
	if (itemfound == NULL)    {
		/* Update pattern count & create new list element */
		(MaxPats)++;
		p = (STmplist *) malloc(sizeof(STmplist));
		if (p == NULL) return(MALLOC_ERROR);

		/* Initialize list element properties */
		p->i = MaxPats;
		strncpy(p->ID,id,MAX_ID);
		p->x = NULL;
		p->y = NULL;
		p->next = Patlist;
		Patlist = p;
		return OK;
	} else return(PATTERN_ID_EXISTS);
}

Network::ErrorCode  Network::addcurve(char *id)
	/*
	**-------------------------------------------------------------
	**  Input:   id = curve ID label
	**  Output:  returns error code
	**  Purpose: adds a new curve to the database
	**--------------------------------------------------------------
	*/
{
	STmplist *c;

	/* Check if ID is same as last one processed */
	if (Curvelist != NULL && strcmp(id,Curvelist->ID) == 0) 
		return(CURVE_ID_EXISTS);

	/* Check that curve was not already created */
	if (findID(id,Curvelist) == NULL)
	{

		/* Update curve count & create new list element */
		(MaxCurves)++;
		c = (STmplist *) malloc(sizeof(STmplist));
		if (c == NULL) return(MALLOC_ERROR);

		/* Initialize list element properties */
		else
		{
			c->i = MaxCurves;
			strncpy(c->ID,id,MAX_ID);
			c->x = NULL;
			c->y = NULL;
			c->next = Curvelist;
			Curvelist = c;
			return(OK);
		}
	} else return(CURVE_ID_EXISTS);

}


Network::STmplist* Network::findID(char *id, STmplist *list) {
	/*
	**-------------------------------------------------------------
	**  Input:   id = ID label
	**           list = pointer to head of a temporary list
	**  Output:  returns list item with requested ID label 
	**  Purpose: searches for item in temporary list
	**-------------------------------------------------------------
	*/
	STmplist *item = NULL;
	for (item = list; item != NULL; item = item->next)     {
		if (strcmp(item->ID,id) == 0)         {
			return item;
		}
	}
	return(NULL);
}


int  Network::findmatch(char *line, char *keyword[])
	/*
	**--------------------------------------------------------------
	**  Input:   *line      = line from input file
	**           *keyword[] = list of NULL terminated keywords
	**  Output:  returns index of matching keyword or
	**           -1 if no match found
	**  Purpose: determines which keyword appears on input line
	**--------------------------------------------------------------
	*/
{
	int i = 0;
	while (keyword[i] != NULL)
	{
		if (match(line,keyword[i])) return(i);
		i++;
	}
	return(-1);
}                        /* end of findmatch */

int  Network::match(char *str, char *substr)
	/*
	**--------------------------------------------------------------
	**  Input:   *str    = string being searched
	**           *substr = substring being searched for
	**  Output:  returns 1 if substr found in str, 0 if not
	**  Purpose: sees if substr matches any part of str
	**
	**      (Not case sensitive)
	**--------------------------------------------------------------
	*/
{
	int i,j;

	/*** Updated 9/7/00 ***/
	/* Fail if substring is empty */
	if (!substr[0]) return(0);

	/* Skip leading blanks of str. */
	for (i=0; str[i]; i++)
		if (str[i] != ' ') break;

	/* Check if substr matches remainder of str. */
	for (i=i,j=0; substr[j]; i++,j++)
		if (!str[i] || ((str[i] | 32) != (substr[j] | 32)))
			return(0);
	return(1);
}                        /* end of match */



void  Network::freeTmplist(STmplist *t)
	/*----------------------------------------------------------------
	**  Input:   t = pointer to start of a temporary list
	**  Output:  none
	**  Purpose: frees memory used for temporary storage
	**           of pattern & curve data
	**----------------------------------------------------------------
	*/
{
	STmplist   *tnext;
	while (t != NULL)
	{
		tnext = t->next;
		freeFloatlist(t->x);
		freeFloatlist(t->y);
		free(t);
		t = tnext;
	}
}


void  Network::freeFloatlist(SFloatlist *f)
	/*----------------------------------------------------------------
	**  Input:   f = pointer to start of list of floats
	**  Output:  none
	**  Purpose: frees memory used for storing list of floats
	**----------------------------------------------------------------
	*/
{
	SFloatlist *fnext;
	while (f != NULL)
	{
		fnext = f->next;
		free(f);
		f = fnext;
	}
}

int  Network::gettokens(char *s)
	/*
	**--------------------------------------------------------------
	**  Input:   *s = string to be tokenized
	**  Output:  returns number of tokens in s
	**  Purpose: scans string for tokens, saving pointers to them
	**           in module global variable Tok[]
	**
	** Tokens can be separated by the characters listed in SEPSTR
	** (spaces, tabs, newline, carriage return) which is defined
	** in TYPES.H. Text between quotes is treated as a single token.
	**--------------------------------------------------------------
	*/
{
	int  len, m, n;
	char *c;

	/* Begin with no tokens */
	for (n=0; n<MAX_TOKS; n++) Tok[n] = NULL;
	n = 0;

	/* Truncate s at start of comment */
	c = strchr(s,';');
	if (c) *c = '\0';
	len = (int)strlen(s);

	/* Scan s for tokens until nothing left */
	while (len > 0 && n < MAX_TOKS)
	{
		m = (int)strcspn(s," \t\n\r");          /* Find token length */
		len -= m+1;                     /* Update length of s */
		if (m == 0) s++;                /* No token found */
		else
		{
			if (*s == '"')               /* Token begins with quote */
			{
				s++;                      /* Start token after quote */
				m = (int)strcspn(s,"\"\n\r");  /* Find end quote (or EOL) */
			}                            
			s[m] = '\0';                 /* Null-terminate the token */
			Tok[n] = s;                  /* Save pointer to token */
			n++;                         /* Update token count */
			s += m+1;                    /* Begin next token */
		}
	}
	return(n);
}                        /* End of gettokens */

int  Network::getfloat(char *s, double *y)
	/*
	**-----------------------------------------------------------
	**  Input:   *s = character string
	**  Output:  *y = floating point number
	**           returns 1 if conversion successful, 0 if not
	**  Purpose: converts string to floating point number
	**-----------------------------------------------------------
	*/
{
	char *endptr;
	*y = (double) strtod(s,&endptr);
	if (*endptr > 0) return(0);
	return(1);
}

TCHAR* Network::typeStr(ComponentType ctype, TCHAR* str) {
    switch (ctype) {
	case JUNC:
        _tcscpy(str, TEXT("Junction"));
        break;
	case RESERV:
        _tcscpy(str, TEXT("Reservoir/water source"));
        break;
	case TANK:
        _tcscpy(str, TEXT("Water Storage/tank"));
        break;
	case CV:
        _tcscpy(str, TEXT("Check valve"));
        break;
	case PIPE:
        _tcscpy(str, TEXT("Pipe"));
        break;
	case PUMP:
        _tcscpy(str, TEXT("Pump"));
        break;
	case PRV:
        _tcscpy(str, TEXT("Pressure reduce valve"));
        break;
	case PSV:
        _tcscpy(str, TEXT("PRessure sustain valve"));
        break;
	case FCV:
        _tcscpy(str, TEXT("Flow control valve"));
        break;
	case TCV:
        _tcscpy(str, TEXT("Throttle control valve"));
        break;
	case GPV:
        _tcscpy(str, TEXT("General purpose valve"));
        break;
	default:
        _tcscpy(str, TEXT("Unknown component"));
	}
    return str;
}

// report network info.
void Network::report() {
	USES_CONVERSION;
	Network* s = singleton;
	if (s == NULL) {
		_tprintf(TEXT("Network not created.\n"));
		return;
	}

	if (s->_ecCur != OK) {
		_tprintf(TEXT("Error loading inp file (%s): %s in Line %d\n"),
			s->_inpfilename,
			s->problemText(s->_ecCur), s->_lineEc);
		return;
	}

	_tprintf(TEXT( "EPANET NETWORK SUMMARY:\n")
		TEXT(" Total nodes: \t %d, among which: \n")
		TEXT("   Junctions: %d [Index 1 ~ %d] \n")
		TEXT("     Junc. w/ emitter (water leaks): %d), \n")
		TEXT("     Junc. w/ demand (water users): %d), \n")
		TEXT("     Tanks/Reservoirs: %d [Index %d ~ %d], listed as: \n")
		, s->MaxNodes, s->MaxJuncs, s->MaxJuncs, s->Nemitters,  s->Nusers,
		s->MaxTanks, s->MaxJuncs+1, s->MaxNodes);

	TCHAR utxt[64]; // unit text 

	// Tank list
	for (int iTank = 1; iTank <= s->MaxTanks; ++ iTank) {
		Stank aTank = s->Tank[iTank];
		int inode = aTank.Node; 
		Snode aNode = s->Node[inode];
		unitText(ELEV, s->Unitsflag, utxt);
		_tprintf(TEXT("        Tank @ Node %d <%s>: ")
			TEXT("Elev: %6.2f %s, Level Range: %6.2f to %6.2f %s\n"),
			inode, A2T(aNode.ID), 
			aNode.El * s->Ucf[ELEV], utxt, 
			aTank.Hmin * s->Ucf[ELEV], 
			aTank.Hmax * s->Ucf[ELEV], utxt
			);
	}
	_tprintf(TEXT("  Channel [L] needed for each tank/reservoir.\n"));

	_tprintf(TEXT( " Total links: \t %d, among which: \n")
		TEXT("   Pipes: %d [Index 1 ~ %d] \n")
		TEXT("   Check Valves, PSVs and PRVs listed as: \n"),
		s->MaxLinks, s->MaxPipes, s->MaxPipes);

	// CRSV list

	for (int iCRSV = 0; iCRSV<s->nCRSV; ++iCRSV) {
		int linkno = s->tCRSV[iCRSV];
		Slink* aLink = &(s->Link[linkno]);
		switch (aLink->Type) {
		case CV:
			_tprintf(TEXT("     Check Valve @ Link %d <%s>, Node %d <%s> to %d <%s>\n"), 
				linkno, A2T(aLink->ID), 
				aLink->N1, A2T(s->Node[aLink->N1].ID),
				aLink->N2, A2T(s->Node[aLink->N2].ID)
				);
			break;

		case PSV:
			unitText(PRESSURE, s->Pressflag, utxt);
			_tprintf(TEXT( "    PSV @ Link %d <%s>: ")
				TEXT("Set Junc. %d's Pres. > %6.2f %s when possible.\n"),
				linkno, A2T(aLink->ID),
				aLink->N1, aLink->Kc *s->Ucf[PRESSURE], utxt
				);
			break;

		case PRV:
			unitText(PRESSURE, s->Pressflag, utxt);
			_tprintf(TEXT( "     PSV @ Link %d <%s>: ")
				TEXT("Set Junc. %d's Pres. < %6.2f %s when possible.\n"),
				linkno, A2T(aLink->ID),
				aLink->N2, aLink->Kc *s->Ucf[PRESSURE], utxt
				);
			break;
		}
	}

	if (s->nCRSV) _tprintf(TEXT("   Pumps and GPVs listed as: \n"));

	// Pump/GPV list
	for (int iPG = 0; iPG<s->nPumpGPV; ++iPG) {
		int linkno = s->tPumpGPV[iPG];
		Slink* aLink = &(s->Link[linkno]);
		switch (aLink->Type) {
		case PUMP:
			//			int pumpindex = int(aLink->Diam); // this is what epanet does this
			_tprintf(TEXT("     Pump @ Link %d <%s>, Node %d <%s> to %d <%s>\n"), 
				linkno, A2T(aLink->ID), 
				aLink->N1, A2T(s->Node[aLink->N1].ID),
				aLink->N2, A2T(s->Node[aLink->N2].ID)
				);
			break;
		case GPV:
			_tprintf(TEXT("     GPV @ Link %d <%s>, Node %d <%s> to %d <%s>\n"), 
				linkno, A2T(aLink->ID), 
				aLink->N1, A2T(s->Node[aLink->N1].ID),
				aLink->N2, A2T(s->Node[aLink->N2].ID)
				);		
			break;
		}
	}

	if (s->nCRSV) _tprintf(TEXT("   Channel [F] and [B] needed for each pump/GPV.\n"));

	// TCV/FCV list
	for (int iTFV = 0; iTFV< s->nTFV; ++iTFV) {
		int linkno = s->tTFV[iTFV];
		Slink* aLink = &(s->Link[linkno]);
		unitText(FLOW, s->Flowflag, utxt);

		switch (aLink->Type) {
		case TCV:
			_tprintf(TEXT("     TCV @ Link %d <%s>:")
				TEXT("minor headloss adjustable.\n"), linkno, A2T(aLink->ID));
			break;
		case FCV:
			_tprintf(TEXT("     FCV @ Link %d <%s>:")
				TEXT("Adjust minor headloss so flowrate maintains %6.2f %s.\n"), 
				linkno, A2T(aLink->ID),
				aLink->Kc * s->Ucf[FLOW], utxt);
			break;
		}
	}

	if (s->nTFV) _tprintf(TEXT("   Channel [V] needed for each TCV/FCV.\n"));

	// report warning information
	s->reportWarns();

}

void Network::reportWarns() {
	TCHAR warnTxt[1024];
	std::list<Warn>::iterator iWarn;

	_tprintf(TEXT("Network Warnings: "));
	for (iWarn=_lsWC.begin(); iWarn!=_lsWC.end(); ++iWarn) {

		_tprintf(TEXT("Line %d: %s \n"), iWarn->first, problemText(iWarn->second));
	}
    _tprintf(TEXT("\n"));
	return;
}
//		"Network Info:\n\tTotal Nodes: %d, Tanks+Reservoirs: %d, Junctions: %d\n\tTotal Links: %d, Pipes: %d, Pumps: %d, Valves: %d\n"),
//		MaxNodes, MaxTanks, MaxJuncs, MaxLinks, MaxPipes, MaxPumps, MaxValves); 


Network::ErrorCode Network::unitText(FieldType type, UnitsType unit, TCHAR* str, AbbrLev al) {
	ErrorCode ec = OK;
	switch (type) {
	case ELEV: case HEAD: case LENGTH: case HEADLOSS:
		switch (unit) {
		case SI: 
			switch (al) {
			case UNIT_SHORT: _tcscpy(str, TEXT("m")); break;
			case UNIT_FULL: _tcscpy(str, TEXT("Meter(s)")); break;
			}
			break;
		case US:
			switch (al) {
			case UNIT_SHORT: _tcscpy(str, TEXT("ft")); break;
			case UNIT_FULL: _tcscpy(str, TEXT("Feet")); break;
			}
			break;
		default: // unkown unit system
			ec = UNKNOWN_UNIT;
		}
		break;

	case DIAM:
		switch (unit) {
		case SI:
			switch (al) {
			case UNIT_SHORT: _tcscpy(str, TEXT("mm")); break;
			case UNIT_FULL: _tcscpy(str, TEXT("Millimeter(s)")); break;
			}
			break;
		case US:
			switch (al) {
			case UNIT_SHORT: _tcscpy(str, TEXT("in")); break;
			case UNIT_FULL: _tcscpy(str, TEXT("Inch(es)")); break;
			}
			break;
		default:
			ec = UNKNOWN_UNIT;
		}
		break;

	case VELOCITY:
		switch (unit) {
		case SI:
			switch (al) {
			case UNIT_SHORT: _tcscpy(str, TEXT("m/s")); break;
			case UNIT_FULL: _tcscpy(str, TEXT("Meters per second")); break;
			}
			break;
		case US:
			switch (al) {
			case UNIT_SHORT: _tcscpy(str, TEXT("fps")); break;
			case UNIT_FULL: _tcscpy(str, TEXT("Feet per second")); break;
			}
			break;
		default:
			ec = UNKNOWN_UNIT;
		}
		break;

	case LSTATUS: case FRICTION: // dimensionless
		_tcscpy(str, TEXT(""));
		break;

	default:// Unknown type of variable
		ec = UNKNOWN_UNIT;
		_tcscpy(str, TEXT(""));

	}
	return ec;
}

Network::ErrorCode Network::unitText(FieldType type, FlowUnitsType unit, TCHAR* str, AbbrLev al) {
	if (type != DEMAND && type != FLOW ) {
		_tcscpy(str, TEXT(""));
		return UNKNOWN_UNIT;
	}

	ErrorCode ec = OK;	
	switch (unit) {
	case CFS:
		switch (al) {
		case UNIT_SHORT: _tcscpy(str, TEXT("cfs")); break;
		case UNIT_FULL: _tcscpy(str, TEXT("Cubic feet per second")); break;
		} 
		break;
	case GPM: 
		switch (al) {
		case UNIT_SHORT: _tcscpy(str, TEXT("gpm")); break;
		case UNIT_FULL: _tcscpy(str, TEXT("Gallons per minute")); break;
		}
		break;
	case MGD:
		switch (al) {
		case UNIT_SHORT: _tcscpy(str, TEXT("mgd")); break;
		case UNIT_FULL: _tcscpy(str, TEXT("Million gallons per day")); break;
		}
		break;
	case IMGD:
		switch (al) {
		case UNIT_SHORT: _tcscpy(str, TEXT("imgd")); break;
		case UNIT_FULL: _tcscpy(str, TEXT("Imperial miilion gallons per day")); break;
		}
		break;

	case AFD:
		switch (al) {
		case UNIT_SHORT: _tcscpy(str, TEXT("afd")); break;
		case UNIT_FULL: _tcscpy(str, TEXT("Acre-feet per day")); break;
		}
		break;

	case LPS:
		switch (al) {
		case UNIT_SHORT: _tcscpy(str, TEXT("L/s")); break;
		case UNIT_FULL: _tcscpy(str, TEXT("Liters per second")); break;
		}
		break;

	case LPM:
		switch (al) {
		case UNIT_SHORT: _tcscpy(str, TEXT("L/m")); break;
		case UNIT_FULL: _tcscpy(str, TEXT("Liters per second")); break;
		}
		break;

	case MLD:
		switch (al) {
		case UNIT_SHORT: _tcscpy(str, TEXT("*1e6 L/d")); break;
		case UNIT_FULL: _tcscpy(str, TEXT("Million liters per day")); break;
		}
		break;

	case CMH:
		switch (al) {
		case UNIT_SHORT: _tcscpy(str, TEXT("m3/h")); break;
		case UNIT_FULL: _tcscpy(str, TEXT("Cublic meters per hour")); break;
		}
		break;

	case CMD:
		switch (al) {
		case UNIT_SHORT: _tcscpy(str, TEXT("m3/d")); break;
		case UNIT_FULL: _tcscpy(str, TEXT("Cubic meters per day")); break;
		}
		break;

	default:
		_tcscpy(str, TEXT(""));
		ec = UNKNOWN_UNIT;
	}
	return ec;
}

Network::ErrorCode Network::unitText(FieldType type, PressUnitsType unit, TCHAR* str, AbbrLev al) {
	if (type != PRESSURE ) {
		_tcscpy(str, TEXT(""));
		return UNKNOWN_UNIT;
	}

	ErrorCode ec = OK;	
	switch (unit) {
	case PSI:
		switch (al) {
		case UNIT_SHORT: _tcscpy(str, TEXT("psi")); break;
		case UNIT_FULL: _tcscpy(str, TEXT("Pounds per square inch")); break;
		}
		break;
	case KPA:
		switch (al) {
		case UNIT_SHORT: _tcscpy(str, TEXT("kPa")); break;
		case UNIT_FULL: _tcscpy(str, TEXT("Kilo Pascals")); break;
		}
		break;
	case METERS:
		switch (al) {
		case UNIT_SHORT: _tcscpy(str, TEXT("m")); break;
		case UNIT_FULL: _tcscpy(str, TEXT("Meter(s)")); break;
		}
		break;
	default:
		_tcscpy(str, TEXT(""));
		ec = UNKNOWN_UNIT;
	}
	return ec;
}



