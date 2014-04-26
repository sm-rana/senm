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

#include <Windows.h>
#include "Network.h"
#include "Solver.h"


Solver::Solver()  {//setdefaults() in input1.c 
	Htol = 0.005;   /* Default head tolerance         */
    Qtol = 1e-4;
	Hacc      = 0.001;            /* Default hydraulic accuracy     */
	MaxIter   = 200;         /* Default max. hydraulic trials  */
	ExtraIter = -1;              /* Stop if network unbalanced     */
	RQtol     = 1E-7;           /* Default hydraulics parameters  */
//	DampLimit = 0; 
    
}

Solver::~Solver() {
	// destroy 
	free(D);
	free(H);
	free(_tabXd);
	free(Aii);
	free(Aij);
	free(F);
	free(E);
	free(P);
	free(Y);
	free(X);
	free(Q);
	free(K);
	free(S);

	free(B);

}


Solver::EWICode Solver::createSolver(
	Network* net, DataSource* ds, unsigned thdN, Solver** outsolver) {
        EWICode err = OK;
        Solver* prec = NULL;

		// check network
		if (net == NULL) {
            reportEWI(err = NETWORK_NOT_EXIST, thdN); goto END;
		} else if (net->_ecCur != OK) {
            reportEWI(err = NETWORK_NOT_READY, thdN); goto END;
		}
		// check datasource
		if (ds == NULL) {
			reportEWI(err = DS_NOT_EXIST, thdN); goto END;
		} 

		int nChan = ds->n_chan;

		if (nChan == 0) { //no channels
			reportEWI(NO_CHANNELS); 
		}

		// load channel availability table
		int N = net->MaxNodes;
		int M = net->MaxLinks;

		Channel::Type *tabCAT = new Channel::Type[M+N+1];
		for (int iTab = 0; iTab<M+N+1; ++iTab) tabCAT[iTab]=Channel::NONE;

		Channel * lsChan = ds->lsChan;
		int nCd = 0; // number of channel [D]

		for (Channel* iChan = lsChan; iChan; iChan = iChan->next) {
			// duplicated channels
			int idx = iChan->mindex;
            int chid = iChan->key;
			if (unsigned(iChan->type) & unsigned(tabCAT[idx])) {
				reportEWI(DUP_CHANNELS, thdN);
				continue;
			}

			switch (iChan->type) {
			case Channel::D: case Channel::P:
				if (idx <=0 || idx > net->MaxJuncs) {
                    reportEWI(err = CHAN_MINDEX_ERR, Network::JUNC, idx, chid, thdN);
				}
				break;
			case Channel::L: 
				if (idx <= net->MaxJuncs || idx > net->MaxNodes) {
					reportEWI(err = CHAN_MINDEX_ERR, Network::TANK, idx, chid, thdN);
				}
				break;
			case Channel::V: 
				if (idx <= net->MaxPipes + net->MaxPumps 
					|| idx > net->MaxLinks) {
                    reportEWI(err = CHAN_MINDEX_ERR, Network::TCV, idx, chid, thdN);
				}
				break;
            case Channel::A: case Channel::B: 
				if (idx <=0 || idx > net->MaxNodes) {
                    reportEWI(err = CHAN_MINDEX_ERR, Network::JUNC, idx, chid, thdN);
				}
                break;
			case Channel::F: 
				if (idx <= net->MaxPipes || 
					idx > net->MaxPipes + net->MaxPumps) {
                    reportEWI(err = CHAN_MINDEX_ERR, Network::PUMP, idx, chid, thdN);
				}
				break;
			case Channel::C: case Channel::Q: 
				if (idx <= 0 || 
					idx > net->MaxLinks) {
                    reportEWI(err = CHAN_MINDEX_ERR, Network::PIPE, idx, chid, thdN);
				}
				break; 
			}
			if (err != OK) goto END;

			switch (iChan->type) {
			case Channel::NONE:
				break; //do nothing
			case Channel::D: //node
				++ nCd;  // go through
			case Channel::L: case Channel::P: case Channel::B: case Channel::A:
				tabCAT[idx] = Channel::Type(
					unsigned(tabCAT[idx]) | unsigned(iChan->type));
				break;
			default: //link
				tabCAT[idx+N] = Channel::Type(
					unsigned(tabCAT[idx+N]) | unsigned(iChan->type));
			}
		}

		//Check data source /channel list integrity
		// 1. Tank/Reservoir must have [L] channel,
		for (int iTank = 1; iTank <= net->MaxTanks; ++ iTank) {
			int tpNode = net->Tank[iTank].Node;
			if (unsigned(tabCAT[tpNode]) & unsigned(Channel::D) ) {
				reportEWI(D_CHANNEL_AT_TANK, Network::TANK, tpNode, thdN);
				--nCd;
			}
			if ((unsigned(tabCAT[tpNode]) & unsigned(Channel::L)) == 0) {
                reportEWI(err = NO_L_TANK, Network::TANK, tpNode, thdN);
                goto END;
			}
		}

		// 2. Pump/GPV must have [F] and [B] channel,
		for (int iPG = 0; iPG< net->nPumpGPV; ++iPG) {
			int tpLink = net->tPumpGPV[iPG];
            int tpN2 = net->Link[tpLink].N2;
			if ( (unsigned(tabCAT[N+tpLink]) & unsigned(Channel::F)) ==0 ) {
                reportEWI(err = NO_F_PUMP_GPV, Network::PUMP, tpLink, thdN);
                goto END;
			}
			if ( (unsigned(tabCAT[tpN2]) & unsigned(Channel::B)) ==0 ) {
                reportEWI(err = NO_B_PUMP_GPV, Network::PUMP, tpLink, thdN);
                goto END;
			}
            if (net->Node[net->Link[tpLink].N2].Ke != 0) {
				reportEWI(B_EMITTER_CONFLICT, Network::PUMP, tpLink, thdN);
			}

		}
		// 3. TCV/TFV must have [V] Channel
		for (int iTFV = 0; iTFV<net->nTFV; ++iTFV) {
			int tpLink = net->tTFV[iTFV];
			if ( (unsigned(tabCAT[N+tpLink]) & unsigned(Channel::V)) ==0 ) {
                reportEWI(err = NO_V_TCV_FCV, Network::TCV, tpLink, thdN); goto END;
			}
		}

		// Build a solver
		prec = new Solver(); //solver precursor
		prec->_net = net;
		prec->_ds = ds;

		prec->N = N; 
		prec->M = M; 
		prec->_tabCAT = tabCAT;
		prec->_nChan = nChan;
		prec->_lsChan = lsChan;
		prec->_nXd = net->MaxJuncs - nCd;


		int n = prec->N + 1;
		int m = prec->M + 1;

		prec->D  = (double *) calloc(n, sizeof(double));
		prec->H =  (double *) calloc(n, sizeof(double));

		prec->_tabXd = (int*) calloc(prec->_nXd, sizeof(int)); 
		prec->_tabIdx2Xd = (int*) calloc(net->MaxJuncs + 1, sizeof(int));
		prec->Aii = (double *) calloc(n,sizeof(double));
		prec->Aij = (double *) calloc(net->Ncoeffs+1,sizeof(double));
		prec->F   = (double *) calloc(n,sizeof(double));
		prec->E   = (double *) calloc(n,sizeof(double));
		prec->P   = (double *) calloc(m,sizeof(double));
		prec->Y   = (double *) calloc(m,sizeof(double));
		prec->X   = (double *) calloc(m>n?m:n,sizeof(double));


		prec->Q    = (double *) calloc(m, sizeof(double));
		prec->K    = (double *) calloc(m, sizeof(double));
		prec->S    = (LinkStatus  *) calloc(m, sizeof(LinkStatus));

		prec->B = (double*) calloc(n, sizeof(double));

		if (prec->D==NULL || prec->H==NULL || prec->Q==NULL || 
			prec->K==NULL || prec->S==NULL) {
                reportEWI(err = MALLOC_ERROR, thdN); goto END;
		}

		// Prepare stochastic demand (xd) indexing array
		int iXd = 0;
		for (int iJunc = 1; iJunc <= net->MaxJuncs; ++ iJunc) {
			Network::Pdemand	pd = net->Node[iJunc].D;

			if (pd!=NULL && pd->Base>0 &&
                 ((unsigned(tabCAT[iJunc]) & unsigned(Channel::D)) ==0 )  )
			{// is a water user and no Channel D
				prec->_tabXd[iXd] = iJunc;
				prec->_tabIdx2Xd[iJunc] = iXd;
				++iXd;
			} else {
				prec->_tabIdx2Xd[iJunc] = -1; // no xd entry
			}
		}			
        // load link availablity
		for (int iLink = 1; iLink <= net->MaxLinks; ++ iLink) {
			if (net->Link[iLink].Stat == Network::OPEN) {
                //enable
                prec->S[iLink] = FULL_OPEN;
			} else {
                prec->S[iLink] = DISABLED; //disable
			}
		}


		// get pointer to snapshot data
		prec->_ss = ds->_datBuf;

		/* Set exponent in head loss equation and adjust flow-resistance tolerance.*/
		if (net->Formflag == Network::HW) prec->Hexp = 1.852;
		else                prec->Hexp = 2.0;

        // set thread number
		prec->threadN = thdN;


		*outsolver = prec;
		return OK;

END:
        if (tabCAT!=NULL) delete[] tabCAT;
        if (prec!=NULL) delete prec;
        *outsolver = NULL;
        return err;
}

//void Solver::fillEWinfo(Error err, TCHAR* txt) {

void Solver::reportEWI(EWICode code, unsigned thdN) {
	return reportEWI(code, Network::UNDETERMINED, 0, 0, thdN);
}
void Solver::reportEWI(EWICode code, Network::ComponentType ctype, int id1, unsigned thdN) {
	return reportEWI(code, Network::UNDETERMINED, id1, 0, thdN);
}
void Solver::reportEWI(EWICode code, Network::ComponentType ctype, int id1, int id2, unsigned thdN){
	TCHAR texts[LAST_DUMMY_EWI][MAX_ERR_STRING_SIZE] = {
		TEXT("OK."), 

		TEXT("The Network is invalid. Load inp file to the network first."),
		TEXT("Can not allocate memory for Solver."), 
		TEXT("The network contains errors."), 
		TEXT("The data source does not exist or is not ready."), 

		TEXT("Tank %d does not have a [L] channel for its current water level."), 
		TEXT("Pump or GPV %d does not have a [F] channel of current flowrate."), 
		TEXT("Pump or GPV %d does not have a [B] channel showing")
                TEXT("discharge-side pressure or pressure head."  ), 
		TEXT("TCV or TFV %d does not have a [V] channel showing")
                TEXT("minor headloss setting."),
        TEXT("Could not find the network index %d for Channel %d."), 

 		TEXT("No channels loaded."),
		TEXT("Multiple channels refer to the same measurement."),
		TEXT("A [D] channel is assigned to a tank/reservoir (index %d)."),
		TEXT("An FCV (index %d) is in the network. will be replaced by a TCV.")
                TEXT("(i.e., changeable minor headloss)"),

		TEXT("Not enough realizations of stochastic demands provided to the ")
                TEXT("solver, will use baseline demand instead."),
		TEXT("Too many demand realizations provided to the solver. only the ")
                TEXT("first nxd of them will be used."),
        TEXT("Channel B (index %d) in the node conflicts with an emitter setting."),
        TEXT("Hydraulic equation can not be solved due to ill-conditions around node %d."),

        TEXT("PRV/PSV causing ill conditionality in the hydraulic equations."),
        TEXT("Hydraulic equations are unable to solve."),
        TEXT("Maximum number of iterations reached. cannot find a solution."),

        TEXT("Potential CV/PSV/PRV cyclic problem.")
	};

	TCHAR txt[MAX_ERR_STRING_SIZE+MAX_ERR_PREFIX_SIZE];
	TCHAR tmptxt[MAX_ERR_STRING_SIZE+MAX_ERR_PREFIX_SIZE];
//	TCHAR ctype_str[MAX_NET_CTYPE_STR_SIZE];

	switch (code) {
	case NO_L_TANK: case NO_F_PUMP_GPV: case NO_B_PUMP_GPV: case NO_V_TCV_FCV:
	case D_CHANNEL_AT_TANK: case FCV_REPLACED_BY_TCV: case B_EMITTER_CONFLICT:
	case EQN_ILL_COND:
        _stprintf(tmptxt, texts[code], id1);
        break;
	case CHAN_MINDEX_ERR:
        _stprintf(tmptxt, texts[code], id1, id2);
        break;
	default:
        _stprintf(tmptxt, texts[code]);
	}

    _stprintf(txt, TEXT("Solver:%s"), tmptxt);

    return ewi(txt, thdN);

}

//
//Solver::ECode Solver::setSnapshot(double* shot) {
//	// Use the info in the (data) snapshot to set 
//	// (1) link open/close status
//	// (2) tank/reservoir level
//	// (3) Known water demands
//	return OK;
//}

void Solver::setLFBVCD() {
	int iiCh = 0;
    int nN1 = 0;  int nN2 = 0;
	for (Channel* iCh = _lsChan; iCh; iCh = iCh->next, ++iiCh) {
		switch (iCh->type) {
		case Channel::L: //set tank levels (alt.)
			H[iCh->mindex] = _ss[iiCh]/_net->Ucf[Network::HEAD];
			break;
		case Channel::F:  
			// cut the link, 
			// add inflow to discharge side
			// add outflow to suction side
			S[iCh->mindex] = DISABLED; //disable
			nN1 = _net->Link[iCh->mindex].N1;
			nN2 = _net->Link[iCh->mindex].N2;
			D[nN1] += _ss[iiCh]/_net->Ucf[Network::FLOW];
			D[nN2] -= _ss[iiCh]/_net->Ucf[Network::FLOW];
			break;
		case Channel::B:
			// Store data into B temporarily
			B[_net->Link[iCh->mindex].N2] = 
				_ss[iiCh]/_net->Ucf[Network::PRESSURE] + 
				_net->Node[iCh->mindex].El;
			break;
		case Channel::V:
			if (_ss[iiCh] == -1) { //FCV or TCV closed
				S[iCh->mindex] = DISABLED;
				K[iCh->mindex] = CBIG;
			} else { //open, set minor head loss coeff
				S[iCh->mindex] = PARTIAL;
				K[iCh->mindex] = _ss[iiCh];
			}
			break;
		case Channel::C:
			//set link status
			if (_ss[iiCh] == 1) { // open
				S[iCh->mindex] = DISABLED;
			} else {
				S[iCh->mindex] = FULL_OPEN;
			}
			break;
		case Channel::D:
			// add demand
			D[iCh->mindex] += _ss[iiCh]/_net->Ucf[Network::FLOW];
			break;
		}
	}
}

void Solver::run(double *xd, int nXd) {
    // reset demands
	for (int idm = 0; idm < _net->Nnodes+1; ++idm) {
        D[idm] = 0.0;
	}

	// set xd using demands 
	if (nXd < _nXd) { // less water demands provided
		for (int ixd = 0; ixd < nXd; ++ixd) {
			D[_tabXd[ixd]] += xd[ixd];
		}
        reportEWI(NOT_ENOUGH_XD, threadN);
	} else {
		for (int ixd = 0; ixd < this->_nXd; ++ixd) {
			D[_tabXd[ixd]] += xd[ixd];
		}
		if (nXd > _nXd) { // more xd provided
            reportEWI(TOO_MANY_XD, threadN);
		}
	}

	// other initialization in inithyd() in hydraul.c
	// init link flows
	for (int i=1; i<M; ++i) {
		/* Initialize status and setting */
		K[i] = _net->Link[i].Kc;

		/* Start active control valves in ACTIVE position,
		unless pressure setting not specified*/                     
		if ( (_net->Link[i].Type == Network::PRV || 
			_net->Link[i].Type == Network::PSV )
			&& (_net->Link[i].Kc != 0) ) S[i] = PARTIAL;                                     

		/* Initialize flows if necessary */
		if (S[i] == CLOSED || S[i] == DISABLED) Q[i] = CSMALL;
		else if (abs(Q[i]) <= CSMALL)
			// init to 1 fps for regular pipes
			Q[i] = _PI*pow(_net->Link[i].Diam, 2)/4.0;
	}

	// set channel data for Channels [L][F][B][V][C][D]
	setLFBVCD();

	/* Initialize emitter flows to 1 fps */
	for (int i=1; i<=_net->Njuncs; i++)
		if (_net->Node[i].Ke > 0.0 ||  // emitter
            unsigned(_tabCAT[i]) & unsigned(Channel::B)) // B channel
			  E[i] = 1.0;

	// run hydraulic simulation, netsolve() in hydraul.c
    int probNode; // link causing ill-conditionality
    double newerr = 0;  // error between subsequent iterations
    int statChange;  // valve change flag
    int iter = 0;

    // iterative solver
    do {
        
        newcoeffs(); //update P, Y, A, F for each component

		probNode = linsolve(_net->Njuncs, Aii, Aij, F); // sparse linear solver

        if (probNode>0) { // check for ill-conditionality
            reportEWI(EQN_ILL_COND, Network::UNDETERMINED, probNode, threadN);
            if (badvalve(_net->Order[probNode])) { // a valve caused it
                reportEWI(VALVE_CAUSE_ILL_COND, threadN);
                continue;
			} else {// trouble in lin solver 
				reportEWI(UNABLE_TO_SOLVE, threadN);
               break;
			}
		}

        // update heads and flows
		for (int i = 1; i<=_net->Njuncs; ++i) H[i] = F[_net->Row[i]];
        newerr = newflows();

        // update cv, prv, psv status
        //statChange = valvestatus() && linkstatus(); 

        ++iter;
        if (iter > MaxIter ) {
			reportEWI(MAX_ITER_REACHED);
            if (statChange) 
				reportEWI(CV_PSV_PRV_PROB);
            break;
		}
	} while (newerr > Hacc);

   /* Add any emitter flows to junction demands */
   for (int i=1; i<=_net->Njuncs; i++) D[i] += E[i];
}

int  Solver::linkstatus()
/*
**--------------------------------------------------------------
**  Input:   none                                                
**  Output:  returns 1 if any link changes status, 0 otherwise   
**  Purpose: determines new status for CVs                  
**--------------------------------------------------------------
*/
{
   int   change = FALSE,             /* Status change flag      */
         k,                          /* Link index              */
         n1,                         /* Start node index        */
         n2;                         /* End node index          */
   double dh;                        /* Head difference         */
   LinkStatus  status;                     /* Current status          */

   /* Examine each CV */
   for (k=1; k<=_net->Npipes; k++)
   {
      n1 = _net->Link[k].N1;
      n2 = _net->Link[k].N2;
      dh = H[n1] - H[n2];

      /* Re-open temporarily closed links (status = XHEAD or TEMPCLOSED) */
      status = S[k];
      if (status != DISABLED) S[k] = FULL_OPEN;

      /* Check for status changes in CVs and pumps */
      if (_net->Link[k].Type == Network::CV) S[k] = cvstatus(S[k],dh,Q[k]);

      if (status != S[k])
      {
         change = TRUE;
      }
   }
   return(change);
}                        /* End of linkstatus */


Solver::LinkStatus  Solver::cvstatus(Solver::LinkStatus s, double dh, double q)
/*
**--------------------------------------------------
**  Input:   s  = current status
**           dh = headloss
**           q  = flow
**  Output:  returns new link status                 
**  Purpose: updates status of a check valve.        
**--------------------------------------------------
*/
{
   /* Prevent reverse flow through CVs */
   if (abs(dh) > Htol)
   {
      if (dh < -Htol)     return(CLOSED);
      else if (q < -Qtol) return(CLOSED);
      else                return(FULL_OPEN);
   }
   else
   {
      if (q < -Qtol) return(CLOSED);
      else           return(s);
   }
}



int  Solver::valvestatus()
/*
**-----------------------------------------------------------------
**  Input:   none                                                
**  Output:  returns 1 if any pressure or flow control valve                   //(2.00.11 - LR)
**           changes status, 0 otherwise                                       //(2.00.11 - LR) 
**  Purpose: updates status for PRVs & PSVs whose status                       //(2.00.12 - LR)
**           is not fixed to OPEN/CLOSED
**-----------------------------------------------------------------
*/
{
   int   change = FALSE,            /* Status change flag      */
         i,k,                       /* Valve & link indexes    */
         n1,n2;                     /* Start & end nodes       */
      LinkStatus s;                         /* Valve status settings   */
   double hset;                     /* Valve head setting      */

   for (i=1; i<=_net->Nvalves; i++)                   /* Examine each valve   */
   {
      k = _net->Valve[i].Link;                        /* Link index of valve  */
	  if (S[k] == DISABLED) continue;            /* Valve status fixed   */
      n1 = _net->Link[k].N1;                          /* Start & end nodes    */
      n2 = _net->Link[k].N2;
      s  = S[k];                                /* Save current status  */

      switch (_net->Link[k].Type)                     /* Evaluate new status: */
      {
	  case Network::PRV:  hset = _net->Node[n2].El + K[k];
                    S[k] = prvstatus(k,s,hset,H[n1],H[n2]);
                    break;
         case Network::PSV:  hset = _net->Node[n1].El + K[k];
                    S[k] = psvstatus(k,s,hset,H[n1],H[n2]);
                    break;

         default:   continue;
      }
      /* Check for status change */
      if (s != S[k])
      {
         change = TRUE;
      }
   }
   return(change);
}                       /* End of valvestatus() */

Solver::LinkStatus  Solver::prvstatus(int k, LinkStatus  s, double hset, double h1, double h2)
/*
**-----------------------------------------------------------
**  Input:   k    = link index                                
**           s    = current status                            
**           hset = valve head setting                        
**           h1   = head at upstream node                     
**           h2   = head at downstream node                   
**  Output:  returns new valve status                         
**  Purpose: updates status of a pressure reducing valve.     
**-----------------------------------------------------------
*/
{
   LinkStatus   status;     /* New valve status */
   double hml;        /* Minor headloss   */
   double htol = Htol;

   status = s;
   if (S[k] == DISABLED) return(status);       /* Status fixed by user */
   hml = _net->Link[k].Km*SQR(Q[k]);                /* Head loss when open  */

/*** Status rules below have changed. ***/                                     //(2.00.11 - LR)

   switch (s)
   {
      case PARTIAL:
         if (Q[k] < -Qtol)            status = CLOSED;
         else if (h1-hml < hset-htol) status = FULL_OPEN;                           //(2.00.11 - LR)
         else                         status = PARTIAL;
         break;
      case FULL_OPEN:
         if (Q[k] < -Qtol)            status = CLOSED;
         else if (h2 >= hset+htol)    status = PARTIAL;                         //(2.00.11 - LR)
         else                         status = FULL_OPEN;
         break;
      case CLOSED:
         if ( h1 >= hset+htol                                                  //(2.00.11 - LR)
           && h2 < hset-htol)         status = PARTIAL;                         //(2.00.11 - LR)
         else if (h1 < hset-htol                                               //(2.00.11 - LR)
               && h1 > h2+htol)       status = FULL_OPEN;                           //(2.00.11 - LR)
         else                         status = CLOSED;
         break;

   }
   return(status);
}


Solver::LinkStatus  Solver::psvstatus(int k, Solver::LinkStatus s, double hset, double h1, double h2)
/*
**-----------------------------------------------------------
**  Input:   k    = link index                                
**           s    = current status                            
**           hset = valve head setting                        
**           h1   = head at upstream node                     
**           h2   = head at downstream node                   
**  Output:  returns new valve status                         
**  Purpose: updates status of a pressure sustaining valve.   
**-----------------------------------------------------------
*/
{
   Solver::LinkStatus  status;       /* New valve status */
   double hml;          /* Minor headloss   */
   double htol = Htol;

   status = s;
   if (S[k] == DISABLED) return(status);       /* Status fixed by user */
   hml = _net->Link[k].Km*SQR(Q[k]);                /* Head loss when open  */

/*** Status rules below have changed. ***/                                     //(2.00.11 - LR)

   switch (s)
   {
      case PARTIAL:
         if (Q[k] < -Qtol)            status = CLOSED;
         else if (h2+hml > hset+htol) status = FULL_OPEN;                           //(2.00.11 - LR)
         else                         status = PARTIAL;
         break;
      case FULL_OPEN:
         if (Q[k] < -Qtol)            status = CLOSED;
         else if (h1 < hset-htol)     status = PARTIAL;                         //(2.00.11 - LR)
         else                         status = FULL_OPEN;
         break;
      case CLOSED:
         if (h2 > hset+htol                                                    //(2.00.11 - LR)
          && h1 > h2+htol)            status = FULL_OPEN;                           //(2.00.11 - LR)
         else if (h1 >= hset+htol                                              //(2.00.11 - LR)
               && h1 > h2+htol)       status = PARTIAL;                         //(2.00.11 - LR)
         else                         status = CLOSED;
         break;

   }
   return(status);
}



double Solver::newflows()
/*
**----------------------------------------------------------------
**  Input:   none                                                
**  Output:  returns solution convergence error                  
**  Purpose: updates link flows after new nodal heads computed   
**----------------------------------------------------------------
*/
{
   double  dh,                    /* Link head loss       */
           dq;                    /* Link flow change     */
   double  dqsum,                 /* Network flow change  */
           qsum;                  /* Network total flow   */
   int   k, n, n1, n2;

   /* Initialize net inflows (i.e., demands) at tanks */
   for (n=_net->Njuncs+1; n <= _net->Nnodes; n++) D[n] = 0.0;

   /* Initialize sum of flows & corrections */
   qsum  = 0.0;
   dqsum = 0.0;

   /* Update flows in all links */
   for (k=1; k<=_net->Nlinks; k++)
   {

      /*
      ** Apply flow update formula:                   
      **   dq = Y - P*(new head loss)                 
      **    P = 1/(dh/dq)                             
      **    Y = P*(head loss based on current flow)   
      ** where P & Y were computed in newcoeffs().   
      */

      n1 = _net->Link[k].N1;
      n2 = _net->Link[k].N2;
      dh = H[n1] - H[n2];
      dq = Y[k] - P[k]*dh;

      Q[k] -= dq;

      /* Update sum of absolute flows & flow corrections */
      qsum += abs(Q[k]);
      dqsum += abs(dq);

      /* Update net flows to tanks */
      if ( S[k] != DISABLED )                                                     //(2.00.12 - LR)
      {
         if (n1 > _net->Njuncs) D[n1] -= Q[k];
         if (n2 > _net->Njuncs) D[n2] += Q[k];
      }

   }

   /* Update emitter flows */
   for (k=1; k<=_net->Njuncs; k++)
   {
      if (_net->Node[k].Ke == 0.0 || 
		  (unsigned(_tabCAT[k]) & unsigned(Channel::B) ) ) // b-emitter-conflict
		  continue;
      dq = emitflowchange(k);
      E[k] -= dq;
      qsum += abs(E[k]);
      dqsum += abs(dq);
   }

   /* Return ratio of total flow corrections to total flow */
   if (qsum > Hacc) return(dqsum/qsum);
   else return(dqsum);

}                        /* End of newflows */


double  Solver::emitflowchange(int i)
/*
**--------------------------------------------------------------
**   Input:   i = node index
**   Output:  returns change in flow at an emitter node                                                
**   Purpose: computes flow change at an emitter node
**--------------------------------------------------------------
*/
{
   double ke, p;
   ke = max(CSMALL, _net->Node[i].Ke);
   p = _net->Qexp*ke*pow(abs(E[i]),(_net->Qexp-1.0));
   if (p < RQtol)
      p = 1/RQtol;
   else
      p = 1.0/p;
   return(E[i]/_net->Qexp - p*(H[i] - _net->Node[i].El));
}


int  Solver::badvalve(int n) {
/*  originally in hydraulic.c
**-----------------------------------------------------------------
**  Input:   n = node index                                                
**  Output:  returns 1 if node n belongs to an active control valve,
**           0 otherwise  
**  Purpose: determines if a node belongs to an active control valve
**           whose setting causes an inconsistent set of eqns. If so,
**           the valve status is fixed open and a warning condition
**           is generated.
**-----------------------------------------------------------------
*/
   int i,k,n1,n2;
   for (i=1; i<=_net->Nvalves; i++) {
      k = _net->Valve[i].Link;
      n1 = _net->Link[k].N1;
      n2 = _net->Link[k].N2;
      if (n == n1 || n == n2) {
         if (_net->Link[k].Type == Network::PRV ||
             _net->Link[k].Type == Network::PSV ) {
                if (S[k] == PARTIAL || S[k] == CLOSED) {
                    S[k] = FULL_OPEN; // force open
                   return(1);
                }
         }
         return(0);
      }
   }
   return(0);
}
   
void   Solver::newcoeffs()
	/*
	**--------------------------------------------------------------
	**  Input:   none                                                
	**  Output:  none                                                
	**  Purpose: computes coefficients of linearized network eqns.   
	**--------------------------------------------------------------
	*/
{
	memset(Aii,0,(N+1)*sizeof(double));   /* Reset coeffs. to 0 */
	memset(Aij,0,(_net->Ncoeffs+1)*sizeof(double));
	memset(F,0,(N+1)*sizeof(double));
	memset(X,0,(N+1)*sizeof(double));
	memset(P,0,(M+1)*sizeof(double));
	memset(Y,0,(M+1)*sizeof(double));

	linkcoeffs();                            /* Compute link coeffs.  */
	emittercoeffs();                         /* Compute emitter coeffs.*/
	nodecoeffs();                            /* Compute node coeffs.  */
	valvecoeffs();                           /* Compute valve coeffs. */
}                        /* End of newcoeffs */


void  Solver::nodecoeffs()
	/*
	**----------------------------------------------------------------
	**  Input:   none                                                
	**  Output:  none                                                
	**  Purpose: completes calculation of nodal flow imbalance (X)   
	**           & flow correction (F) arrays                        
	**----------------------------------------------------------------
	*/
{
	int   i;

	/* For junction nodes, subtract demand flow from net */
	/* flow imbalance & add imbalance to RHS array F.    */
	for (i=1; i<=_net->Njuncs; i++)
	{
		X[i] -= D[i];
		F[_net->Row[i]] += X[i];
	}
}                        /* End of nodecoeffs */


void  Solver::valvecoeffs()
	/*
	**--------------------------------------------------------------
	**   Input:   none                                                
	**   Output:  none                                                
	**   Purpose: computes matrix coeffs. for PRVs, PSVs & FCVs       
	**            whose status is not fixed to OPEN/CLOSED            
	**--------------------------------------------------------------
	*/
{
	int i,k,n1,n2;

	for (i=1; i<=_net->Nvalves; i++)                   /* Examine each valve   */
	{
		k = _net->Valve[i].Link;                        /* Link index of valve  */
		if (S[k] == DISABLED ) continue;            /* Valve status fixed   */
		n1 = _net->Link[k].N1;                          /* Start & end nodes    */
		n2 = _net->Link[k].N2; 
		switch (_net->Link[k].Type)                     /* Call valve-specific  */
		{                                         /*   function           */
		case Network::PRV:  prvcoeff(k,n1,n2); break;
		case Network::PSV:  psvcoeff(k,n1,n2); break;
			//         case FCV:  fcvcoeff(k,n1,n2); break;
		default:   continue;
		}
	}
}                        /* End of valvecoeffs */

void  Solver::prvcoeff(int k, int n1, int n2)
	/*
	**--------------------------------------------------------------
	**   Input:   k    = link index                                   
	**            n1   = upstream node of valve                       
	**            n2   = downstream node of valve                       
	**   Output:  none                                                
	**   Purpose: computes solution matrix coeffs. for pressure       
	**            reducing valves                                     
	**--------------------------------------------------------------
	*/
{
	int   i,j;                       /* _net->Rows of solution matrix */
	double hset;                      /* Valve head setting      */
	i  = _net->Row[n1];                    /* Matrix rows of nodes    */
	j  = _net->Row[n2];
	hset   = _net->Node[n2].El + K[k];     /* Valve setting           */

	if (S[k] != DISABLED)
	{
		/*
		Set coeffs. to force head at downstream 
		node equal to valve setting & force flow (when updated in       
		newflows()) equal to flow imbalance at downstream node. 
		*/
		P[k] = 0.0;
		Y[k] = Q[k] + X[n2];       /* Force flow balance   */
		F[j] += (hset*CBIG);       /* Force head = hset    */
		Aii[j] += CBIG;            /*   at downstream node */
		if (X[n2] < 0.0) F[i] += X[n2];
		return;
	}

	/* 
	For OPEN, CLOSED, or XPRESSURE valve
	compute matrix coeffs. using the valvecoeff() function.                  //(2.00.11 - LR)
	*/
	valvecoeff(k);  /*pipecoeff(k);*/                                           //(2.00.11 - LR)
	Aij[_net->Ndx[k]] -= P[k];
	Aii[i] += P[k];
	Aii[j] += P[k];
	F[i] += (Y[k]-Q[k]);
	F[j] -= (Y[k]-Q[k]);
}                        /* End of prvcoeff */


void  Solver::psvcoeff(int k, int n1, int n2)
	/*
	**--------------------------------------------------------------
	**   Input:   k    = link index                                   
	**            n1   = upstream node of valve                       
	**            n2   = downstream node of valve                       
	**   Output:  none                                                
	**   Purpose: computes solution matrix coeffs. for pressure       
	**            sustaining valve                                    
	**--------------------------------------------------------------
	*/
{
	int   i,j;                       /* _net->Rows of solution matrix */
	double hset;                      /* Valve head setting      */
	i  = _net->Row[n1];                    /* Matrix rows of nodes    */
	j  = _net->Row[n2];
	hset   = _net->Node[n1].El + K[k];     /* Valve setting           */

	if (S[k] != DISABLED)
	{
		/*
		Set coeffs. to force head at upstream 
		node equal to valve setting & force flow (when updated in       
		newflows()) equal to flow imbalance at upstream node. 
		*/
		P[k] = 0.0;
		Y[k] = Q[k] - X[n1];              /* Force flow balance   */
		F[i] += (hset*CBIG);              /* Force head = hset    */
		Aii[i] += CBIG;                   /*   at upstream node   */
		if (X[n1] > 0.0) F[j] += X[n1];
		return;
	}

	/* 
	For OPEN, CLOSED, or XPRESSURE valve
	compute matrix coeffs. using the valvecoeff() function.                  //(2.00.11 - LR)
	*/
	valvecoeff(k);  /*pipecoeff(k);*/                                           //(2.00.11 - LR)
	Aij[_net->Ndx[k]] -= P[k];
	Aii[i] += P[k];
	Aii[j] += P[k];
	F[i] += (Y[k]-Q[k]);
	F[j] -= (Y[k]-Q[k]);
}                        /* End of psvcoeff */



void  Solver::linkcoeffs()
	/*
	**--------------------------------------------------------------
	**   Input:   none                                                
	**   Output:  none                                                
	**   Purpose: computes solution matrix coefficients for links     
	**--------------------------------------------------------------
	*/
{
	int   k,n1,n2;

	/* Examine each link of network */
	for (k=1; k<=M; k++)
	{
		n1 = _net->Link[k].N1;           /* Start node of link */
		n2 = _net->Link[k].N2;           /* End node of link   */

		/* Compute P[k] = 1 / (dh/dQ) and Y[k] = h * P[k]   */
		/* for each link k (where h = link head loss).      */
		/* FCVs, PRVs, and PSVs with non-fixed status       */
		/* are analyzed later.                              */

		switch (_net->Link[k].Type)
		{
		case Network::CV:
		case Network::PIPE:  pipecoeff(k); break;
		case Network::PUMP:  
		case Network::GPV:  pumpcoeffB(k); break;
		case Network::TCV:   tcvcoeff(k);  break;
//		case Network::PRV:
//		case Network::PSV:   
//			/* If valve status fixed then treat as pipe */
//			/* otherwise ignore the valve for now. */
//			if (S[k] != Network::ACTIVE) // valve pressure setting not specified.
//				valvecoeff(k);  //pipecoeff(k);      //(2.00.11 - LR)    
//			else continue;
//			break;
		default:    continue;                  
		}                                         

		/* Update net nodal inflows (X), solution matrix (A) and RHS array (F) */
		/* (Use covention that flow out of node is (-), flow into node is (+)) */
		X[n1] -= Q[k];
		X[n2] += Q[k];
		Aij[_net->Ndx[k]] -= P[k];              /* Off-diagonal coeff. */
		if (n1 <= _net->Njuncs)                 /* Node n1 is junction */
		{
			Aii[_net->Row[n1]] += P[k];          /* Diagonal coeff. */
			F[_net->Row[n1]] += Y[k];            /* RHS coeff.      */
		}
		else F[_net->Row[n2]] += (P[k]*H[n1]);  /* Node n1 is a tank   */
		if (n2 <= _net->Njuncs)                 /* Node n2 is junction */
		{
			Aii[_net->Row[n2]] += P[k];          /* Diagonal coeff. */
			F[_net->Row[n2]] -= Y[k];            /* RHS coeff.      */
		}
		else  F[_net->Row[n1]] += (P[k]*H[n2]); /* Node n2 is a tank   */
	}
}                        /* End of linkcoeffs */

void Solver::pumpcoeffB(int k) {
	// cut pump link 
	P[k] = 1.0/CBIG;
	Y[k] = Q[k];   

	//and add a reservoir to the discharge side for [B] channels
	// similar to emitter 
//	double  ke=CSMALL;
//	double  p;
//	double  q;
//	double  y;
//	double  z;

//    int n2 = _net->Link[k].N2;
//	q = E[n2];  // E gets updated in the previous iteration
//	z = ke*pow(abs(q),_net->Qexp);
//	p = _net->Qexp*z/abs(q);
//	if (p < RQtol) p = 1.0/RQtol;
//	else p = 1.0/p;
//	y = SGN(q)*z*p;

//	Aii[_net->Row[n2]] += p;
//    // set the head of the fictitious reservoir to be B channel data;
//	F[_net->Row[n2]] += y + p*B[n2]; 
//	X[n2] -= q;

}


void  Solver::pipecoeff(int k)
	/*
	**--------------------------------------------------------------
	**   Input:   k = link index                                      
	**   Output:  none                                                
	**  Purpose:  computes P & Y coefficients for pipe k              
	**                                                              
	**    P = inverse head loss gradient = 1/(dh/dQ)                
	**    Y = flow correction term = h*P                            
	**--------------------------------------------------------------
	*/
{
	double  hpipe,     /* Normal head loss          */
		hml,       /* Minor head loss           */
		ml,        /* Minor loss coeff.         */
		p,         /* q*(dh/dq)                 */
		q,         /* Abs. value of flow        */
		r,         /* Resistance coeff.         */
		r1,        /* Total resistance factor   */
		f,         /* D-W friction factor       */
		dfdq;      /* Derivative of fric. fact. */

	/* For closed pipe use headloss formula: h = CBIG*q */
	if (S[k] == DISABLED)
	{
		P[k] = 1.0/CBIG;
		Y[k] = Q[k];
		return;
	}

	/* Evaluate headloss coefficients */
	q = abs(Q[k]);                         /* Absolute flow       */
	ml = _net->Link[k].Km;                       /* Minor loss coeff.   */
	r = _net->Link[k].R;                         /* Resistance coeff.   */
	f = 1.0;                               /* D-W friction factor */
	if (_net->Formflag == Network::DW) f = DWcoeff(k,&dfdq);   
	r1 = f*r+ml;

	/* Use large P coefficient for small flow resistance product */
	if (r1*q < RQtol)
	{
		P[k] = 1.0/RQtol;
		Y[k] = Q[k]/Hexp;
		return;
	}

	/* Compute P and Y coefficients */
	if (_net->Formflag == Network::DW)                  /* D-W eqn. */
	{
		hpipe = r1*pow(q,2);                /* Total head loss */
		p = 2.0*r1*q;                     /* |dh/dQ| */
		/* + dfdq*r*q*q;*/                 /* Ignore df/dQ term */
		p = 1.0/p;
		P[k] = p;
		Y[k] = SGN(Q[k])*hpipe*p;
	}
	else                                 /* H-W or C-M eqn.   */
	{
		hpipe = r*pow(q,Hexp);            /* Friction head loss  */
		p = Hexp*hpipe;                   /* Q*dh(friction)/dQ   */
		if (ml > 0.0)
		{
			hml = ml*q*q;                  /* Minor head loss   */
			p += 2.0*hml;                  /* Q*dh(Total)/dQ    */
		}
		else  hml = 0.0;
		p = Q[k]/p;                       /* 1 / (dh/dQ) */
		P[k] = abs(p);
		Y[k] = p*(hpipe + hml);
	}
}                        /* End of pipecoeff */
void  Solver::tcvcoeff(int k)
	/*
	**--------------------------------------------------------------
	**   Input:   k = link index                                      
	**   Output:  none                                                
	**   Purpose: computes P & Y coeffs. for throttle control valve   
	**--------------------------------------------------------------
	*/
{
	double km, tpkm, p;

	km = K[k];

	/* If valve not fixed OPEN or CLOSED, compute its loss coeff. */
	//   if (K[k] != Network::_TINY)
	if (S[k] != DISABLED) {// active
		tpkm = 0.02517*K[k]/(SQR(km)*SQR(km));
		p = 2.0*km*fabs(Q[k]);
		if ( p < RQtol ) p = RQtol;
		P[k] = 1.0/p;
		Y[k] = Q[k]/2.0;
	} else {
		P[k] = 1.0/CBIG;
		Y[k] = Q[k];
	}
}                        /* End of tcvcoeff */

void Solver::valvecoeff(int k)
	/*
	**--------------------------------------------------------------
	**   Input:   k    = link index                                   
	**   Output:  none                                                
	**   Purpose: computes solution matrix coeffs. for a completely
	**            open, or closed valve.                                               
	**--------------------------------------------------------------
	*/
{
	double p;

	// Valve is closed. Use a very small matrix coeff.
	if (S[k] == DISABLED)
	{
		P[k] = 1.0/CBIG;
		Y[k] = Q[k];
		return;
	}

	// Account for any minor headloss through the valve
	if (_net->Link[k].Km > 0.0)
	{
		p = 2.0*_net->Link[k].Km*fabs(Q[k]);
		if ( p < RQtol ) p = RQtol;
		P[k] = 1.0/p;
		Y[k] = Q[k]/2.0;
	}
	else // valve is open
	{
		P[k] = 1.0/RQtol;
		Y[k] = Q[k];
	}
}
void  Solver::emittercoeffs()
	/*
	**--------------------------------------------------------------
	**   Input:   none                                                
	**   Output:  none                                                
	**   Purpose: computes matrix coeffs. for emitters
	**
	**   Note: Emitters consist of a fictitious pipe connected to
	**         a fictitious reservoir whose elevation equals that
	**         of the junction. The headloss through this pipe is
	**         Ke*(Flow)^Qexp, where Ke = emitter headloss coeff.
	**--------------------------------------------------------------
	*/
{
	int   i;
	double  ke;
	double  p;
	double  q;
	double  y;
	double  z;
	for (i=1; i<=_net->Njuncs; i++)
	{
		if (_net->Node[i].Ke == 0.0 ||  // not an emitter
			(unsigned(_tabCAT[i]) & unsigned(Channel::B))  )  // has a B channel here
			continue;
		ke = max(CSMALL, _net->Node[i].Ke);
		q = E[i];
		z = ke*pow(abs(q),_net->Qexp); // nodal head loss 
		p = _net->Qexp*z/abs(q);
		if (p < RQtol) p = 1.0/RQtol;
		else p = 1.0/p;
		y = SGN(q)*z*p;
		Aii[_net->Row[i]] += p;
		F[_net->Row[i]] += y + p*_net->Node[i].El; 
		X[i] -= q;
	}
}


double Solver::DWcoeff(int k, double *dfdq)
	/*
	**--------------------------------------------------------------
	**   Input:   k = link index                                      
	**   Output:  returns Darcy-Weisbach friction factor              
	**   Purpose: computes Darcy-Weisbach friction factor             
	**                                                              
	**    Uses interpolating polynomials developed by               
	**    E. Dunlop for transition flow from 2000 < Re < 4000.      
	**
	**   df/dq term is ignored as it slows convergence rate.
	**--------------------------------------------------------------
	*/
{
	double q,             /* Abs. value of flow */
		f;             /* Friction factor    */
	double x1,x2,x3,x4,
		y1,y2,y3,
		fa,fb,r;
	double s,w;

	*dfdq = 0.0;
	if (_net->Link[k].Type > Network::PIPE) return(1.0); /* Only apply to pipes */
	q = abs(Q[k]);
	s = Viscos*_net->Link[k].Diam;
	w = q/s;                       /* w = Re(Pi/4) */
	if (w >= A1)                   /* Re >= 4000; Colebrook Formula */
	{
		y1 = A8/pow(w,0.9);
		y2 = _net->Link[k].Kc/(3.7*_net->Link[k].Diam) + y1;
		y3 = A9*log(y2);
		f = 1.0/pow(y3,2);
		/*  *dfdq = (2.0+AA*y1/(y2*y3))*f; */   /* df/dq */
	}
	else if (w > A2)              /* Re > 2000; Interpolation formula */
	{
		y2 = _net->Link[k].Kc/(3.7*_net->Link[k].Diam) + AB;
		y3 = A9*log(y2);
		fa = 1.0/pow(y3,2);
		fb = (2.0+AC/(y2*y3))*fa;
		r = w/A2;
		x1 = 7.0*fa - fb;
		x2 = 0.128 - 17.0*fa + 2.5*fb;
		x3 = -0.128 + 13.0*fa - (fb+fb);
		x4 = r*(0.032 - 3.0*fa + 0.5*fb);
		f = x1 + r*(x2 + r*(x3+x4));
		/*  *dfdq = (x1 + x1 + r*(3.0*x2 + r*(4.0*x3 + 5.0*x4)));  */
	}
	else if (w > A4)              /* Laminar flow: Hagen-Poiseuille Formula */
	{
		f = A3*s/q;
		/*  *dfdq = A3*s; */
	}
	else
	{
		f = 8.0;
		*dfdq = 0.0;
	}
	return(f);
}                        /* End of DWcoeff */

//
//void  Solver::nodecoeffs()
//	/*
//	**----------------------------------------------------------------
//	**  Input:   none                                                
//	**  Output:  none                                                
//	**  Purpose: completes calculation of nodal flow imbalance (X)   
//	**           & flow correction (F) arrays                        
//	**----------------------------------------------------------------
//	*/
//{
//	int   i;

//	/* For junction nodes, subtract demand flow from net */
//	/* flow imbalance & add imbalance to RHS array F.    */
//	for (i=1; i<=_net->Njuncs; i++)
//	{
//		X[i] -= D[i];
//		F[_net->Row[i]] += X[i];
//	}
//}                        /* End of nodecoeffs */

int  Solver::linsolve(int n, double *Aii, double *Aij, double *B)
/*
**--------------------------------------------------------------
** Input:   n    = number of equations                          
**          Aii  = diagonal entries of solution matrix          
**          Aij  = non-zero off-diagonal entries of matrix      
**          B    = right hand side coeffs.                      
** Output:  B    = solution values                              
**          returns 0 if solution found, or index of            
**          equation causing system to be ill-conditioned       
** Purpose: solves sparse symmetric system of linear            
**          equations using Cholesky factorization              
**                                                              
** NOTE:   This procedure assumes that the solution matrix has  
**         been symbolically factorized with the positions of   
**         the lower triangular, off-diagonal, non-zero coeffs. 
**         stored in the following integer arrays:              
**            _net->XLNZ  (start position of each column in _net->NZSUB)    
**            _net->NZSUB (row index of each non-zero in each column) 
**            _net->LNZ   (position of each _net->NZSUB entry in Aij array) 
**                                                              
**  This procedure has been adapted from subroutines GSFCT and  
**  GSSLV in the book "Computer Solution of Large Sparse        
**  Positive Definite Systems" by A. George and J. W-H Liu      
**  (Prentice-Hall, 1981).                                      
**--------------------------------------------------------------
*/
{
   int    *link, *first;
   int    i, istop, istrt, isub, j, k, kfirst, newk;
   int    errcode = 0;
   double bj, diagj, ljk;
   double *temp;

   temp = (double *) calloc(n+1, sizeof(double));
   link = (int *) calloc(n+1,sizeof(int));
   first = (int *) calloc(n+1,sizeof(int));
   //ERRCODE(MEMCHECK(temp));
   //ERRCODE(MEMCHECK(link));
   //ERRCODE(MEMCHECK(first));
   if (errcode)
   {
      errcode = -errcode;
      goto ENDLINSOLVE;
   }
   memset(temp,0,(n+1)*sizeof(double));
   memset(link,0,(n+1)*sizeof(int));

   /* Begin numerical factorization of matrix A into L */
   /*   Compute column L(*,j) for j = 1,...n */
   for (j=1; j<=n; j++)
   {
      /* For each column L(*,k) that affects L(*,j): */
      diagj = 0.0;
      newk = link[j];
      k = newk;
      while (k != 0)
      {

         /* Outer product modification of L(*,j) by  */
         /* L(*,k) starting at first[k] of L(*,k).   */
         newk = link[k];
         kfirst = first[k];
         ljk = Aij[_net->LNZ[kfirst]];
         diagj += ljk*ljk;
         istrt = kfirst + 1;
         istop = _net->XLNZ[k+1] - 1;
         if (istop >= istrt)
         {

	     /* Before modification, update vectors 'first' */
	     /* and 'link' for future modification steps.   */
            first[k] = istrt;
            isub = _net->NZSUB[istrt];
            link[k] = link[isub];
            link[isub] = k;

	    /* The actual mod is saved in vector 'temp'. */
            for (i=istrt; i<=istop; i++)
            {
               isub = _net->NZSUB[i];
               temp[isub] += Aij[_net->LNZ[i]]*ljk;
            }
         }
         k = newk;
      }

      /* Apply the modifications accumulated */
      /* in 'temp' to column L(*,j).         */
      diagj = Aii[j] - diagj;
      if (diagj <= 0.0)        /* Check for ill-conditioning */
      {
         errcode = j;
         goto ENDLINSOLVE;
      }
      diagj = sqrt(diagj);
      Aii[j] = diagj;
      istrt = _net->XLNZ[j];
      istop = _net->XLNZ[j+1] - 1;
      if (istop >= istrt)
      {
         first[j] = istrt;
         isub = _net->NZSUB[istrt];
         link[j] = link[isub];
         link[isub] = j;
         for (i=istrt; i<=istop; i++)
         {
            isub = _net->NZSUB[i];
            bj = (Aij[_net->LNZ[i]] - temp[isub])/diagj;
            Aij[_net->LNZ[i]] = bj;
            temp[isub] = 0.0;
         }
      }
   }      /* next j */

   /* Foward substitution */
   for (j=1; j<=n; j++)
   {
      bj = B[j]/Aii[j];
      B[j] = bj;
      istrt = _net->XLNZ[j];
      istop = _net->XLNZ[j+1] - 1;
      if (istop >= istrt)
      {
         for (i=istrt; i<=istop; i++)
         {
            isub = _net->NZSUB[i];
            B[isub] -= Aij[_net->LNZ[i]]*bj;
         }
      }
   }

   /* Backward substitution */
   for (j=n; j>=1; j--)
   {
      bj = B[j];
      istrt = _net->XLNZ[j];
      istop = _net->XLNZ[j+1] - 1;
      if (istop >= istrt)
      {
         for (i=istrt; i<=istop; i++)
         {
            isub = _net->NZSUB[i];
            bj -= Aij[_net->LNZ[i]]*B[isub];
         }
      }
      B[j] = bj/Aii[j];
   }

ENDLINSOLVE:
   free(temp);
   free(link);
   free(first);
   return(errcode);
}                        /* End of linsolve */
