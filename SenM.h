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

/* class SenM is the central class for running sensor-driven
water distribution system models*/

#include <Windows.h>
#include <tchar.h>
#include <stdio.h>
#include <stdlib.h>

#include "vtkincludes.h"
#include "Sensor.h"
#include "DataSource.h"
#include "Viewer.h"
#include "Population.h"
#include "Varima.h"

//Number of working threads, should = no. of processers in the computer
#define N_WORKING_THREADS 4
#define MAX_ZLINE_ARRAY_SIZE 9  //in _vfPhi view, at most 9 variables are shown

class SenM;

class SenMInterationStyle : public vtkInteractorStyleTrackballCamera {
	// iteration style for handling mouse events
protected:
	int _debug;
	vtkPropPicker* _pk;
	vtkCellPicker* _pkCell;
	SenM* _hl; // event handler
	CRITICAL_SECTION* _propProtector;

public:
	enum Err { OK, NO_RENDERER, NO_EVENT_HANDLER };
	
	static SenMInterationStyle* New();
	vtkTypeMacro(SenMInterationStyle, vtkInteractorStyleTrackballCamera);

	// set the mouse event handler (updateSelection) and critical section for
	// protecting rendering
	Err init(SenM* handler_in, CRITICAL_SECTION* propProtector) ;

	
	virtual void OnMouseMove() ;
	//have to redefine the following for multi-threading

	virtual void 	OnLeftButtonDown ();

	virtual void 	OnLeftButtonUp ();

	virtual void 	OnMiddleButtonDown ();

	virtual void 	OnMiddleButtonUp ();

	virtual void 	OnRightButtonDown ();

	virtual void 	OnRightButtonUp ();

	virtual void 	OnMouseWheelForward ();

	virtual void 	OnMouseWheelBackward ();


	virtual void OnTimer();
};


class SenM {

public:

enum Err {
	OK,
	NET_FILE_ERROR,
	DB_CONN_ERROR,
	DB_DATA_ERROR,
	ARIMA_ERROR
};

SenM(): _dbcstr(NULL),_netfile(NULL),
		_source(NULL), _net(NULL), _dmodel(NULL),
		_elev_exf(1.0/4), _para_exf(1.0/20)
{
	_pcPara = 0.9; //by default only 10% of the node phi - cov parameters are visualized
	_lastSel = _lastSelCov = -1;
	_lastSelSensor = NULL;
	_covTh = 0.3;
	_senm_delay = 1;
	_debug = 0;
	_hasDB = 1;
	InitializeCriticalSection(&_access_sensor);
}


public:
	Err init(LPCTSTR network_file_name,
			 LPCTSTR database_connection_string, 
			 CTime*	 senm_start_time, /* if NULL use current time*/
			 unsigned senm_time_quantum_in_seconds,
			 Varima* demand_model);

	Err updateSelection(SenMInterationStyle*, vtkPropAssembly*);
	Err updateSelection(SenMInterationStyle*, vtkMapper* selMp, vtkIdType cellid);
	Err run();

protected:
	//three event loops starting functions for parameter, scada, and demand viewing
	static unsigned WINAPI _startParaEl(void* obj); 
	static unsigned WINAPI _startScadaEl(void* obj);
	static unsigned WINAPI _startDemandEl(void* obj);
	
	//init subroutine
	//Err initOpenGL();
	Err initPara();
	Err initScada();

protected:
	int _debug; //debug flag

	LPCTSTR _dbcstr;  //db connection string
	int _hasDB;  //if using a database, if not, the program is only
	//for testing graphic performance, no estimation can be done.

	LPCTSTR _netfile;  //network file name

	DataSource* _source;  //db datasource
	Network* _net;  //the epanet network model
	Varima* _dmodel;  //demand model

	int _n; //number of water users, also dimension of the unkown water demands
			
	//current populations of water demands
	//Population* _pops[N_WORKING_THREADS];  

	CTime _senm_time;  //time of the modeling system
	unsigned _dt; // db time quantum (in sec)

	int _senm_delay; // senm data fetch time delay (for testing) (in sec)

	//viewer objects
	vtkRenderWindow *_vwT; //parameter viewer windows

	vtkPolyData *_net2d, *_net3d, *_net3dex; //2d and 3d base network dataset, 3d network w/ elevation exaggerated
	vtkPolyDataMapper *_pmNodes;
	vtkActor* _acNodes;
	double bounds[6];  //3d network space bounds
	vtkPolyDataMapper *_pmNet2d, *_pmNet3d;
	vtkActor *_aNet2d, *_aNet3d;  // base network actors
	vtkPlaneSource* _dsBase; //base plane
	vtkPolyDataMapper *_pmBase;
	vtkActor *_aBase;
	double _z_scaling; // z-axis exaggeration factor
	double _xy_unit; // x,y axis unit length

	vtkRenderer *_vrPhi, *_vrCov; // parameter (phi, cov) view
	double _pcPara; //parameter display percentiles
	double _covTh; //covariance display threshold
	int _nPara; //number of linear parameters
	vtkCubeSource *_dsPhi; //polydata for parameters, dim = _nPara * nnodes
	vtkPolyData **_dsCov, **_dsVar; //polydata of covariance
	vtkPolyDataMapper *_pmPhi, **_pmCov, **_pmVar; //data mapper array for phi and data mapper of covariance
	vtkActor **_acPhi, **_acCov, **_acVar;
	vtkPropAssembly **_aPhi; // all bars for a node makes a group
	std::map<vtkPropAssembly*, vtkIdType> _hashPhi; 
	
	int _lastSel, _lastSelCov; //last selected user
	vtkCubeSource *_sHLBox; //highlight box for selection
	vtkPolyData  *_sHLCovBox;
	vtkPolyDataMapper *_pmHLBox, *_pmHLCovBox;
	vtkActor *_acHLBox;
	vtkActor *_acHLCovBox;

	char _tttPhi[1024]; //tooltip text for phi
	vtkTextSource* _tPhi; 
	
	vtkCaptionActor2D *_acPhiT2;

	vtkTextMapper *_tmPhi, *_tmCov;  //tooltip texts

	vtkRenderWindowInteractor *_viT;
	SenMInterationStyle *_sisPhi;
	vtkCamera* _camP; //single camera for para spatial views

	//const arrays for align zline array for _vrphi's columns
	static const double zlineArrayOffset[MAX_ZLINE_ARRAY_SIZE+1][MAX_ZLINE_ARRAY_SIZE][2];
	static double specColor[MAX_ZLINE_ARRAY_SIZE][3]; //distinct rgb color


	vtkRenderWindow *_vwS; //scada viewer windows
	vtkRenderer *_vrScada; // scada spatial view
	vtkPolyDataMapper* _pmNetS;
	vtkActor* _acNetS;
	Sensor **_sensors; 
	vtkPropAssembly **_aSensors; // array of sensors
	std::map<vtkPropAssembly*, Sensor*> _hashScada; //map propassembly to sensor
	CRITICAL_SECTION _access_sensor;
	vtkPropAssembly *_lastSelSensor; //sensor previously selected
	int _nSensors;
	vtkTextMapper* _tmTimeBanner;  // time banner
	vtkActor2D* _acTimeBanner;

	//vtkRenderer *_vrPlc;  // scada head/tank level/control temporal view
	//vtkRenderer *_vrQd;  //scada flowrate/demand temporal view

	vtkRenderWindowInteractor *_viS; 
	SenMInterationStyle *_sisS;

	vtkRenderWindow* _vwD;  //demand viewer windows
	vtkRenderer *_vrDs, *_vrDt; // demand temporal and spatial views
	vtkRenderWindowInteractor* _viD;

	



	double _elev_exf; //elevation exaggration factor (apprx. ratio of elevation range/x-y range)
	double _para_exf; //parameter exaggration factor (apprx. ratio of parameter range/x-y range)


	

};


	