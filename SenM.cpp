#include <process.h>
#include <atltime.h>
#include "DataSource.h"
#include "SenM.h"

vtkStandardNewMacro(SenMInterationStyle);

SenMInterationStyle::Err SenMInterationStyle::init(SenM* handler_in, CRITICAL_SECTION* propProtector) { 
	_debug = 0;
	if (handler_in == NULL) return NO_EVENT_HANDLER;
	_pk = vtkPropPicker::New();
	_pkCell = vtkCellPicker::New();
	_hl = handler_in;
	_propProtector = propProtector;
	UseTimersOn();

	return OK;
}

void SenMInterationStyle::OnMouseMove() {
    vtkRenderWindow* rw = this->GetInteractor()->GetRenderWindow();
    int x,y;
    this->GetInteractor()->GetEventPosition(x, y);
    //_pk->PickProp(x, y, rw->GetRenderers()->GetFirstRenderer());
    //vtkPropAssembly* ac = _pk->GetPropAssembly();

    rw->GetRenderers()->InitTraversal();
    vtkRenderer* left = rw->GetRenderers()->GetNextItem();
    vtkRenderer* right = rw->GetRenderers()->GetNextItem();

    if (left && left->IsInViewport(x,y)) {
        if (_debug) _tprintf(TEXT("Int-style (%p), Try picking from left (%p) at (%d,%d) using picker (%p).\n"), this, left, x, y, _pkCell);
        //_pk->Pick(x, y, 0, left);
        if (_pkCell->Pick(x, y, 0, left)) {
            //if (_debug) _tprintf(TEXT("Done picking.\n"));

            vtkPropAssembly* pa = _pkCell->GetPropAssembly();
            if (_debug) _tprintf(TEXT("Done getting propassembly %p.\n"), pa);
            
            //if (pa != NULL) {
            //vtkIdType cid = _pk->GetCellId();
            
            rw->SetCurrentCursor(VTK_CURSOR_HAND);
            if (pa) _hl->updateSelection(this, pa);
            
            OnTimer();
        } else {
            rw->SetCurrentCursor(VTK_CURSOR_ARROW);
        }
    }

    if (right && right->IsInViewport(x,y)) {
        _pkCell->Pick(x, y, 0, right);
        vtkMapper* mp=vtkMapper::SafeDownCast(_pkCell->GetMapper());
        vtkIdType cid = _pkCell->GetCellId();

        if (mp != NULL){
            rw->SetCurrentCursor(VTK_CURSOR_HAND);
            _hl->updateSelection(this, mp, cid);
            //rw->Render();
            OnTimer();
        } else {
        rw->SetCurrentCursor(VTK_CURSOR_ARROW);
        }
    }
    
    if (_propProtector!=NULL) EnterCriticalSection(_propProtector);
    vtkInteractorStyleTrackballCamera::OnMouseMove();
    if (_propProtector!=NULL) LeaveCriticalSection(_propProtector);
}

void SenMInterationStyle::OnTimer() {
	if (_propProtector!=NULL) EnterCriticalSection(_propProtector);
	this->GetInteractor()->GetRenderWindow()->Render();
	vtkInteractorStyleTrackballCamera::OnTimer();
	if (_propProtector!=NULL) LeaveCriticalSection(_propProtector);
}

void SenMInterationStyle::OnLeftButtonDown() {
	if (_propProtector!=NULL) EnterCriticalSection(_propProtector);
	vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
	if (_propProtector!=NULL) LeaveCriticalSection(_propProtector);
}

void SenMInterationStyle::OnLeftButtonUp() {
	if (_propProtector!=NULL) EnterCriticalSection(_propProtector);
	vtkInteractorStyleTrackballCamera::OnLeftButtonUp();
	if (_propProtector!=NULL) LeaveCriticalSection(_propProtector);
}

void SenMInterationStyle::OnMiddleButtonDown() {
	if (_propProtector!=NULL) EnterCriticalSection(_propProtector);
	vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
	if (_propProtector!=NULL) LeaveCriticalSection(_propProtector);
}

void SenMInterationStyle::OnMiddleButtonUp() {
	if (_propProtector!=NULL) EnterCriticalSection(_propProtector);
	vtkInteractorStyleTrackballCamera::OnMiddleButtonUp();
	if (_propProtector!=NULL) LeaveCriticalSection(_propProtector);
}

void SenMInterationStyle::OnRightButtonDown() {
	if (_propProtector!=NULL) EnterCriticalSection(_propProtector);
	vtkInteractorStyleTrackballCamera::OnRightButtonDown();
	if (_propProtector!=NULL) LeaveCriticalSection(_propProtector);
}

void SenMInterationStyle::OnRightButtonUp() {
	if (_propProtector!=NULL) EnterCriticalSection(_propProtector);
	vtkInteractorStyleTrackballCamera::OnRightButtonUp();
	if (_propProtector!=NULL) LeaveCriticalSection(_propProtector);
}

void SenMInterationStyle::OnMouseWheelForward() {
	if (_propProtector!=NULL) EnterCriticalSection(_propProtector);
	vtkInteractorStyleTrackballCamera::OnMouseWheelForward();
	if (_propProtector!=NULL) LeaveCriticalSection(_propProtector);
}

void SenMInterationStyle::OnMouseWheelBackward() {
	if (_propProtector!=NULL) EnterCriticalSection(_propProtector);
	vtkInteractorStyleTrackballCamera::OnMouseWheelBackward();
	if (_propProtector!=NULL) LeaveCriticalSection(_propProtector);
}


SenM::Err SenM::init(LPCTSTR netn, LPCTSTR dcs, CTime* st, unsigned dt, Varima* dm) {
	//show  time 
	CTime cur_time = CTime::GetCurrentTime();
	_tprintf(TEXT("[0] Init SenM... Current system time: %s, current local time: %s\n"), 
		cur_time.FormatGmt(TEXT("Year %Y, Day %j at %H:%M:%S")), 
		cur_time.Format(TEXT("%c")) );
	
	// fetch time
	_senm_time = *st;


	//init database and network file
	_tprintf(TEXT("[0] Loading network (.inp) file...\n"));
	if (Network::getNetwork(netn, &_net)) {
		_ftprintf(stderr, TEXT("Cannot open network file %s  \n"), netn);
		return NET_FILE_ERROR;
	}
	_tprintf(_net->report());

	if (dcs == NULL) {
		_hasDB = 0; //no database, (testing only)
		_tprintf(TEXT("[0] No Database...\n"));
	} else {
		_tprintf(TEXT("[0] Connecting Scada Database...\n"));
		_tprintf(TEXT("\tAccess string: %s\n"), dcs);
		Provider* db_conn = new Provider;
		if (db_conn->connect(dcs)) {
			_ftprintf(stderr, TEXT("Cannot connect DB %s \n"), dcs);
			return DB_CONN_ERROR;
		}

		_tprintf(TEXT("[0] Loading channels...\n"));
		if (dt == 0) {
			_ftprintf(stderr, TEXT("Time quantum must be greater than zero. \n"));
		}
		_source = new DataSource(dt);
		if (_source->addProvider(db_conn)) {
			_ftprintf(stderr, TEXT("Database data Error \n"));
			return DB_DATA_ERROR;
		}
		_dt = dt;
		_tprintf(TEXT("\tTime quantum: %d, Number of Channels: %d\n"), _dt, _source->getNChan());
	}

	_tprintf(TEXT("[0] Building base network...\n"));
	_net2d = vtkPolyData::New();
	_net3d = vtkPolyData::New();
	
	_net->get2d3dNet(_net2d, _net3d);
	// compute data bounds for better vis/ exeggration
	_net3d->ComputeBounds();
	_net3d->GetBounds(bounds);
	vtkTransformFilter* etff = vtkTransformFilter::New();
	vtkTransform* etf = vtkTransform::New();
	_z_scaling = (bounds[1]-bounds[0]+bounds[3]-bounds[2])/2*_elev_exf/(bounds[5]-bounds[4]);
	_xy_unit = (bounds[1]-bounds[0]+bounds[3]-bounds[2])/400;
	etf->Scale(1,1, _z_scaling); //elevation scale
	etff->SetTransform(etf);
	etff->SetInput(_net3d);

	_pmNet2d = vtkPolyDataMapper::New();
	_pmNet2d->SetInput(_net2d);
	_pmNet3d = vtkPolyDataMapper::New();
	_net3dex = vtkPolyData::SafeDownCast(etff->GetOutput());
	_net3dex->Update();
	_net3dex->ComputeBounds();
	_pmNet3d->SetInput(_net3dex);
	
	//check demand model
	_tprintf(TEXT("[0] Checking demand model...\n"));
	if (dm == NULL)  {

		_ftprintf(stderr, TEXT("ARIMA model Error \n"));
		return ARIMA_ERROR;
	}
	_dmodel = dm;

	Err err = OK;
	err = initPara();
	//err = initOpenGL();
	
	_ftprintf(stderr, TEXT("About to init scada data\n"));

	if (_hasDB) err = initScada();

	return err;
}
//
//SenM::Err SenM::initOpenGL() {
//
//	vtkOpenGLExtensionManager *extensions = vtkOpenGLExtensionManager::New();
//	extensions->SetRenderWindow(_vwT);
//	_vwT->Render();
//
//
//	const char *gl_vendor =
//		reinterpret_cast<const char *>(glGetString(GL_VENDOR));
//	const char *gl_version =
//		reinterpret_cast<const char *>(glGetString(GL_VERSION));
//	const char *gl_renderer =
//		reinterpret_cast<const char *>(glGetString(GL_RENDERER));
//
//	cout << endl;
//	cout << "GL_VENDOR: " << (gl_vendor ? gl_vendor : "(null)") << endl;
//	cout << "GL_VERSION: " << (gl_version ? gl_version : "(null)") << endl;
//	cout << "GL_RENDERER: " << (gl_renderer ? gl_renderer : "(null)") << endl;
//
//	cout << endl;
//	_vwT->Print(cout);
//
//	cout << "LoadSupportedExtension..." << endl;
//	int supported=extensions->ExtensionSupported("GL_VERSION_1_2");
//	int loaded=0;
//	if(supported)
//	{
//		cout << "Driver claims to support OpenGL 1.2" <<endl;
//		loaded=extensions->LoadSupportedExtension("GL_VERSION_1_2");
//		if(loaded)
//		{
//			cout << "OpenGL 1.2 features loaded." <<endl;
//		}
//		else
//		{
//			cout << "Failed to load OpenGL 1.2 features!" <<endl;
//		}
//	}
//	supported=extensions->ExtensionSupported("GL_VERSION_1_3");
//	if(supported)
//	{
//		cout << "Driver claims to support OpenGL 1.3" <<endl;
//		loaded=extensions->LoadSupportedExtension("GL_VERSION_1_3");
//		if(loaded)
//		{
//			cout << "OpenGL 1.3 features loaded." <<endl;
//		}
//		else
//		{
//			cout << "Failed to load OpenGL 1.3 features!" <<endl;
//		}
//	}
//	supported=extensions->ExtensionSupported("GL_VERSION_1_4");
//	if(supported)
//	{
//		cout << "Driver claims to support OpenGL 1.4" <<endl;
//		loaded=extensions->LoadSupportedExtension("GL_VERSION_1_4");
//		if(loaded)
//		{
//			cout << "OpenGL 1.4 features loaded." <<endl;
//		}
//		else
//		{
//			cout << "Failed to load OpenGL 1.4 features!" <<endl;
//		}
//	}
//	supported=extensions->ExtensionSupported("GL_VERSION_1_5");
//	if(supported)
//	{
//		cout << "Driver claims to support OpenGL 1.5" <<endl;
//		loaded=extensions->LoadSupportedExtension("GL_VERSION_1_5");
//		if(loaded)
//		{
//			cout << "OpenGL 1.5 features loaded." <<endl;
//		}
//		else
//		{
//			cout << "Failed to load OpenGL 1.5 features!" <<endl;
//		}
//	}
//	supported=extensions->ExtensionSupported("GL_VERSION_2_0");
//	if(supported)
//	{
//		cout << "Driver claims to support OpenGL 2.0" <<endl;
//		loaded=extensions->LoadSupportedExtension("GL_VERSION_2_0");
//		if(loaded)
//		{
//			cout << "OpenGL 2.0 features loaded." <<endl;
//		}
//		else
//		{
//			cout << "Failed to load OpenGL 2.0 features!" <<endl;
//		}
//	}
//	supported=extensions->ExtensionSupported("GL_VERSION_2_1");
//	if(supported)
//	{
//		cout << "Driver claims to support OpenGL 2.1" <<endl;
//		loaded=extensions->LoadSupportedExtension("GL_VERSION_2_1");
//		if(loaded)
//		{
//			cout << "OpenGL 2.1 features loaded." <<endl;
//		}
//		else
//		{
//			cout << "Failed to load OpenGL 2.1 features!" <<endl;
//		}
//	}
//	cout << "GetExtensionsString..." << endl;
//	cout << extensions->GetExtensionsString() << endl;
//	return OK;
//}



SenM::Err SenM::initPara() {

	if (_dmodel->getN() != _net->getNusers()) {
		_ftprintf(stderr, TEXT(
			"ARIMA model dimension doesn't match number of network water users.\n"));
		return ARIMA_ERROR;
	}
	_n = _dmodel->getN();
	_tprintf(TEXT("\tNetwork water users: %d\n"), _n);

	// vis pipeline for parameter Phi
	_tprintf(TEXT("[0] Building AR parameter visualization...\n"));
	vtkFloatArray* users = vtkFloatArray::New();
	users->DeepCopy(vtkFloatArray::SafeDownCast(
		_net2d->GetPointData()->GetScalars()));
	int nPhi = _dmodel->getNphi();
	int nTheta = _dmodel->getNtheta(); 
	_nPara = nPhi + nTheta;
	if (_nPara > MAX_ZLINE_ARRAY_SIZE) _nPara = MAX_ZLINE_ARRAY_SIZE; //limit the no. of parameters shown
	
	int i; //parameter iterator
	vtkIdType j; //node iterator
	vtkIdType k; //user iterator
	int Nnode = users->GetNumberOfTuples();
	_vrPhi = vtkRenderer::New();
	_acPhi = new vtkActor*[_nPara * _n];
	_aPhi = new vtkPropAssembly*[_n];
	_dsPhi = vtkCubeSource::New();
	_pmPhi = vtkPolyDataMapper::New();
	_pmPhi->SetInput(_dsPhi->GetOutput());
	_pmPhi->ScalarVisibilityOff();

	double extent_ref = (bounds[1]-bounds[0]+bounds[3]-bounds[2])/2;
	
	double a = extent_ref /200; //x,y length of the column
	double dth = _net->getBaseDemandPercentile(1-_pcPara);
	
	for (j=0, k=0; j<Nnode; ++j) {
		
		if (users->GetValue(j) > 0) {//water user
			_aPhi[k] = vtkPropAssembly::New();
			for (i=0; i<_nPara; ++i) {
				_acPhi[i*_n + k] = vtkActor::New();
				_acPhi[i*_n + k]->SetMapper(_pmPhi);

				float h = 0;
				if (i+1>nPhi)  //MA param
					h = (float)(_dmodel->getTheta(k, i+1-nPhi));
				else  //AR
					h = (float)(_dmodel->getPhi(k, i+1));

				double coords[3];
				_net2d->GetPoint(j, coords);

				//build the parameter bars
				
				double tph = (bounds[1]-bounds[0]+bounds[3]-bounds[2])/2*_para_exf * h;
				double zmin, zmax;
				if (tph>0) {
					zmin = 0; zmax = tph;
				} else {
					zmin = tph; zmax = 0;
				}
				_acPhi[i*_n + k]->SetScale(a, a, tph);
				_acPhi[i*_n + k]->SetPosition(
					coords[0]+a*zlineArrayOffset[_nPara][i][0],
					coords[1]+a*zlineArrayOffset[_nPara][i][1],
					tph/2);

				_acPhi[i*_n + k]->GetProperty()->SetColor(
					specColor[i][0]/255,specColor[i][1]/255,specColor[i][2]/255);
				_acPhi[i*_n + k]->GetProperty()->SetInterpolationToFlat();
				_aPhi[k]->AddPart(_acPhi[i*_n + k]);
			}
			_vrPhi->AddActor(_aPhi[k]);
			if (users->GetValue(j) < dth) _aPhi[k]->VisibilityOff(); //hide small users
			_hashPhi[_aPhi[k]] = k;
			++k;
		}
		
	}


	_sHLBox = NULL; 
	_aNet2d = vtkActor::New();
	//_aNet3d = vtkActor::New();
	_pmNet2d->ScalarVisibilityOff();
	_aNet2d->SetMapper(_pmNet2d);
	_aNet2d->GetProperty()->SetColor(0,1,1);
	//_aNet3d->SetMapper(_pmNet3d);

	_dsBase = vtkPlaneSource::New();
	_dsBase->SetOrigin(bounds[0]-extent_ref/100, bounds[2]-extent_ref/100, 0);
	_dsBase->SetPoint1(bounds[0]-extent_ref/100, bounds[3]+extent_ref/100, 0);
	_dsBase->SetPoint2(bounds[1]+extent_ref/100, bounds[2]-extent_ref/100, 0);
	_pmBase = vtkPolyDataMapper::New();
	_pmBase->SetInput(_dsBase->GetOutput());
	_aBase = vtkActor::New();
	_aBase->SetMapper(_pmBase);
	_aBase->GetProperty()->SetOpacity(0.3);
	_aBase->PickableOff();
	_vrPhi->AddActor(_aBase);

	vtkVertexGlyphFilter* vgf = vtkVertexGlyphFilter::New();
	vgf->SetInput(_net2d);
	_pmNodes = vtkPolyDataMapper::New();
	_pmNodes->SetInput(vgf->GetOutput());
	_pmNodes->ScalarVisibilityOff();
	_acNodes = vtkActor::New();
	_acNodes->SetMapper(_pmNodes);
	_acNodes->GetProperty()->SetPointSize(6);
	_acNodes->GetProperty()->SetColor(0,1,0);

	_vrPhi->AddActor(_aNet2d);

	// use a green box to highlight current selection
	_sHLBox = vtkCubeSource::New();
	_pmHLBox = vtkPolyDataMapper::New();
	_acHLBox = vtkActor::New();   
	_pmHLBox->SetInput(_sHLBox->GetOutput());
	_acHLBox->SetMapper(_pmHLBox);
	_acHLBox->GetProperty()->SetRepresentationToWireframe();
	_acHLBox->GetProperty()->SetColor(0, 1, 0);
	_acHLBox->VisibilityOff();
	_vrPhi->AddActor(_acHLBox);

	//tool tip text 
	//alternative tooltip
	_acPhiT2 = vtkCaptionActor2D::New();
	//_tPhi = vtkTextSource::New();
	//_tPhi->BackingOn();
	//_tPhi->SetBackgroundColor(0.2, 0.2, 0.2);
	//_tPhi->SetForegroundColor(1,1,1);
	
	_acPhiT2->GetCaptionTextProperty()->SetFontSize(40);
	_acPhiT2->GetCaptionTextProperty()->SetFontFamilyToArial();
	_acPhiT2->GetCaptionTextProperty()->ItalicOff();
	_acPhiT2->GetCaptionTextProperty()->BoldOff();
	_acPhiT2->GetCaptionTextProperty()->ShadowOff();
	
	_acPhiT2->VisibilityOff();
	_acPhiT2->LeaderOn();
	_acPhiT2->BorderOff();
	_vrPhi->AddActor2D(_acPhiT2);
	
	// build the covariance vis
	_tprintf(TEXT("[0] Building covariance matrix Visualization...\n"));
	_vrCov = vtkRenderer::New();
	int ii, jj; //iterator of covariance
	int kkV = 0; //cell counter
	int kkC = 0;

	_dsVar = new vtkPolyData*[_n];
	_pmVar = new vtkPolyDataMapper*[_n];
	_acVar = new vtkActor*[_n];

	_dsCov = new vtkPolyData*[_n];
	_pmCov = new vtkPolyDataMapper*[_n];
	_acCov = new vtkActor*[_n];
	double covRange = _dmodel->getMaxCov() - _dmodel->getMinCov();
	for (ii=0; ii<_n; ++ii) {
		double tpVar = _dmodel->getCov(ii, ii);
		_dsCov[ii] = vtkPolyData::New();
		_pmCov[ii] = vtkPolyDataMapper::New();
		_acCov[ii] = vtkActor::New();
		_dsVar[ii] = vtkPolyData::New();
		_pmVar[ii] = vtkPolyDataMapper::New();
		_acVar[ii] = vtkActor::New();
		vtkPoints* ptsC = vtkPoints::New(); //covariance
		ptsC->SetNumberOfPoints(_n*_n);
		vtkCellArray* linesC = vtkCellArray::New();
		linesC->SetNumberOfCells(_n);
		vtkFloatArray* attC = vtkFloatArray::New();
		attC->SetNumberOfTuples(_n);
		vtkPoints* ptsV = vtkPoints::New(); //variance
		ptsC->SetNumberOfPoints(_n*2);
		vtkCellArray* linesV = vtkCellArray::New();
		linesC->SetNumberOfCells(1);
		vtkFloatArray* attV = vtkFloatArray::New();
		attV->SetNumberOfTuples(1);


		kkC = 0; kkV = 0;
		double coords[3];
		_net2d->GetPoint(_net->userId2NodeId(ii), coords); // get x,y

		for (jj=0; jj<=ii; ++jj) {
			double tpCov = _dmodel->getCov(ii, jj);
			coords[2] =  tpCov / covRange * extent_ref*_para_exf*5;
			vtkIdType p[2];
			if (ii==jj) { // vertical lines showing variance
				
					p[0] = _n*_n + ii; //point number for line end 1, 
					p[1] = _n*_n + ii + _n; //line end 2
					ptsV->InsertPoint(p[1], coords[0], coords[1], -coords[2]);
					ptsV->InsertPoint(p[0], coords[0], coords[1], coords[2]);
					linesV->InsertNextCell(2, p);
					attV->InsertTuple1(kkV++, coords[2]);
				
			} else { // horizontal lines showing covariance
				if (abs(tpCov)/tpVar >= _covTh  //ignore weak dependencies
					&& users->GetValue(jj) >= dth){  //ignore small users
					double coords2[3];
					_net2d->GetPoint(_net->userId2NodeId(jj), coords2);
					coords2[2] = _dmodel->getCov(jj, jj) / covRange * extent_ref*_para_exf*5;

					p[0] = ii*_n + jj;
					p[1] = jj*_n + ii;
					//if (users->GetValue(_net->userId2NodeId(ii)) >= dth) {//hide the small users

					ptsC->InsertPoint(p[0], coords[0], coords[1], coords[2]);
					ptsC->InsertPoint(p[1], coords2[0], coords2[1], coords[2]);
					linesC->InsertNextCell(2, p);
					attC->InsertTuple1(kkC++, coords[2]);


					ptsV->InsertPoint(p[0], coords2[0], coords2[1], coords2[2]);
					ptsV->InsertPoint(p[1], coords2[0], coords2[1], -coords2[2]);
					linesV->InsertNextCell(2, p);
					attV->InsertTuple1(kkV++, coords2[2]);

				}
				//}
			}
		}
		_dsCov[ii]->SetPoints(ptsC); ptsC->Delete();
		_dsCov[ii]->SetLines(linesC); linesC->Delete();
		_dsCov[ii]->GetCellData()->SetScalars(attC); attC->Delete();

		_dsVar[ii]->SetPoints(ptsV); ptsV->Delete();
		_dsVar[ii]->SetLines(linesV); linesV->Delete();
		_dsVar[ii]->GetCellData()->SetScalars(attV); attV->Delete();

		vtkTubeFilter *tub1 = vtkTubeFilter::New();
		vtkTubeFilter *tub2 = vtkTubeFilter::New();
		double tube_r = (bounds[1]-bounds[0]+bounds[3]-bounds[2])/1000;
		tub1->SetNumberOfSides(4);
		tub1->SetRadius(tube_r);
		tub2->SetNumberOfSides(4);
		tub2->SetRadius(tube_r);
		tub1->SetInput(_dsCov[ii]);
		tub2->SetInput(_dsVar[ii]);
		vtkTriangleFilter* trft = vtkTriangleFilter::New();
		trft->SetInput(tub1->GetOutput());
		vtkStripper* stp = vtkStripper::New();
		vtkStripper* stp2 = vtkStripper::New();
		stp->SetInput(trft->GetOutput());
		stp2->SetInput(tub2->GetOutput());
		

		_pmCov[ii]->SetInput(stp->GetOutput());
		_acCov[ii]->SetMapper(_pmCov[ii]);
		_acCov[ii]->VisibilityOff();
		_acCov[ii]->GetProperty()->SetColor(0,0,1);
		_pmVar[ii]->SetInput(stp2->GetOutput());
		_pmVar[ii]->ScalarVisibilityOff();
		//_pmVar[ii]->Update();
		//_pmVar[ii]->StaticOn();
		_acVar[ii]->SetMapper(_pmVar[ii]);
		_acVar[ii]->VisibilityOff();
		_acVar[ii]->GetProperty()->SetColor(1,1,1);

		_vrCov->AddActor(_acCov[ii]);
		_vrCov->AddActor(_acVar[ii]);

	}

	// use a green box to highlight selection
	//_sHLCovBox = vtkCubeSource::New();
	_sHLCovBox = vtkPolyData::New();
	_pmHLCovBox = vtkPolyDataMapper::New();
	_acHLCovBox = vtkActor::New();   
	_pmHLCovBox->SetInput(_sHLCovBox);
	_acHLCovBox->SetMapper(_pmHLCovBox);
	_acHLCovBox->GetProperty()->SetRepresentationToWireframe();
	_acHLCovBox->GetProperty()->SetColor(1, 1, 0);
	_acHLCovBox->VisibilityOff();
	_vrCov->AddActor(_acHLCovBox);


	_vrCov->AddActor(_aNet2d);
	_vrCov->AddActor(_aBase);
	//_vrCov->AddActor(_acNodes);
	return OK;

}

SenM::Err SenM::initScada() {
	//initialize scada data view
	_tprintf(TEXT("[0] Building Scada channel visuals...\n"));
	_vrScada = vtkRenderer::New();
		
	_nSensors = _source->getNChan(); 
	Channel* clist = _source->getChanList();
	_source->dumpChannelsInfo();

	Sensor::setCommonParameters(_net3dex, _z_scaling, _xy_unit); 
	// allocate mem, build sensor geometry
	_sensors = new Sensor*[_nSensors];
	_aSensors = new vtkPropAssembly*[_nSensors];
	for (int i=0; i<_nSensors; ++i) {
		_sensors[i] = Sensor::makeSensor(
			_net->netId2CellId(clist->name, clist->mtype), clist);
		clist = clist->next;
		if (_sensors[i] == NULL) {
			_tprintf(TEXT("\tChannel No. %d can't be shown...\n"), i);
			_aSensors[i] = NULL;
			continue;
		}
		_aSensors[i] = _sensors[i]->getGeom();
		_hashScada[_aSensors[i]] = _sensors[i];
		_vrScada->AddActor(_aSensors[i]);
	}

	_tmTimeBanner = vtkTextMapper::New();
	_tmTimeBanner->SetInput("Network");
	_tmTimeBanner->GetTextProperty()->SetFontSize(15);

	_acTimeBanner = vtkActor2D::New();
	_acTimeBanner->SetMapper(_tmTimeBanner);
	_acTimeBanner->GetProperty()->SetColor(1, 1, 0);

	return OK;
	
}



unsigned WINAPI SenM::_startParaEl(void* obj) {//parameter view
	SenM* s = static_cast<SenM*>(obj);

	_tprintf(TEXT("[1] Building parameter visualization window...\n"));
	s->_vwT = vtkRenderWindow::New();
	s->_vwT->SetMultiSamples(0);
	s->_vwT->SetSize(800,400);
	
	s->_vwT->AddRenderer(s->_vrPhi); 
	s->_vwT->AddRenderer(s->_vrCov);

	vtkLight *aLight = s->_vrPhi->MakeLight();
	aLight->SetLightTypeToHeadlight();
	s->_vrPhi->AddLight(aLight);
	s->_vrPhi->SetLightFollowCamera(1);
	s->_vrCov->AddLight(aLight);
	s->_vrCov->SetLightFollowCamera(1);

	s->_vrPhi->SetViewport(0, 0, 0.5, 1);
	s->_vrPhi->SetBackground(0.05,0.08,0.05);
	s->_vrCov->SetViewport(0.5, 0, 1, 1);
	s->_vrPhi->SetBackground(0.08,0.05,0);
	s->_camP = vtkCamera::New();
	s->_camP->SetClippingRange(10, 1000);
	//_cam->SetParallelProjection(1);
	s->_vrPhi->SetActiveCamera(s->_camP);
	s->_vrCov->SetActiveCamera(s->_camP);
	s->_vrPhi->ResetCamera();
	s->_vrCov->ResetCamera();
	
	s->_viT = vtkRenderWindowInteractor::New();
	s->_sisPhi = SenMInterationStyle::New();
	
	s->_viT->SetInteractorStyle(s->_sisPhi);
	s->_viT->SetRenderWindow(s->_vwT);

	s->_viT->CreateRepeatingTimer(50);
	s->_viT->Initialize();
	s->_sisPhi->init(s, NULL); // no protection for multi-threading
	
	//s->initOpenGL();
	s->_vwT->Render();
	s->_vwT->SetWindowName("VARIMA Model Parameter View");
	
	_tprintf(TEXT("[1] Parameter viewer event loop started...\n"));
		
	s->_viT->Start();
	return 0;
}

unsigned WINAPI SenM::_startScadaEl(void* obj) {//scada view
	SenM* s = static_cast<SenM*>(obj);

		
	_tprintf(TEXT("[2] Building Scada visualization window...\n"));
	//initialize scada view;
	s->_vwS = vtkRenderWindow::New();
	s->_vwS->SetSize(500,500);

	s->_pmNetS = vtkPolyDataMapper::New();
		s->_pmNetS->SetInput(s->_net3dex);
		s->_pmNetS->ScalarVisibilityOff();
		
	s->_acNetS = vtkActor::New();
		s->_acNetS->SetMapper(s->_pmNetS);
		//s->_acNetS->GetProperty()->SetColor(30.0/255,144.0/255,1); //blue pipelines
		s->_acNetS->GetProperty()->SetColor(0, 1, 0);
		s->_acNetS->GetProperty()->SetPointSize(4);
	
	s->_vrScada->AddActor(s->_acNetS);

	//debug only, show cell ids
	if (s->_debug) {
		vtkIdFilter* idf = vtkIdFilter::New();
		idf->SetInput(s->_net3dex);
		idf->CellIdsOn();

		vtkCellCenters* cc = vtkCellCenters::New();
		cc->SetInput(idf->GetOutput());

		vtkLabeledDataMapper* ldm = vtkLabeledDataMapper::New();
		ldm->SetInput(cc->GetOutput());
		ldm->SetLabelModeToLabelFieldData();
		ldm->GetLabelTextProperty()->SetColor(1,1,1);

		vtkActor2D* celllabels = vtkActor2D::New();
		celllabels->SetMapper(ldm);

		s->_vrScada->AddActor2D(celllabels);
	}

	

	s->_vwS->AddRenderer(s->_vrScada);

	s->_viS = vtkRenderWindowInteractor::New();
	s->_sisS = SenMInterationStyle::New();
	
	s->_viS->SetInteractorStyle(s->_sisS);
	s->_vwS->SetInteractor(s->_viS);
	
	EnterCriticalSection(&(s->_access_sensor));
	s->_vwS->Render();
	LeaveCriticalSection(&(s->_access_sensor));
	s->_vwS->SetWindowName("SCADA View");

	s->_acTimeBanner->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
	s->_acTimeBanner->GetPositionCoordinate()->SetValue(0.1, 0.9);

	s->_vrScada->AddActor2D(s->_acTimeBanner);

	_tprintf(TEXT("[2] Scada viewer event loop started...\n"));
	s->_viS->CreateRepeatingTimer(50);
	s->_viS->Initialize();
	s->_sisS->init(s, &(s->_access_sensor));
	s->_viS->Start();

	return 0;
}

SenM::Err SenM::run() {
	// main thread

	//spawn threads for visualization
	_beginthreadex(NULL, 0, &SenM::_startParaEl, this, 0, NULL);
	
	if (_hasDB) {
		_beginthreadex(NULL, 0, &SenM::_startScadaEl, this, 0, NULL);

		if (_senm_time == NULL) 
			_senm_time = CTime::GetCurrentTime();

		_tprintf(TEXT("[0] SenM starts at Scada system time: %s, local time %s\n"), 
			_senm_time.FormatGmt(TEXT("Year %Y, Day %j at %H:%M:%S")), 
			_senm_time.Format(TEXT("%c")) );

		CTimeSpan dts(0, 0, 0, _dt);
		int nChan = _source->getNChan();
		double* snapshot = new double[nChan];

		for (;; _senm_time += dts) { // main scada loading cycle
			_source->getASnapshot(_senm_time, snapshot);
			for (int i=0; i<nChan; ++i) {
				if (_sensors[i] != NULL) {
					EnterCriticalSection(&_access_sensor);
					_sensors[i]->setData(snapshot[i]);
					LeaveCriticalSection(&_access_sensor);
				}
			}

			EnterCriticalSection(&_access_sensor);
			_tmTimeBanner->SetInput(CT2A(
				_senm_time.FormatGmt(TEXT("(%Y)Day %j (%H:%M:%S)"))));
			LeaveCriticalSection(&_access_sensor);

			//Sleep(_senm_delay * 1000);
		}

		delete[] snapshot;
	} else {
		Sleep(INFINITE);
	}

	
	_tprintf(TEXT("[0] Main thread is ending...\n"));
	return OK;
}

SenM::Err SenM::updateSelection(SenMInterationStyle* sis, vtkPropAssembly* ac) {
	if (sis == _sisPhi) { //called by phi window interator
		if (_debug) _tprintf(TEXT("update highlight box\n"));
		double* extents;
		if (extents = ac->GetBounds()) {//get a bound;
			_sHLBox->SetBounds(extents);
			_acHLBox->VisibilityOn();
		}

		if (_debug) _tprintf(TEXT("highlight box on. update banner\n"));
		int userid = _hashPhi[ac];
		int nodeid = _net->userId2NodeId(userid);
		char *nodename = NULL;
		char tpstr[512];
		_net->nodeId2NodeName(nodeid+1, &nodename);
		sprintf(_tttPhi, 
			"EPANET Node(%d): %s\nWater User ID: %d\nGraphics ID: %d\n",
			nodeid+1, nodename, userid, nodeid);
		for (int i=0; i<_nPara; ++i) {
			_dmodel->getParaDesc(i, tpstr);
			sprintf(tpstr, "%s%g\n", tpstr, _dmodel->getPara(userid, i));
			strcat(_tttPhi, tpstr);
		}
		sprintf(tpstr, "Variance: %g", _dmodel->getCov(userid, userid));
		strcat(_tttPhi, tpstr);
		if (_debug) _tprintf(TEXT("banner text done. show banner\n"));

		_acPhiT2->SetCaption(_tttPhi);
		_acPhiT2->SetAttachmentPoint(
			(extents[0]+extents[1])/2, (extents[2]+extents[3])/2, extents[5]);

		_acPhiT2->SetPosition(20,20);
		_acPhiT2->SetHeight(100);
		//_acPhiT2->SetPosition2(70,100);
		_acPhiT2->VisibilityOn();
		if (_debug) _tprintf(TEXT("banner shown\n"));
		// control covariance shown
		if (_lastSel != -1) {
			_acCov[_lastSel]->VisibilityOff();
			_acVar[_lastSel]->VisibilityOff();
		}
		_acCov[userid]->VisibilityOn();
		_acVar[userid]->VisibilityOn();
		_lastSel = userid;
		if (_debug) _tprintf(TEXT("update visibility\n"));

	} else if (sis == _sisS) {//scada viewer
		if (_lastSelSensor != NULL) _hashScada[_lastSelSensor]->unselect();
		_hashScada[ac]->select();
		_lastSelSensor = ac;
	}
	return OK;
}

SenM::Err SenM::updateSelection(SenMInterationStyle* sis, vtkMapper* mp, vtkIdType cellid) {
	//double extents[6];
	//mp->GetInput()->GetCellBounds(cellid, extents);
	//_sHLCovBox->SetBounds(extents);

	if (_lastSel == -1) return OK;

	vtkCellArray* tp = vtkCellArray::New();
	tp->InsertNextCell(mp->GetInput()->GetCell(cellid));
	_sHLCovBox->SetPoints(vtkPolyData::SafeDownCast(mp->GetInput())->GetPoints());
	_sHLCovBox->SetStrips(tp);
	_acHLCovBox->VisibilityOn();
	

	return OK;
}

const double SenM::zlineArrayOffset[MAX_ZLINE_ARRAY_SIZE+1][MAX_ZLINE_ARRAY_SIZE][2] = {
	// the unit of the following values are relative to tube radius
	// and are {x,y} offsets to the center of the node
	{{0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}}, //no v-bars
	{{0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}}, // 1 bar in the center
	{{-0.5,0}, {0.5,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}}, 
	{{-1,0}, {0,0}, {1,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}}, 
	{{-1.5,0}, {-0.5,0}, {0.5,0}, {1.5,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}}, 
	{{-1,0.5}, {0,0.5}, {1,0.5}, {-1,-0.5}, {0,-0.5}, {0,0}, {0,0}, {0,0}, {0,0}}, // 5 bars in 2 rows
	{{-1,0.5}, {0,0.5}, {1,0.5}, {-1,-0.5}, {0,-0.5}, {1,-0.5}, {0,0}, {0,0}, {0,0}},
	{{-1.5,0.5}, {-0.5,0.5}, {0.5,0.5}, {1.5,0.5}, {-1.5,-0.5}, {-0.5,-0.5}, {0.5,-0.5}, {0,0}, {0,0}},
	{{-1.5,0.5}, {-0.5,0.5}, {0.5,0.5}, {1.5,0.5}, {-1.5,-0.5}, {-0.5,-0.5}, {0.5,-0.5}, {1.5,-0.5}, {0,0}},
	{{-1,1}, {0,1}, {1,1}, {-1,0}, {0,0}, {1,0}, {-1,-1}, {0,-1}, {1,-1}}
};

double SenM::specColor[MAX_ZLINE_ARRAY_SIZE][3] = {
	///*Light Pink */	{255 ,	182 	,193 },//	#ffb6c1 
	///*Light Coral*/	{240 ,	128 ,	128}, 
	///*Light Salmon */	{255 ,	160 	,122 },//	#ffa07a 
	///*Light Gray */	{211, 	211 	,211 },//	#d3d3d3 	   
	///*Light Green */	{144 ,	238 ,	144 },//	#90ee90 
	//	 	   
	///*Light Cyan */	{224 ,	255 ,	255}, 	 	   
	///*Light Goldenrod*/ 	{238 ,	221 ,	130 },		   
	///*Light Goldenrod Yellow */	{250 ,	250 	,210 },//	#fafad2 	      
	///*Light Sea Green */	{32 ,	178 ,	170 } //	#20b2aa 	    


	/*Dark Magenta*/	{139.0, 	0, 	139.0},	
	/*Dark Orange*/ 	{255, 	140.0, 	0} 	,
	/*Dark Red */	{139.0, 0, 0},
	/*Dark Olive Green*/ 	{85.0,	107.0, 	47.0},
	/*Dark Orchid*/ 	{153.0, 	50.0, 	204.0},
	/*Dark Goldenrod*/{184.0, 134.0, 11.0},
	/*Dark Gray*/ 	{169.0,	169.0, 169.0 },		   
	/*Dark Green*/	{0 ,	100.0, 	0 },

	/*Dark Khaki*/ 	{189.0,	183.0, 	107.0}
};