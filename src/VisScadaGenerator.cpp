#include <float.h>
#include "VisScadaGenerator.h"

#ifdef _DEBUG
#define debug 1
#else
#define debug 0
#endif

VSG_ERR VSG_new(TCHAR* inpfile_path, char const* ini_path, int seed, Varima* md, int timespan, VisScadaGenerator** vsg_out) {
	VSG_ERR err;

	VisScadaGenerator* prec = (VisScadaGenerator*)calloc(1,sizeof(VisScadaGenerator));
	if (prec == NULL) {
		VSG_ewi(err = VSG_MEM_NOT_ALLOCED);
		goto END;
	}

	// a ScadaGenerator is used internally
	SG_ERR sgerr;
	ScadaGenerator* sg;
	sgerr = SG_new(inpfile_path, ini_path, seed, md, &sg);
	if (sgerr != SG_OK) {
		VSG_ewi(err = VSG_SG_ERR); goto END;}
	sg->extObj = prec;  //register visualization hook in the sg
	sg->extUpdate = &VSG_vupdate;
	prec->_sg = sg;


	// build vtk visuals
	// 1. Build base 3d and 2d network
	int iNet;
	for (iNet =0 ;iNet < VSG_N_NETPOLY; ++iNet) {
		prec->_baseNet[iNet] = vtkPolyData::New();
	}

	double coords[3];
	prec->maxx = prec->maxy =  prec->maxz = -DBL_MAX;
	prec->minx = prec->miny =  prec->minz = DBL_MAX;

	vtkPoints* pts3 = vtkPoints::New();
	vtkPoints* pts2 = vtkPoints::New();
	vtkCellArray* ca = vtkCellArray::New();  // cells for nodes
	vtkCellArray* cl = vtkCellArray::New();  // cells for links

	int iNode;
	for (iNode=1; iNode<=sg->nnode; ++iNode) {
		coords[0] = Node[iNode].x;
		if (coords[0] > prec->maxx ) prec->maxx = coords[0]; //update bands
		if (coords[0] < prec->minx ) prec->minx = coords[0];
		coords[1] = Node[iNode].y;
		if (coords[1] > prec->maxy ) prec->maxy = coords[1];
		if (coords[1] < prec->miny ) prec->miny = coords[1];
		coords[2] = Node[iNode].El;
		if (coords[2] > prec->maxz ) prec->maxz = coords[2];
		if (coords[2] < prec->minz ) prec->minz = coords[2];

		pts3->InsertPoint(iNode-1, coords);
		vtkIdType icell = iNode-1;
		ca->InsertNextCell(1, &icell);

		// flat (2d) map
		coords[2] = 0;
		pts2->InsertPoint(iNode-1, coords);

	}

	prec->_baseNet[0]->SetPoints(pts3);
	prec->_baseNet[1]->SetPoints(pts2); 
	prec->_baseNet[0]->SetVerts(ca);
	prec->_baseNet[1]->SetVerts(ca);

	int iLink;
	for (iLink=1; iLink<sg->nlink; ++iLink) {
		vtkIdType jcell[2];
		jcell[0] = Link[iLink].N1 - 1; // to vis index
		jcell[1] = Link[iLink].N2 - 1;
		cl->InsertNextCell(2, jcell);
	}

	prec->_baseNet[0]->SetLines(cl);
	prec->_baseNet[1]->SetLines(cl);

	// set up net skelecton 3d and 2d
	vtkPolyDataMapper* net3d = vtkPolyDataMapper::New();
	vtkPolyDataMapper* net2d = vtkPolyDataMapper::New();

	// set mapper to initial poly datasets
	net3d->SetInput(prec->_baseNet[0]);
	net2d->SetInput(prec->_baseNet[1]);

	vtkActor *actn3d = vtkActor::New();
	vtkActor *actn2d = vtkActor::New();
	actn3d->SetMapper(net3d);
	actn2d->SetMapper(net2d);


	pts3->Delete(); 
	pts2->Delete();
	ca->Delete();
	cl->Delete();


	//2. Build reference zline's
	prec->zlineP = vtkPolyData::New();
	prec->zlineD = vtkPolyData::New();
	vtkPoints* ptsP = vtkPoints::New();
	vtkPoints* ptsD = vtkPoints::New();
	ptsP->InsertPoint(0, 0, 0, 0);
	ptsD->InsertPoint(0, 0, 0, 0);
	ptsP->InsertPoint(1, 0, 0, 
		(prec->maxx-prec->minx+prec->maxy-prec->miny)/4/(prec->maxz-prec->minz)); //set scale factor
	ptsD->InsertPoint(1, 0, 0, 
		(prec->maxx-prec->minx+prec->maxy-prec->miny)/4/(sg->maxd-sg->mind));
	vtkCellArray* cla = vtkCellArray::New();
	vtkIdType jcell[] = {0, 1};
	cla->InsertNextCell(2, jcell);
	prec->zlineP->SetPoints(ptsP);
	prec->zlineD->SetPoints(ptsD);
	prec->zlineP->SetLines(cla);
	prec->zlineD->SetLines(cla);

	//3. set up p and d visuals
	prec->_mVisP = vtkPolyDataMapper::New(); 
	prec->_mVisD = vtkPolyDataMapper::New(); 
	prec->_mVisP->SetScalarVisibility(0);
	prec->_mVisD->SetScalarVisibility(0);
	vtkActor *pact = vtkActor::New();
	vtkActor *dact = vtkActor::New();
	pact->SetMapper(prec->_mVisP);
	pact->GetProperty()->SetColor(1, 1, 0);
	pact->GetProperty()->SetOpacity(0.5);
	dact->SetMapper(prec->_mVisD);
	dact->GetProperty()->SetColor(0.2, 0.5, 1);
	dact->GetProperty()->SetOpacity(1);

	// Allocate p, d visuals storage
	int ntstep = sg->ntstep;
	int dur = sg->dur;
	prec->_timespan = timespan;
	int size = ntstep * (timespan/dur + 1);
	prec->_netStore = (vtkPolyData**)calloc(size* 2, sizeof(vtkPolyData*));
	prec->_dispStep = (int*)calloc(size, sizeof(int));
	//each time step has two datasets - pressure and demand
	for (int i=0; i< size * 2; ++i) {
		prec->_netStore[i] = vtkPolyData::New();
	}
	prec->_nNetStore = size;

	//4. set up banner
	prec->_timeTxtM = vtkTextMapper::New();
	prec->_timeTxtM->SetInput("Network");
	prec->_timeTxtM->GetTextProperty()->SetFontSize(14);

	vtkActor2D* banAct = vtkActor2D::New();
	banAct->SetMapper(prec->_timeTxtM);
	banAct->GetProperty()->SetColor(1, 1, 0);
	banAct->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
	banAct->GetPositionCoordinate()->SetValue(0.1, 0.9);

	//5. set up camera, renders, render windows, interactor and interactor 
	// styles, add actors to scene 
	//there is only one camera (for both renderer)
	vtkCamera *cam = vtkCamera::New();

	vtkRenderer *renp = vtkRenderer::New(); //3d view
	renp->AddActor(pact);
	renp->AddActor(actn3d);
	renp->AddActor2D(banAct);

	vtkRenderer *rend = vtkRenderer::New(); //2d flat view
	rend->AddActor(dact);
	rend->AddActor(actn2d);

	prec->_renWin = vtkRenderWindow::New();
	prec->_renWin->AddRenderer(renp); 
	renp->SetViewport(0, 0.5, 1, 1);
	prec->_renWin->AddRenderer(rend); 
	rend->SetViewport(0, 0, 1, 0.5);

	prec->_iren = vtkRenderWindowInteractor::New();

	ScadaGenVisInteratorStyle* ctrl =  ScadaGenVisInteratorStyle::New();
	ctrl->setVSG(prec);

	prec->_iren->SetInteractorStyle(ctrl);
	prec->_iren->SetRenderWindow(prec->_renWin);

	renp->SetActiveCamera(cam);
	renp->ResetCamera();
	renp->SetBackground(0,0,0);
	rend->SetActiveCamera(cam);
	rend->SetBackground(0,0,0);

	prec->_renWin->SetSize(800,600);

	// interact with data

	prec->_iren->Initialize();
	prec->_iren->CreateRepeatingTimer(50);

	*vsg_out = prec;
	return VSG_OK;

END:
	free(prec);
	*vsg_out = NULL;
	return err;
}

void VSG_delete(VisScadaGenerator* vsg) {
    //work needed here
}

void VSG_moveto(VisScadaGenerator* vsg, int facestep) {
	if (facestep<0 || facestep>*vsg->_datastep) return; //not in range
	vsg->_mVisP->SetInput(vsg->_netStore[facestep*2]);
	vsg->_mVisP->Update();
	vsg->_mVisD->SetInput(vsg->_netStore[facestep*2+1]);
	int tp = vsg->_dispStep[facestep];
	sprintf(vsg->timeTxt, "Day %02d, %02d:%02d:%02d", 
		tp/(3600*24),
		tp%(3600*24)/3600,
		tp%3600/60,
		tp%60);
	vsg->_timeTxtM->SetInput(vsg->timeTxt);
}

void VSG_vupdate(void* vsg_in, ScadaGenerator* sg) {
	VisScadaGenerator* vsg = (VisScadaGenerator*)vsg_in;

	//build new set of visuals

	//vis, to run this part, _net_array must be allocated
	vtkDoubleArray* pressures = vtkDoubleArray::New();
	vtkDoubleArray* dmds = vtkDoubleArray::New();
	for (int ii=0; ii<sg->nnode; ++ii) {
		pressures->InsertTuple1(ii, sg->vpressures[ii]);
		dmds->InsertTuple1(ii, sg->vdemands[ii]);
	}

	vtkPolyData* tmpoly1 = vtkPolyData::New(); 
	vtkPolyData* tmpoly2 = vtkPolyData::New();
	tmpoly1->DeepCopy(vsg->_baseNet[0]);
	tmpoly1->GetPointData()->SetScalars(pressures);
	tmpoly2->DeepCopy(vsg->_baseNet[1]);
	tmpoly2->GetPointData()->SetScalars(dmds);

    //turn referece z-lines to zlines at junctions
	vtkGlyph3D* gly = vtkGlyph3D::New();
	vtkGlyph3D* gly2 = vtkGlyph3D::New();
	gly->SetInput(tmpoly1);
	gly2->SetInput(tmpoly2);
	gly->SetSource(vsg->zlineP);
	gly2->SetSource(vsg->zlineD);
	gly->SetScaleModeToScaleByScalar();
	gly2->SetScaleModeToScaleByScalar();

    //clena the polydata
    vtkCleanPolyData *cleanFilter1 = vtkCleanPolyData::New();
    vtkCleanPolyData *cleanFilter2 = vtkCleanPolyData::New();
	cleanFilter1->SetInput(gly->GetOutput());
	cleanFilter2->SetInput(gly2->GetOutput());

	vtkTubeFilter* tubP = vtkTubeFilter::New();
	vtkTubeFilter* tubD = vtkTubeFilter::New();
	tubP->SetInputConnection(cleanFilter1->GetOutputPort());
	tubP->SetNumberOfSides(10);
	tubP->SetRadius((vsg->maxx-vsg->minx)/400);
	tubD->SetInputConnection(cleanFilter2->GetOutputPort());
	tubD->SetNumberOfSides(10);
	tubD->SetRadius((vsg->maxx-vsg->minx)/100);
	tubP->SetCapping(1); 
	tubD->SetCapping(1); 

	tubP->Update();
	tubD->Update();

	// add to repository visuals and time stamps
	vsg->_netStore[(*vsg->_datastep+1)*2] = tubP->GetOutput();
	vsg->_netStore[(*vsg->_datastep+1)*2+1] = tubD->GetOutput();
	vsg->_dispStep[*vsg->_datastep+1] = sg->elapTime;

	//update complete
	InterlockedExchangeAdd(vsg->_datastep, 1);

	pressures->Delete();
	dmds->Delete();
}



void VSG_ewi(VSG_ERR err) {
	TCHAR errTxts[VSG_DUMMY_LAST][MAX_ERR_STRING_SIZE] = {
		TEXT("OK."),
		TEXT("Memory could not be allocated for the VisScadaGenerator."),
		TEXT("Could not create the ScadaGenerator."),
	};

	TCHAR errTxt[MAX_ERR_STRING_SIZE+MAX_ERR_PREFIX_SIZE];
	_stprintf(errTxt, TEXT("<GEN>VisScadaGenerator: %s"), errTxts[err]);

	ewi(errTxt);
}

void VSG_run(VisScadaGenerator* vsg, Tstamp* t_start) {
	SG_MakeTArg* arguments = (SG_MakeTArg*)calloc(1, sizeof(SG_MakeTArg));
	arguments->sg = vsg->_sg;
	arguments->t = t_start;
	arguments->timespan = vsg->_timespan;

    // start work thread
	_beginthreadex(NULL, 0, &SG_makeT, arguments, 0, NULL);

    vsg->_renWin->Render();
    // start visualizer event loop
    vsg->_iren->Start();
}


void ScadaGenVisInteratorStyle::setVSG(VisScadaGenerator* vsg_in) {
	if (vsg_in == NULL) return;
	vsg = vsg_in;
	vsg->_datastep = &datastep;
	//GlobalWarningDisplayOff(); //don't display warnings (vtkOutput)
	UseTimersOn();
	printf("VSG controls:\n [Right Arrow] Go the next calculated time step;"
		"\n [Left Arrow] Go the the previous calculated time step;"
		"\n [G or g] Catch up with the most recently cacluated time step.\n");
}

void ScadaGenVisInteratorStyle::OnKeyPress() 
{
	// Get the keypress
	vtkRenderWindowInteractor *rwi = this->Interactor;
	std::string kc = rwi->GetKeySym();

	if (kc == "Down" || kc == "Right") {
		if (desiredstep == -1) desiredstep=facestep; //enter step mode
		else desiredstep = facestep+1;
	}
	if (kc == "Left" || kc == "Up" ) {
		if (desiredstep == -1) desiredstep=facestep;
		else desiredstep = facestep>=1?facestep-1:0;
	}
	if (kc == "G" || kc == "g") {
		desiredstep = -1;  //enter sync mode
	}
	OnTimer();
	// Forward events
	vtkInteractorStyleTrackballCamera::OnKeyPress();
}

void ScadaGenVisInteratorStyle::OnTimer() {
	long oldface = facestep;
	long olddata = datastep;
	if (desiredstep == -1) 
		facestep = olddata;  //sync face with data
	else 
		facestep = desiredstep>olddata? olddata: desiredstep; 
	// sync face with desired, if not available, sync with data

	if (facestep != oldface) { //needs update
		VSG_moveto(vsg, facestep);
	}
	this->Interactor->GetRenderWindow()->Render();
    //if (debug) fprintf(stderr, "on timer render\n");
	vtkInteractorStyleTrackballCamera::OnTimer();
}

vtkStandardNewMacro(ScadaGenVisInteratorStyle);

