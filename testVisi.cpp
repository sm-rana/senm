// Use Vtk to visualize network demands and hydraulics
#include <process.h>
#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkRenderer.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkObjectFactory.h"
#include "vtkProperty.h"
#include <vtkTextMapper.h>
#include <vtkTextProperty.h>
#include <vtkActor2D.h>
#include <vtkProperty2D.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include "vtkRenderer.h"
#include "ScadaGenerator.h"


class ScadaGenVisInteratorStyle: public vtkInteractorStyleTrackballCamera
{
	
protected:
	ScadaGenVisInteratorStyle():
		datastep(-1),
		facestep(-1),
		desiredstep(-1),
		vtkInteractorStyleTrackballCamera(){};
	volatile long datastep;
	long facestep;
	long desiredstep;
	int* dispstep; // display time step in seconds
	vtkPolyData** _netinfo; //polydata for p and d
	vtkPolyDataMapper* _netvis[2]; //handle of data mapper
	vtkTextMapper* _timeTxtM; 
	char timeTxt[256];
	ScadaGenerator* _sg; 
	Varima* _md;
	Tstamp* _t_init;
	int _t_sim;

public:
	static ScadaGenVisInteratorStyle* New();
	vtkTypeMacro(ScadaGenVisInteratorStyle, vtkInteractorStyleTrackballCamera);

	void initPDdata(
		ScadaGenerator* sg, 
		Varima* md, 
		vtkPolyDataMapper* netvis[2],
		vtkTextMapper* tm,
		Tstamp* t_init,
		int t_sim) {

		if (sg == NULL) return;

		int ntstep = sg->getNTstep();
		int dur = sg->getDur();
		int size = ntstep * (t_sim/dur + 1);
		_netinfo = new vtkPolyData* [size * 2];
		dispstep = new int[size];
		//each time step has two datasets - pressure and demand
		for (int i=0; i< size * 2; ++i) {
			_netinfo[i] = vtkPolyData::New();
		}
		_netvis[0] = netvis[0];
		_netvis[0]->SetScalarVisibility(0);
		_netvis[1] = netvis[1];
		_netvis[1]->SetScalarVisibility(0);
		_timeTxtM = tm;
		_sg = sg;
		_md = md;
		_t_init = t_init;
		_t_sim = t_sim;

		GlobalWarningDisplayOff(); //don't display warnings (vtkOutput)
		UseTimersOn();
		_sg->makeTd(_t_init, _t_sim, _md, _netinfo, &datastep, dispstep);
	}

	virtual void OnEnter() {
	//when mouse enter the window area
		
		vtkInteractorStyleTrackballCamera::OnEnter();

	}

	virtual void OnKeyPress() 
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

	virtual void OnTimer() {
		long oldface = facestep;
		long olddata = datastep;
		if (desiredstep == -1) 
			facestep = olddata;  //sync face with data
		else 
			facestep = desiredstep>olddata? olddata: desiredstep; 
		// sync face with desired, if not available, sync with data

		if (facestep != oldface) {
			//needs update
			_netvis[0]->SetInput(_netinfo[facestep*2]);
			_netvis[1]->SetInput(_netinfo[facestep*2+1]);
			int tp = dispstep[facestep];
			sprintf(timeTxt, "Day %02d, %02d:%02d:%02d", 
				tp/(3600*24),
				tp%(3600*24)/3600,
				tp%3600/60,
				tp%60);
			_timeTxtM->SetInput(timeTxt);

		}
		this->Interactor->GetRenderWindow()->Render();
		vtkInteractorStyleTrackballCamera::OnTimer();
	}


}; 
vtkStandardNewMacro(ScadaGenVisInteratorStyle);


int main() {
	LPCTSTR inp = TEXT("Data\\c-town_true_network.inp");
	//LPCTSTR inp = TEXT("Data\\net3.inp");
	LPCTSTR dsn = TEXT(
		"DSN=jcx;\
		DESCRIPTION={Jinduan's Senm Database};\
		SERVER=localhost;\
		UID=rtx_db_agent;\
		PWD=rtx_db_agent;\
		DATABASE=rtx_demo_db;\
		PORT=3306;\
		FOUND_ROWS=1");
	Tstamp t0 = { //generation start time
	2013, 4, 22,  
	4, 7, 0, 0
	};
	int timespan = 3600*24*7;


	// make iid lognormal demands
	ScadaGenerator* psg = new ScadaGenerator();

	vtkPolyData* net[2];

	//psg->init(inp, NULL, 1234, net); // init, pull netvis info
	psg->init(inp, dsn, 1234, net);

	// create a varima model
	int nuser = psg->getNuser(); //388 in this example

	Varima* varima = new Varima(nuser, 2, 0);  //p, q, d
	//the network has 15-min demand interval so daily period = 24*4
	Varima* varima_diurnal = new Varima(varima, 0, 0, 24*4, 1); //P,Q,S,D
	Varima* varima_weekly = new Varima(varima_diurnal, 0, 0, 24*4*7, 1);

	double variance = 0.0003;  // variance of log(demand)
	double* phi_in = new double[nuser*2];
	double* phi_in2 = new double[nuser];
	double* phi_in3 = new double[nuser];
	double* cov = new double[nuser*nuser];

	Rvgs* rgp = new Rvgs(1111); //random number generator for parameters
	for (int i=0; i<nuser; ++i) {
		switch (i%3) {
		case 0: 
			phi_in[2*i] = 1.28;
			phi_in[2*i+1] = -0.4;
			break;
		case 1:
			phi_in[2*i] = 0.7;
			phi_in[2*i+1] = 0;
			break;
		case 2:
			phi_in[2*i] = 1;
			phi_in[2*i+1] = -0.5;
			break;
		}

		switch (i%2) {
		case 0:
			phi_in2[i]=-0.2;
			break;
		case 1:
			phi_in2[i]=0.12;
			break;
		}

		switch (i%5) {
		case 0: phi_in3[i] = 0.2; break;
		case 1: phi_in3[i] = -0.3; break;
		case 2: phi_in3[i] = rgp->Normal(0.1, 0.1); break;
		case 3: phi_in3[i] = -0.14; break;
		case 4: phi_in3[i] = rgp->Normal(0, 0.2); break;
		}

		for (int j=i; j<nuser; ++j) {
			if (j==i) cov[i*nuser + j] = cov[j*nuser+i] = variance;
			else 
				cov[i*nuser + j]=cov[j*nuser+i] = 
				  variance * rgp->Uniform(-0.5, 0.5);
		}
	}

	// test varima and generator
	varima->setPara(phi_in, NULL, cov);
	varima_diurnal->setPara(phi_in2, NULL, NULL);
	varima_weekly->setPara(phi_in3, NULL, NULL);


	// set up net skelecton 3d and 2d
	vtkPolyDataMapper* net3d = vtkPolyDataMapper::New();
	vtkPolyDataMapper* net2d = vtkPolyDataMapper::New();

	// set mapper to initial poly datasets
	net3d->SetInput(net[0]);
	net2d->SetInput(net[1]);

	vtkActor *actn3d = vtkActor::New();
	vtkActor *actn2d = vtkActor::New();
	actn3d->SetMapper(net3d);
	actn2d->SetMapper(net2d);

	// set up p and d visuals
	vtkPolyDataMapper* netvis[2]; //handle of data mapper
	netvis[0] = vtkPolyDataMapper::New(); //3d geom
	netvis[1] = vtkPolyDataMapper::New(); //2d geom
	vtkActor *pact = vtkActor::New();
	vtkActor *dact = vtkActor::New();
	pact->SetMapper(netvis[0]);
	pact->GetProperty()->SetColor(1, 1, 0);
	pact->GetProperty()->SetOpacity(0.5);
	dact->SetMapper(netvis[1]);
	dact->GetProperty()->SetColor(0.2, 0.5, 1);
	dact->GetProperty()->SetOpacity(1);

	// set up banner
	vtkTextMapper* banner = vtkTextMapper::New();
	banner->SetInput("Network");
	banner->GetTextProperty()->SetFontSize(14);

	vtkActor2D* banAct = vtkActor2D::New();
	banAct->SetMapper(banner);
	banAct->GetProperty()->SetColor(1, 1, 0);
	banAct->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
	banAct->GetPositionCoordinate()->SetValue(0.1, 0.9);

	// there is only one camera (for both renderer)
	vtkCamera *cam = vtkCamera::New();

	vtkRenderer *renp = vtkRenderer::New(); //3d view
	renp->AddActor(pact);
	renp->AddActor(actn3d);
	renp->AddActor2D(banAct);

	vtkRenderer *rend = vtkRenderer::New(); //2d flat view
	rend->AddActor(dact);
	rend->AddActor(actn2d);
	
	vtkRenderWindow *renWin = vtkRenderWindow::New();
	renWin->AddRenderer(renp); 
	renp->SetViewport(0, 0.5, 1, 1);
	renWin->AddRenderer(rend); 
	rend->SetViewport(0, 0, 1, 0.5);

	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	ScadaGenVisInteratorStyle* ctrl =  ScadaGenVisInteratorStyle::New();
	
	iren->SetInteractorStyle(ctrl);

	iren->SetRenderWindow(renWin);


	renp->SetActiveCamera(cam);
	renp->ResetCamera();
	renp->SetBackground(0,0,0);
	rend->SetActiveCamera(cam);
	rend->SetBackground(0,0,0);

	renWin->SetSize(800,600);

	// interact with data


	iren->Initialize();
	iren->CreateRepeatingTimer(50);

	ctrl->initPDdata(psg, varima_weekly, netvis, banner, &t0, timespan);
	renWin->Render();
	iren->Start();

	//netvis->Delete();
	pact->Delete();
	cam->Delete();
	renp->Delete();
	dact->Delete();

	rend->Delete();
	renWin->Delete();
	iren->Delete();
	delete psg;

	return 0;

}