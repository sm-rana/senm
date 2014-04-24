#include "SenmCoreIncs.h"
#include "SenmVtkIncs.h"
#include "ScadaGenerator.h"

#define VSG_N_NETPOLY 2
#define VSG_TIME_BANNER_SIZE 512

/// 3-d visualized SCADA generator
/** Stores and manages the various vtk components*/
struct VisScadaGenerator {
	ScadaGenerator* _sg; ///>internal scada generator
    int _timespan; ///>total generation time in seconds

	vtkRenderWindow* _renWin;
	vtkRenderWindowInteractor* _iren;

    // 3d [0] and 2d [1] base network polydatas
	vtkPolyData* _baseNet[VSG_N_NETPOLY];

	vtkPolyData** _netStore; //polydata storage for p and d, one in each time step
	int* _dispStep; // displayed elapsed time - array
    int _nNetStore; ///> size of _netStore

	volatile LONG* _datastep; ///> pointer to most recent step of p/d availability

    /// unit z-axis line segment
	vtkPolyData *zlineP, *zlineD;

    // Visual limits
    double minx, miny, minz, maxx, maxy, maxz;

    // Mappers for p and d visuals
	vtkPolyDataMapper *_mVisP, *_mVisD; 

    // text banner 
	vtkTextMapper* _timeTxtM; 
	char timeTxt[VSG_TIME_BANNER_SIZE];

};

/// visual scada generator errors
enum VSG_ERR {
    VSG_OK,
    VSG_MEM_NOT_ALLOCED,
    VSG_SG_ERR,

	VSG_DUMMY_LAST};

/// VSG error reporting
void VSG_ewi(VSG_ERR err);

/// VSG object factory
VSG_ERR VSG_new(TCHAR* inpfile_path, char const* ini_path, int seed, Varima* md, int timespan, VisScadaGenerator** vsg_out);

/// VSG visual update function - to be assigned to ScadaGenerator's extUpdate pointer
void VSG_vupdate(void* vsg, ScadaGenerator* sg);

/// Move P and D visuals to a time step
void VSG_moveto(VisScadaGenerator* vsg, int facestep);

/// VSG run: run scada generator with real-time 3d visuals
void VSG_run(VisScadaGenerator* vsg, Tstamp* t_start);

/// destroy
void VSG_delete(VisScadaGenerator* vsg);


///InteratorStyle for ScadaGenerator
class ScadaGenVisInteratorStyle: public vtkInteractorStyleTrackballCamera {
	
protected:
	volatile long datastep; // where generator has done its work and v/d visuals available
	long facestep; // what is shown currently
	long desiredstep; // desired step to be shown

    VisScadaGenerator* vsg;

    //hidden constructor
	ScadaGenVisInteratorStyle():
		datastep(-1), facestep(-1), desiredstep(-1), vsg(NULL),
		vtkInteractorStyleTrackballCamera(){};
public:
	static ScadaGenVisInteratorStyle* New();
	vtkTypeMacro(ScadaGenVisInteratorStyle, vtkInteractorStyleTrackballCamera);

	void setVSG(VisScadaGenerator* vsg_in); 

	virtual void OnKeyPress(); 
	virtual void OnTimer(); 

}; 

