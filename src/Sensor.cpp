#include <Windows.h>
#include <tchar.h>
#include <math.h>
#include "SenmVtkIncs.h"
#include "Sensor.h"

//default values
//
//
vtkPolyData* Sensor::_net = NULL;
double Sensor::_z_scaling_factor = 0;
double Sensor::_xy_unit_len = 0;
vtkLookupTable* Sensor::_cmap = NULL;

void Sensor::setCommonParameters(vtkPolyData* net, double z_scaling, double unit) {
	if (net == NULL ) {
		_ftprintf(stderr, TEXT("Base network not provided for sensor building.\n"));
		return;
	}

	if (z_scaling <= 0 || unit <= 0) {
		_ftprintf(stderr, TEXT("Illegal parameters for sensor factory.\n"));
	}
	_net = net;
	_z_scaling_factor = z_scaling;
	_xy_unit_len = unit;
	
	vtkLookupTable* _cmap = vtkLookupTable::New();
	_cmap->SetNumberOfColors(64);
	//_cmap->SetHueRange(0.001, 108.0/360); // red to green
	_cmap->SetHueRange(0.999, 248/360);  //red to blue
	_cmap->Build();

}


Sensor* Sensor::makeSensor(vtkIdType cellid, Channel* chan) 
	{
		if (_net == NULL || _z_scaling_factor == 0 || _xy_unit_len == 0) 
			return NULL;
		
		vtkIdType nCell = _net->GetNumberOfCells();
		if (cellid < 0 || cellid >= nCell) {			
			_ftprintf(stderr, TEXT("Cellid out of range. can't build sensor.\n"));
			return NULL;
		}

		switch (chan->type) {
		case Channel::L:
			return SensorTank::make(cellid, chan);
		case Channel::P:
			return SensorHead::make(cellid, chan);
		case Channel::Q:
			return SensorFlow::make(cellid, chan);
		case Channel::C:
			return SensorControl::make(cellid, chan);
		default:
			return NULL;
		}
	};


SensorTank* SensorTank::make(vtkIdType cellid, Channel* chan) {
	//produce the geometry for water source visualization
	
	double ll = chan->lower_lim; 
	double ul = chan->upper_lim;

	SensorTank* re = new SensorTank();
	re->_geom = vtkPropAssembly::New();
	re->_waterColumn = vtkActor::New();
	re->_caption = vtkCaptionActor2D::New();

	re->_lowlevel = ll;
	re->_highlevel = ul;

	//position of the sensor
	double bounds[6];
	_net->GetCellBounds(cellid, bounds);
	re->_x = bounds[0]; 
	re->_y = bounds[2]; 
	re->_z = bounds[4];

	//re->_caption->GetCaptionTextProperty()->SetFontSize(14);
	//re->_caption->GetCaptionTextProperty()->SetFontFamilyToCourier();
	re->_caption->GetCaptionTextProperty()->ItalicOff();
	re->_caption->GetCaptionTextProperty()->BoldOff();
	re->_caption->GetCaptionTextProperty()->ShadowOff();
	re->_caption->GetCaptionTextProperty()->SetLineSpacing(1.1);
	re->_caption->BorderOff();
	re->_caption->SetHeight(0.05);
	re->_caption->LeaderOff();
	re->_caption->AttachEdgeOnlyOn();
	re->_caption->SetPosition(-90,6);
	re->_channel = chan;

	re->_geom->AddPart(re->_waterColumn);
	re->_geom->AddPart(re->_caption);

	if (ll == 0 && ul == 0) {//water source - no tank body

		vtkCubeSource *cs = vtkCubeSource::New();
		vtkPolyDataMapper *pm = vtkPolyDataMapper::New();
		pm->SetInput(cs->GetOutput());
		re->_waterColumn->SetMapper(pm);

		pm->ScalarVisibilityOff();

		//temporary value for water column
		re->_waterColumn->SetScale(
			_xy_unit_len*8, _xy_unit_len*8, _xy_unit_len*8);
		re->_waterColumn->SetPosition(
			re->_x, re->_y, re->_z +_xy_unit_len*2); 
	} else { //tank
		
		//tank shell
		vtkCylinderSource* ys = vtkCylinderSource::New();
		ys->SetResolution(10);
		re->_tankBody = vtkActor::New();
		vtkPolyDataMapper *ym = vtkPolyDataMapper::New();

		ym->SetInput(ys->GetOutput());
		ym->Update();
		ym->ScalarVisibilityOff();
		re->_tankBody->SetMapper(ym);

		re->_tankBody->SetScale(_xy_unit_len*8, 
			ul*_z_scaling_factor - re->_z, _xy_unit_len*8);
		re->_tankBody->SetOrientation(90,0,0);
		re->_tankBody->SetPosition(re->_x, re->_y, 
			(ul*_z_scaling_factor + re->_z)/2);

		re->_tankBody->GetProperty()->SetColor(1, 1, 0);
		re->_tankBody->GetProperty()->SetOpacity(0.4);

		re->_geom->AddPart(re->_tankBody);

		//water column
		vtkCylinderSource *cs = vtkCylinderSource::New();
		vtkPolyDataMapper *pm = vtkPolyDataMapper::New();
		cs->SetCapping(1);
		cs->SetResolution(20);
		cs->SetCenter(0, 0.5, 0);
		pm->SetInput(cs->GetOutput());
		re->_waterColumn->SetMapper(pm);
		re->_waterColumn->RotateX(90);

		//temporary value for water column
		re->_waterColumn->SetScale(
			_xy_unit_len*4, _xy_unit_len*4, _xy_unit_len*4);
		re->_waterColumn->SetPosition(
			re->_x, re->_y, re->_lowlevel * _z_scaling_factor); 

		pm->ScalarVisibilityOff();
		
	}

	re->_waterLevel = 0; //temporary
	

	re->_waterColumn->GetProperty()->SetColor(0,0,1); //blue
	re->_geom->VisibilityOff();

	return re;
}

void SensorTank::setData(double newlevel) {
	// set up a new water level
	
	if (_highlevel != 0 || _lowlevel != 0) {//tank
	_waterColumn->SetScale(
		_xy_unit_len*7.9, (newlevel-_lowlevel)*_z_scaling_factor , _xy_unit_len*7.9);
	}

	char string[128];
	sprintf(string, "[%s] %2.f%% full\n%7.2f (%+3.1f)", _channel->name,
		(newlevel-_lowlevel)/(_highlevel-_lowlevel)*100,
		 newlevel, newlevel-_waterLevel);
	_caption->SetCaption(string);
	_caption->SetAttachmentPoint(_x, _y, newlevel*_z_scaling_factor);
	
	_waterLevel = newlevel;
	_geom->VisibilityOn();
	
}

void SensorTank::select() {
	_waterColumn->GetProperty()->SetColor(0, 1, 1); //cyan
	_geom->VisibilityOn();
}

void SensorTank::unselect() {
	_waterColumn->GetProperty()->SetColor(0, 0, 1); //blue
	_geom->VisibilityOn();
}

SensorHead* SensorHead::make(vtkIdType cellid, Channel* chan) {
	
	SensorHead* re = new SensorHead();
	re->_channel = chan;
	re->_freeHead = 0; //temp
	re->_geom = vtkPropAssembly::New();
	re->_waterHead = vtkActor::New();
	
	//set up caption
	re->_caption = vtkCaptionActor2D::New();
	re->_caption->GetCaptionTextProperty()->ItalicOff();
	re->_caption->GetCaptionTextProperty()->BoldOff();
	re->_caption->GetCaptionTextProperty()->ShadowOff();
	re->_caption->GetCaptionTextProperty()->SetLineSpacing(1.1);
	re->_caption->BorderOff();
	re->_caption->SetHeight(0.028);
	re->_caption->SetPosition(0,0);
	re->_caption->LeaderOff();
	
	re->_caption->PickableOff();
	re->_geom->AddPart(re->_caption);

	//position of the sensor
	double bounds[6];
	_net->GetCellBounds(cellid, bounds);
	re->_x = bounds[0]; 
	re->_y = bounds[2]; 
	re->_z = bounds[4];

	// set up a cylinder showing free head
	vtkCylinderSource *cs = vtkCylinderSource::New();
		cs->SetCenter(0, 0.5, 0);
		cs->CappingOn();
		cs->SetResolution(3);
	vtkPolyData* tmp = cs->GetOutput();
	tmp->Update();
	re->_cellAtt = vtkFloatArray::New();
	re->_cellAtt->SetNumberOfTuples(tmp->GetNumberOfCells());
	for (vtkIdType i = 0; i<tmp->GetNumberOfCells(); ++i) {
		re->_cellAtt->InsertTuple1(i, 0);  //init the array for cell attribute
	}
	tmp->GetCellData()->SetScalars(re->_cellAtt);
	vtkPolyDataMapper *pm = vtkPolyDataMapper::New();
		pm->SetInput(tmp);
		pm->SetLookupTable(_cmap);
		pm->UseLookupTableScalarRangeOff();
		pm->SetScalarRange(0, 100); //the max value should consider channel's unit type in the future
		pm->SetColorModeToMapScalars();
		pm->SetScalarModeToUseCellData();

	re->_waterHead->SetMapper(pm);
	re->_waterHead->RotateX(90);
	re->_waterHead->SetPosition(re->_x, re->_y, re->_z);

	re->_geom->AddPart(re->_waterHead);
	re->_geom->VisibilityOff();

	return re;
}

void SensorHead::setData(double newhead) {
	
	if (newhead != 0) {//update geometry and color
		_waterHead->SetScale(
			_xy_unit_len*1.7, 
			_z_scaling_factor*newhead,
			_xy_unit_len*1.7);
		for (vtkIdType i = 0; i<_cellAtt->GetNumberOfTuples(); ++i) {
			_cellAtt->SetTuple1(i, newhead);
		}
		_geom->VisibilityOn();
	}
	
	char string[128];
	sprintf(string, "[%s] %7.2f (%+3.1f)", _channel->name,
		 newhead, newhead-_freeHead);
	_caption->SetCaption(string);
	_caption->SetAttachmentPoint(_x, _y, _z + newhead*_z_scaling_factor);
	
	if (_freeHead == 0) {// first time set data
		_caption->VisibilityOff();
	}
	
	_freeHead = newhead;
	
}

void SensorHead::select() {
	_waterHead->GetMapper()->ScalarVisibilityOff();
	_waterHead->GetProperty()->SetColor(0, 1, 1);
	_caption->VisibilityOn();
	_geom->VisibilityOn();
}

void SensorHead::unselect() {
	_waterHead->GetMapper()->ScalarVisibilityOn();
	_caption->VisibilityOff();
	_geom->VisibilityOn();
}


SensorControl* SensorControl::make(vtkIdType cellid, Channel* chan) {
	SensorControl* re = new SensorControl();

	re->_channel = chan;
	re->_on = 1; //temp
	re->_geom = vtkPropAssembly::New();
	
	//set up caption
	re->_caption = vtkCaptionActor2D::New();
	re->_caption->GetCaptionTextProperty()->ItalicOff();
	re->_caption->GetCaptionTextProperty()->BoldOff();
	re->_caption->GetCaptionTextProperty()->ShadowOff();
	re->_caption->BorderOff();
	re->_caption->SetHeight(0.028);
	re->_caption->SetPosition(3,-10);
	re->_caption->VisibilityOff();
	re->_caption->LeaderOff();
	
	re->_caption->PickableOff();
	re->_geom->AddPart(re->_caption);

	//get the two points' positions
	vtkIdType npts, *ptsid;
	_net->GetCellPoints(cellid, npts, ptsid);
	_net->GetPoint(ptsid[0], re->_pos1);
	_net->GetPoint(ptsid[1], re->_pos2);

	//build geometry, in the future fancier widgets can be added
	vtkLineSource* ls = vtkLineSource::New();
		ls->SetPoint1((float)(re->_pos1[0]), (float)(re->_pos1[1]),(float)(re->_pos1[2]));
		ls->SetPoint2((float)(re->_pos2[0]), (float)(re->_pos2[1]),(float)(re->_pos2[2]));

	vtkTubeFilter* tf = vtkTubeFilter::New();
		tf->SetInput(ls->GetOutput());
		tf->SetRadius(_xy_unit_len*0.4);
		tf->CappingOn();
		tf->SetNumberOfSides(6);

	vtkPolyDataMapper* pm = vtkPolyDataMapper::New();
		pm->SetInput(tf->GetOutput());
		pm->ScalarVisibilityOff();

	re->_acCtrl = vtkActor::New();
		re->_acCtrl->SetMapper(pm);
		re->_acCtrl->GetProperty()->SetColor(0, 1, 0);//default green

	re->_geom->AddPart(re->_acCtrl);
	re->_geom->VisibilityOff();

	return re;

}

void SensorControl::setData(double newstatus) {
	if (newstatus == 1) {
		_acCtrl->GetProperty()->SetColor(0, 1, 0);
	} else if (newstatus == 0) {
		_acCtrl->GetProperty()->SetColor(1, 0, 0);
	}

	_on = newstatus;

	char string[128];
	sprintf(string, "[%s] %s", _channel->name, _on?"On":"Off");
	_caption->SetCaption(string);
	_caption->SetAttachmentPoint(
		0.9*_pos1[0]+0.1*_pos2[0],
		0.9*_pos1[1]+0.1*_pos2[1],
		0.9*_pos1[2]+0.1*_pos2[2]
		);

	_geom->VisibilityOn();

}

void SensorControl::select() {
	_caption->VisibilityOn();
}

void SensorControl::unselect() {
	_caption->VisibilityOff();
}

double SensorFlow::finCoords[12][2] = 
{{0,0}, {-1, 2}, {-2, 2}, {-1,0}, {-2,-2}, {-1, -2}, //
 {0,0}, {1, 2}, {2, 2}, {1,0}, {2,-2}, {1, -2} };  //fin reversed

vtkIdType SensorFlow::finPtsId[] = {0, 1, 2, 3, 4, 5};
vtkIdType SensorFlow::finPtsIdR[] = {6, 7, 8, 9, 10, 11};

double SensorFlow::pi = 3.141592653589793238463;

SensorFlow* SensorFlow::make(vtkIdType cellid, Channel* chan) {
	SensorFlow* re = new SensorFlow();

	re->_channel = chan;
	re->_flowRate = 0; //temp value
	if (chan->lower_lim == 0 && chan->upper_lim == 0) {
		re->_maxAbsQ = 200; //default max. needs to consider unit in the future
	} else { //set max abs value
		re->_maxAbsQ = abs(chan->lower_lim) > abs(chan->upper_lim)? 
			abs(chan->lower_lim):abs(chan->upper_lim);
	}

	re->_geom = vtkPropAssembly::New();
	
	//set up caption
	re->_caption = vtkCaptionActor2D::New();
	re->_caption->GetCaptionTextProperty()->ItalicOff();
	re->_caption->GetCaptionTextProperty()->BoldOff();
	re->_caption->GetCaptionTextProperty()->ShadowOff();
	re->_caption->BorderOff();
	re->_caption->SetHeight(0.028);
	re->_caption->SetPosition(3,10);
	re->_caption->VisibilityOff();
	re->_caption->LeaderOff();
	
	re->_caption->PickableOff();
	re->_geom->AddPart(re->_caption);

	//get the two points' positions
	vtkIdType npts, *ptsid;
	_net->GetCellPoints(cellid, npts, ptsid);
	_net->GetPoint(ptsid[0], re->_pos1);
	_net->GetPoint(ptsid[1], re->_pos2);

	//build shark fins
	vtkPolyData* aFin = vtkPolyData::New();
	vtkPoints* finPts = vtkPoints::New();
		
		for (int i =0;i<12;++i) {
			finPts->InsertPoint(i, finCoords[i][0], finCoords[i][1], 0);
		}
	
	vtkCellArray* finCell = vtkCellArray::New();
		finCell->InsertNextCell(6, finPtsId); //has only one cell
	aFin->SetPoints(finPts); finPts->Delete();
	aFin->SetPolys(finCell); finCell->Delete();

	for (int i=0; i<4; ++i) {
		re->_fins[i] = vtkActor::New();
		vtkPolyDataMapper *mp = vtkPolyDataMapper::New();
			mp->SetInput(aFin);
			mp->ScalarVisibilityOff();
		re->_fins[i]->SetMapper(mp);
		re->_fins[i]->SetScale(_xy_unit_len);

		//line vector (vx, vy, vz) shows the direction
		double vx = re->_pos2[0] - re->_pos1[0];
		double vy = re->_pos2[1] - re->_pos1[1];
		double vz = re->_pos2[2] - re->_pos1[2];
		double norm2 = sqrt(vx*vx + vy*vy + vz*vz);

		//rotate the fin to orient the right direction
		double roty = atan(vz/sqrt(vx*vx+vy*vy)) * 180/pi * (vz>0?-1:1); //vtk RotateY seems to have a +/- problem
		double rotz = atan(vy/vx) * 180/pi + (vx<0?180:0);
		re->_fins[i]->RotateY(roty);
		re->_fins[i]->RotateZ(rotz);

		//line center
		double cx = (re->_pos2[0] + re->_pos1[0])/2;
		double cy = (re->_pos2[1] + re->_pos1[1])/2;
		double cz = (re->_pos2[2] + re->_pos1[2])/2;

		re->_fins[i]->SetPosition(cx + vx/norm2 * (-2+2*i)*_xy_unit_len, 
			cy + vy/norm2 * (-2+2*i)*_xy_unit_len, 
			cz + vz/norm2 * (-2+2*i)*_xy_unit_len); //shift the 4 fins

		re->_fins[i]->GetProperty()->SetColor(1,1,0); //yellow
		re->_fins[i]->GetProperty()->SetLineWidth(2);
		re->_geom->AddPart(re->_fins[i]);
	}

	re->_geom->VisibilityOff();
	
	return re;
}

double SensorFlow::helper[] = {0, 2, 5, 9, 14}; //helper[i] is the max number classes representable by i fins

void SensorFlow::setData(double flowrate) {
	for (int i = 0; i<4; ++i) {
		vtkCellArray* finCell = vtkCellArray::New();
		if (flowrate < 0) {//reverse flow
			finCell->InsertNextCell(6, finPtsIdR); 
		} else { //normal
			finCell->InsertNextCell(6, finPtsId);
		}
		vtkPolyData::SafeDownCast(_fins[i]->GetMapper()->GetInput())->SetPolys(finCell);
		finCell->Delete();

		_fins[i]->GetProperty()->SetRepresentationToWireframe();
		_fins[i]->VisibilityOff();
		
	}

	int nBox = abs(flowrate)/_maxAbsQ * 14; //which level
	int i;
	for (i=0; i<4; ++i) {
		if (nBox > helper[i]) {
			_fins[i]->VisibilityOn();
		} else {
			break;
		}
	}
	if (i>0) for (int j=0; j<nBox-helper[i-1]-1; ++j) {
		_fins[j]->GetProperty()->SetRepresentationToSurface();
	}

	char string[128];
	sprintf(string, "[%s] %7.2f (%+3.1f)", _channel->name, 
		abs(flowrate), abs(flowrate-_flowRate));
	_caption->SetCaption(string);
	_caption->SetAttachmentPoint(
		0.5*_pos1[0]+0.5*_pos2[0],
		0.5*_pos1[1]+0.5*_pos2[1],
		0.5*_pos1[2]+0.5*_pos2[2]
		);
	
	_flowRate = flowrate;
	_geom->VisibilityOn();
}

void SensorFlow::select() {
	_caption->VisibilityOn();
	for (int i = 0; i < 4; ++i ) {
		_fins[i]->GetProperty()->SetColor(0, 1, 1);
	}

}

void SensorFlow::unselect() {
	_caption->VisibilityOff();
	for (int i = 0; i < 4; ++i ) {
		_fins[i]->GetProperty()->SetColor(1, 1, 0);
	}
}
