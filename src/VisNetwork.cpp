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


#include "VisNetwork.h"

VisNetwork::VisNetwork(): 
		_tabIdx2cellN (NULL), _tabIdx2cellL(NULL),
        _net(NULL), _net2d(NULL), _net3d(NULL) 
{
	Network::getNetwork(NULL, &_net);
    if (_net == NULL) return;  // Network not created

	double coords[3];

	_tabIdx2cellN = new vtkIdType[_net->MaxNodes+1];
	_tabIdx2cellL = new vtkIdType[_net->MaxLinks+1];

	vtkPoints* pts3 = vtkPoints::New();
	vtkPoints* pts2 = vtkPoints::New();
	vtkCellArray* ca = vtkCellArray::New();  // cells for nodes
	vtkCellArray* cl = vtkCellArray::New();  // cells for links
	vtkFloatArray* scalar_data = vtkFloatArray::New();
	vtkPolyData* tpn1 = vtkPolyData::New();

	int i, iUsers;
	for (i=1, iUsers=0; i<=_net->MaxNodes; ++i) {
		coords[0] = _net->Node[i].x;
		coords[1] = _net->Node[i].y;
		coords[2] = _net->Node[i].El * _net->Ucf[Network::ELEV];

		pts3->InsertPoint(i-1, coords);
		vtkIdType icell = i-1;
		_tabIdx2cellN[i] = ca->InsertNextCell(1, &icell);

		// 2d
		coords[2] = 0;
		pts2->InsertPoint(i-1, coords);

		if (_net->Node[i].D == NULL) {
			scalar_data->InsertTuple1(i-1, 0);
		} else {
			scalar_data->InsertTuple1(i-1, _net->Node[i].D->Base);
//			if (Node[i].D->Base >0) {// is it a water user?
//				iUsers++;
//			}
		}
	}

	for (i=1; i<=_net->MaxLinks; ++i) {
		vtkIdType jcell[2];
		jcell[0] = _net->Link[i].N1 - 1; // to vis index
		jcell[1] = _net->Link[i].N2 - 1;
		_tabIdx2cellL[i] = cl->InsertNextCell(2, jcell) + 
			ca->GetNumberOfCells(); 
		// this is very important to make each line cell id 
		// right when using the table to fetch
	}

	_net3d = vtkPolyData::New();
	_net2d = vtkPolyData::New();

	if(_net3d) {
		_net3d->SetPoints(pts3);
		_net3d->SetVerts(ca);
		_net3d->SetLines(cl);
		_net3d->GetPointData()->SetScalars(scalar_data);
		_net3d->BuildCells();
	}
	if (_net2d) {
		_net2d->SetPoints(pts2);
		_net2d->SetVerts(ca);
		_net2d->SetLines(cl);
		_net2d->GetPointData()->SetScalars(scalar_data);
		_net2d->BuildCells();
	}
	pts3->Delete();
	pts2->Delete();
	ca->Delete();
	cl->Delete();
	scalar_data->Delete();

}

VisNetwork::~VisNetwork() {
	delete[] _tabIdx2cellN;
	delete[] _tabIdx2cellL;
	_net2d->Delete();
	_net3d->Delete();
}

void VisNetwork::get2d3dNet(vtkPolyData* net2d, vtkPolyData* net3d) {
	if (net2d) net2d->DeepCopy(_net2d);
	if (net3d) net3d->DeepCopy(_net3d);
}

vtkIdType VisNetwork::index2CellId(int index, Channel::Type type) {
	//lookup cell ids
	if (_net == NULL) return -1;
	switch (type) {
	case Channel::L:
	case Channel::P:
	case Channel::D:
		if (index >0 && index <= _net->MaxNodes) 
			return _tabIdx2cellN[index];
		break;
	default:
		if (index >0 && index <= _net->MaxLinks) 
			return _tabIdx2cellL[index];
	}
	return -1;
}

vtkIdType VisNetwork::index2CellId(int index, Network::FieldType type) {
	//lookup cell ids
	if (_net == NULL) return -1;
	if (type <= Network::PRESSURE) {//lookup node cellid
		if (index >0 && index <= _net->MaxNodes) 
			return _tabIdx2cellN[index];
	} else {
		if (index >0 && index <= _net->MaxLinks) 
			return _tabIdx2cellL[index];
	}
	return -1;
}


vtkIdType VisNetwork::netId2CellId(char* netId, Network::FieldType type) {
	//lookup cell id
	if (type <= Network::PRESSURE) {
		Network::HTIt jit = _net->Nht.find(netId);
		if (jit == _net->Nht.end() ) return -1; //can't find the id
		return index2CellId(jit->second, type);
	} else {
		Network::HTIt jit = _net->Lht.find(netId);
		if (jit == _net->Lht.end() ) return -1; //can't find the id
		return index2CellId(jit->second, type);
	}
	return -1;
}
vtkIdType VisNetwork::netId2CellId(char* netId, Channel::Type type) {
	//lookup cell id
	Network::HTIt jit;
	switch (type) {
	case Channel::L:
	case Channel::P:
	case Channel::D:
	case Channel::B:
	case Channel::A:
		jit = _net->Nht.find(netId);
		if (jit == _net->Nht.end() ) return -1; //can't find the id
		return index2CellId(jit->second, type);
	default:
		jit = _net->Lht.find(netId);
		if (jit == _net->Lht.end() ) return -1; //can't find the id
		return index2CellId(jit->second, type);
	}
	return -1;
}



