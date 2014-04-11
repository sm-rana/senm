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
#include "Network.h"
#include "DataSource.h"
#include "vtkincludes.h"

/// Visuals of EPANET network infrastructure 
/**The vtkPolyData-typed 2d and 3d network base map can be created. 
 Cell ID:  the cell id of the 2d 3d vtkPolyData for visualization. Starts from 0.
*/

class VisNetwork{
protected:
	Network * _net; // The epanet network
	///@{
	//look up table for net2d/3d cellid
	vtkIdType *_tabIdx2cellN;
	vtkIdType *_tabIdx2cellL; 
	///@}

	vtkPolyData* _net2d; //2d base map
	vtkPolyData* _net3d; //3d base map


public:
	///> use the singleton network to initialize 
	VisNetwork();

	///> destroy everything except for the network
	~VisNetwork();

	///> Create vtk base network
	/** build a 2d and a 3d vtkNetwork polydata, 
	net2d and net3d must have been vtk-New()-ed,
	the scalars of the output vtkPolyData (attribute) at nodes (vertex-typed cells) 
	are the baseline water demands
	\param [in,out]  net2d         Pointer to the 2D base network (no elevation info), 
	memory must be allocated before calling this function.
	\param[in, out]  net3d        Pointer to the 3d base network,
	memory must be allocated before calling this function.
	*/
	void get2d3dNet(vtkPolyData* net2d, vtkPolyData* net3d);

	///> Find the vtkId (cell id) of an EPANET index
	/** Lookup the net2d/3d graphic cell id from a EPANETnetwork index.
	Return -1 if the component can not be found.
	*/
	vtkIdType index2CellId(int index, Network::FieldType type) ;
	vtkIdType index2CellId(int index, Channel::Type type) ;

	///> Find the vtkId (cell id) of an EPANET id string (name).
	/**look up the graphic cell id from network component id (string)
	Return -1 if the component can not be found.
	*/
	vtkIdType netId2CellId(char* id, Network::FieldType type) ;
	vtkIdType netId2CellId(char* id, Channel::Type type) ;




}; 