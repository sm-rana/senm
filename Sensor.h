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


/** \file
	Visual widgets representing sensors/monitors in a distribution system.
*/
#pragma once
#include "vtkincludes.h"
#include "Network.h"
#include "DataSource.h"

///  Widget representing a scada channel/sensor in the vtk scenes.
/**	 Sensor class is the super class of all sensors. Currently
  * five types of sensors are provided as sub classes of Sensor:
  *	- Tank. water surface level (L) monitored, min and max possible level fixed
  *	- Reservoir/water source. level (L) monitored
  *	- Piezometer. Free water head (P) monitored
  *	- Flow rate meter. Flow rate (Q) monitored
  *	- Pump/valve on/off switch.  Status (C) monitored
*/

class Sensor {
public:

	/// Set common parameters for all sensors. Must be called before using the factory method.
	/** \param [in] net			Base vtk 3D network
	*	\param [in] z_scaling	z-axis scaling factor
	*	\param [in] unit		x,y-axis unit length for 3D visuals
	*/
	static void setCommonParameters(vtkPolyData* net, double z_scaling, double unit) ;

	/// Factory method for producing Sensor's sub class. Must call setCommonParameters() in advance.
	/** \param [in] cellid		Vtk cell id of where the sensor is located. may be a point or a line
	*	\param [in] chan		Pointer to the Channel struct for the Sensor
		\return The pointer to the produced sensor object, if the method fails, return NULL.
	*/
	static Sensor* makeSensor(vtkIdType cellid , Channel* chan) ;

	/// Get the vtk visuals.
	/** \return The visuals in assembly.
	*/
	virtual vtkPropAssembly* getGeom() {return _geom;};

	/// @{ 
	/** Set/get the sensor's data
	*/
	virtual void setData(double data_in) = 0; 
	virtual double getData() = 0;
	/// @}

	///@{
	/**Change geometry when selected/unselected. The operation will not change data. */
	virtual void select() =0 ;
	virtual void unselect() = 0;
	///@}


protected:
	vtkPropAssembly *_geom;  ///< geometry of the vis
	vtkCaptionActor2D* _caption;  ///< Context-based text caption/tooltip showing sensor information

	Channel* _channel; ///< the channel represented by the sensor
	Sensor():_geom(NULL) {};

	// parameters common to all sensors (for a network)
	static vtkPolyData* _net;  ///< 3D base network
	static double _z_scaling_factor; ///< z-axis scaling/exaggeration factor
	static double _xy_unit_len; ///< reference unit length for visuals
	static vtkLookupTable* _cmap; ///< color map for free head coloring
	
private:
	// copy operator and constructor are disabled
	Sensor(const Sensor&);
	void operator=(const Sensor&);

};


/// Widget of control device
/** A widget visualizing a on/off switchable controllable device.
	A SensorControl instance consists of a cylinder and 
	a single status variable (_on) indicating the 
	on/off status. The geometry is colored green/red for visualization
	*/

class SensorControl : public Sensor {
public:
	///Factory method
	/** \param [in] cellid		id of the line cell in the vtk network
		\param [in] chan		data channel for the sensor
		*/
	static SensorControl* make(vtkIdType cellid, Channel* chan);

	/// @{ 
	///Control status setter/getter
	virtual void setData(double status);
	virtual double getData() {return _on;};
	///@}

	/// @{
	/**Re-color the geometry when selected/unselected. The operation will not change data. */ 
	virtual void select();
	virtual void unselect();
	/// @}

protected:
	SensorControl() {};
	int _on;  ///> the on/off status
	double _pos1[3], _pos2[3]; ///> Start and end node coordinates
	vtkActor* _acCtrl; ///> visual for the control device
	
};

/// Widget of flow meter
/** A widget visualizing a flowrate measurement device.
	A SensorFlow instance consists of a collection of 3d arrows (shark fins),
	the number of which represents the relative flow rate.
	*/
class SensorFlow : public Sensor {
public:
	///Factory method
	/** \param [in] cellid		id of the line cell in the vtk network
		\param [in] chan		data channel for the sensor
		*/
	static SensorFlow* make(vtkIdType cellid, Channel* chan);

	/// @{ 
	///flow rate setter/getter
	virtual void setData(double flowrate);
	virtual double getData() {return _flowRate;};
	/// @}

	/// @{
	/**Re-color the geometry when selected/unselected. The operation will not change data. */ 
	virtual void select();
	virtual void unselect();
	/// @}

protected:
	SensorFlow() {};

	double _flowRate; ///> negative value means reverse flow (based on the link's direction)
	double _maxAbsQ; ///> maximum absolute flowrate
	double _pos1[3], _pos2[3]; ///> coordinates of the two ends
	vtkActor* _fins[4]; ///> Visuals (shark fins)

	static double finCoords[12][2]; ///> relative coordinates for the control points of a fin
	static vtkIdType finPtsId[6]; ///> Indices into finCoords[] for forward fins
	static vtkIdType finPtsIdR[6]; ///> Indices into finCoords[] for reverse fins
	static double helper[]; ///> Helper array for determining number of fins
	static double pi;
};


/// Widget of Pressure transducer
/** A widget visualizing a pressure head measurement device.
	A SensorHead instance consists of a tall cylinder,
	the height of which is in scale with the water head.
	the z-axis scaling of the cylinder is consistent with that of the
	elevations (coordinate z) used in the base map.
	*/
class SensorHead : public Sensor {
public:
	///Factory method
	/** \param [in] cellid		id of the line cell in the vtk network
		\param [in] chan		data channel for the sensor
		*/
	static SensorHead* make(vtkIdType cellid, Channel* chan);

	/// @{ 
	///Free head setter/getter
	virtual void setData(double newhead);
	virtual double getData() {return _freeHead;};
	/// @}

	/// @{
	/**Re-color the geometry when selected/unselected. The operation will not change data. */ 
	virtual void select();
	virtual void unselect();
	///@}

protected:
	SensorHead() {};
	double _freeHead;  ///> Free water head
	double _x, _y, _z;  ///> coordinates of the cylindric bottom.
	vtkActor* _waterHead;  ///> visual
	vtkFloatArray* _cellAtt;

};


/// Widget of Tanks
/** A widget visualizing a water storage facility.
	A SensorTank instance shows a vertical cylinder partly filled and partly transparent.
	The bottom and the top of the cylinder show the min/max water levels.
	the height of the filled part shows the current water level in the tank.
	the z-axis scaling of the cylinder is consistent with that of the
	elevations (coordinate z) used in the base map.
	*/
class SensorTank : public Sensor {
public:
	///Factory method
	/** \param [in] cellid		id of the line cell in the vtk network
		\param [in] chan		data channel for the sensor
		*/
	static SensorTank* make(
		vtkIdType cellid /*cell id where the tank is located*/,
		Channel* chan);

	/// @{ 
	///Current water level setter/getter
	virtual void setData(double) ; 
	virtual double getData() {return _waterLevel;};
	/// @}

	/// @{
	/**Re-color the geometry when selected/unselected. The operation will not change data. */ 
	virtual void select();
	virtual void unselect();
	///@}


protected:
	SensorTank() {};
	double _waterLevel; ///> Current water level
	
	/**  Min/max water levels of the tank. 
	They may be different from the elevation of the tank.
	   - (_lowlevel <= _z) represent a underground water tank
	   - (_lowlevel >= _z) represent a water tower
    */
	double _lowlevel, _highlevel; 
	
	double _x, _y, _z; ///> Tank coordinates
	vtkActor* _waterColumn; ///> The visual showing current water level
	
	vtkActor* _tankBody; ///> The visual showing the tank itself
};
