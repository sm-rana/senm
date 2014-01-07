#pragma once
#include "vtkincludes.h"
#include "Network.h"
#include "DataSource.h"

///  A widget representing a scada channel/sensor in the vtk scenes.
/**	 Sensor class is the super class of all sensors. Currently
  * five types of sensors are provided as sub class of Sensor:
  *	- Tank. water surface level (L) monitored, min and max possible level fixed
  *	- Reservoir/water source. level (L) monitored
  *	- Piezometer. Free water head (P) monitored
  *	- Flow rate meter. Flow rate (Q) monitored
  *	- Pump/valve on/off switch.  Status (C) monitored
*/

class Sensor {
public:

	/// Set common parameters for all sensors.
	/** \param [in] Base vtk 3D network
	*	\param [in] z-axis scaling factor
	*	\param [in] x,y-axis unit length for 3D visuals
	*/
	static void setCommonParameters(vtkPolyData* net, double z_scaling, double unit) ;

	/// Factory method for producing Sensor's sub class
	/** \param [in] Vtk cell id of where the sensor is located. can be point or line
	*	\param [in] Pointer to the Channel struct for the Sensor
	*/
	static Sensor* makeSensor(vtkIdType cellid , Channel* chan) ;

	//get the geometry for viewing
	virtual vtkPropAssembly* getGeom() {return _geom;};

	//set/get the sensor data and update geometry
	virtual void setData(double) = 0; 
	virtual double getData() = 0;

	//select/unselect
	virtual void select() =0 ;
	virtual void unselect() = 0;


protected:
	vtkPropAssembly *_geom;  //geometry of the vis
	vtkCaptionActor2D* _caption;  // text caption/tooltip

	Channel* _channel; //the channel the sensor represent
	Sensor():_geom(NULL) {};

	// parameters common to all sensors (for a network)
	static vtkPolyData* _net;
	static double _z_scaling_factor; //z-axis scaling/exaggeration factor
	static double _xy_unit_len; //unit length for visuals
	static vtkLookupTable* _cmap; //color map for free head coloring
	
private:
	// copy operator and constructor are disabled
	Sensor(const Sensor&);
	void operator=(const Sensor&);

};

class SensorControl : public Sensor {
public:
	//factory
	static SensorControl* make(vtkIdType cellid, Channel* chan);

	virtual void setData(double status);
	virtual double getData() {return _on;};

	//hightlight the geometry
	virtual void select();
	virtual void unselect();
protected:
	SensorControl() {};
	int _on;  //the status
	double _pos1[3], _pos2[3]; //node coordinates
	vtkActor* _acCtrl;
	
};

class SensorFlow : public Sensor {
public:
	//factory
	static SensorFlow* make(vtkIdType cellid, Channel* chan);

	virtual void setData(double flowrate);
	virtual double getData() {return _flowRate;};

	virtual void select();
	virtual void unselect();

protected:
	SensorFlow() {};
	double _flowRate; //negative value means reverse flow
	double _maxAbsQ; //maximum absolute flowrate
	double _pos1[3], _pos2[3]; //ends
	vtkActor* _fins[4]; //shark fins
	static double finCoords[12][2]; //fin coordinates
	static vtkIdType finPtsId[6];
	static vtkIdType finPtsIdR[6]; //reverse fin cell id
	static double helper[]; //helper array for fin ranking
	static double pi;
};

class SensorHead : public Sensor {
public:
	//factory
	static SensorHead* make(vtkIdType cellid, Channel* chan);

	virtual void setData(double newhead);
	virtual double getData() {return _freeHead;};

	//hightlight the geometry
	virtual void select();
	virtual void unselect();
protected:
	SensorHead() {};
	double _freeHead;
	double _x, _y, _z; //node coordinates
	vtkActor* _waterHead;
	vtkFloatArray* _cellAtt;

};

class SensorTank : public Sensor {
public:
	static SensorTank* make(
		vtkIdType cellid /*cell id where the sensor is located*/,
		Channel* chan);

	//set/get the sensor data and update geometry
	virtual void setData(double) ; 
	virtual double getData() {return _waterLevel;};

	//hightlight the geometry
	virtual void select();
	virtual void unselect();

protected:
	SensorTank() {};
	double _waterLevel; //water level of the Tank
	double _lowlevel, _highlevel; // bottom elevation of the tank (may be different from the elevation of the node)
	double _x, _y, _z; //nodal elevation
	vtkActor* _waterColumn;
	
	vtkActor* _tankBody;
};
