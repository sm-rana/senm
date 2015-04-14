#pragma once

#include "SenMCoreIncs.h"
#include "Network.h"

struct Provider ;

/// Data channel information (metadata)
/** The struct stores descriptive information about a data channel.
Data channels must be created and destroyed by Provider::loadChannels(),
*/
struct Channel
{	
	int			key;		///> channel key (e.g., database primary key)

	/// Channel Types
	enum Type {
		NONE = 0, ///>No channel
		L = 0x1, ///> Water level
		F = 0x2, ///> Pump flow rate
		V = 0x4, ///> Minor headloss coefficient for TCV, FCV
		B = 0x8, ///> Discharge-side pressure
		A = 0x10, ///> Suction-side pressure
		Q = 0x20, ///> Pipe flow rate 
		C = 0x40, ///> Pipe/valve on/off
		P = 0x80, ///> Nodal pressure
		D = 0x100, ///> Real-time demand (flow meter)

	};
	Type    type;  ///> type of the channel

	/// return asset type in an EPANET model
	Network::FieldType  mtype() {
		switch (type) {
		case L: return Network::HEAD;  
		case V: return Network::SETTING;
		case B: case A: case P: return Network::PRESSURE;
		case C: return Network::LSTATUS;
		case F: case Q: case D: return Network::FLOW;
		default:
			return Network::UNKNOWN;
		}
	}

	int			mindex;  ///> index in the hydraulic network

	Provider	*provider;  ///> data provider

	char		name[MAX_NET_ID_LEN];  ///>  name of the component, consistent with network component id string

	//DataSource::UnitsType  unit; ///> unit of measurement, work is needed here
	//right now assume channel unit is consistent with network unit


	double		stde;  ///>  standard error of measurements (assuming normal)

	double		lower_lim;  ///> lower limits of the possible measurements
	double		upper_lim;  ///> upper limits

	double		error_default;   
	///>  if the the sensor has anormaly, measurements will default to this value

	enum Status {OK, DATA_ERR, DATA_OUTBOUND};
	Status		status;  ///>  working status of this channel

	Channel*	next;   ///> next channel in a list

	/// fetch the type string of a channel 
	/**  str must have MAX_COMP_TYPE_STR_SIZE tchars allocated */
	static TCHAR* ctype(Type ch_type, TCHAR* str) {
		switch (ch_type) {
		case A:
			_tcscpy(str, TEXT("Suction-side pressure of a pump/GPV"));
			break;
		case B:  
			_tcscpy(str, TEXT("Discharge-side pressure of a pump/GPV"));
			break;
		case C: 
			_tcscpy(str, TEXT("On/off status of a pipe/valve"));
			break;
		case D:
			_tcscpy(str, TEXT("Real-time water demand of a junction"));
			break;
		case F:
            _tcscpy(str, TEXT("Flow-rate of a pump"));
            break;
		case L: 
			_tcscpy(str, TEXT("Water levels in a tank/reservoir"));
			break;
		case Q:
			_tcscpy(str, TEXT("Flow rate inside a pipe"));
			break;
		case P: 
			_tcscpy(str, TEXT("Pressure at a junction"));
			break;
		case V: 
			_tcscpy(str, TEXT("Minor headloss coefficient of a TCV/FCV"));
			break;
		default:
			_tcscpy(str, TEXT("Unknown Component"));
		}
		return str;
	}


};

// for bitwise op to work...
inline Channel::Type operator|(Channel::Type a, Channel::Type b) {
	return static_cast<Channel::Type>(
		static_cast<unsigned>(a) | static_cast<unsigned>(b)
		);
}

inline Channel::Type operator&(Channel::Type a, Channel::Type b) {
	return static_cast<Channel::Type>(
		static_cast<unsigned>(a) & static_cast<unsigned>(b)
		);
}

inline Channel::Type& operator|=(Channel::Type& a, Channel::Type b) {
	return a = a|b;
}

inline Channel::Type& operator&=(Channel::Type& a, Channel::Type b) {
	return a = a&b;
}
