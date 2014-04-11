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


#include <Windows.h>

class Viewer {
	// the singleton class for the async viewer

public:
	enum UItype {
		// the types of UI to be used, can be combined
		MATLAB = 0x1, 
		NATIVE_W32 = 0x2
	};

	//get the viewer
	static Viewer* getViewer(UItype uit);

	//show message
	void msg(const char*);

	~Viewer();


private:
	Viewer();
	UItype uitype;
	static Viewer* singleton; // the singleton instance

};