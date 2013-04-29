#include "Viewer.h"

Viewer* Viewer::singleton = NULL;

Viewer* Viewer::getViewer(Viewer::UItype uit) {
	if (!singleton) 
		singleton = new Viewer;

	return singleton;
}

Viewer::Viewer() {
}

Viewer::~Viewer() {
}
