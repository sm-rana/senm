#include <stdio.h>
#include "SenmCoreIncs.h"
/// Error reporting facilities

/// global log file
FILE* log_file = NULL;

void ewi(TCHAR const* error_text, int thread_no) {
	USES_CONVERSION;
	time_t rawtime;
	tm * timeinfo;
	char buffer[32];

	time (&rawtime);
	timeinfo = localtime (&rawtime);
	strftime(buffer, 32, "%X", timeinfo);
	_ftprintf(stderr, TEXT("%s [%d] %s\n"), A2T(buffer), thread_no, error_text);
	if (log_file != NULL) {
		_ftprintf(log_file, TEXT("%s [%d] %s\n"), A2T(buffer), thread_no, error_text);
	}

}

void init_ewi(TCHAR const * log_filepath) {
	log_file = _tfopen(log_filepath, TEXT("w"));
    if (log_file == NULL) 
		ewi(TEXT("Could not open log file"));

    TCHAR init_info[128];

    _stprintf(init_info, TEXT("Number of Workers: %d"), N_WORKERS);
    ewi(init_info);
}

void end_ewi() {
    if (log_file != NULL) {  
        ewi(TEXT("Close log file."));
        fclose(log_file);
	}
}

