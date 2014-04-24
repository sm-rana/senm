#ifndef __vtkTextAnalysisInstantiator_h
#define __vtkTextAnalysisInstantiator_h

#include "vtkInstantiator.h"



class VTK_TEXT_ANALYSIS_EXPORT vtkTextAnalysisInstantiator
{
  public:
  vtkTextAnalysisInstantiator();
  ~vtkTextAnalysisInstantiator();
  private:
  static void ClassInitialize();
  static void ClassFinalize();
  static unsigned int Count;
}; 

static vtkTextAnalysisInstantiator vtkTextAnalysisInstantiatorInitializer;

#endif
