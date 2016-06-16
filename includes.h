#include <list>
#include <cmath>
#include <cstdio>
#include <random>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <algorithm>

#include <vtkGlyph3D.h>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkImageFlip.h>
#include <vtkDataArray.h>
#include <vtkContourFilter.h>
#include <vtkDataObject.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>
#include <vtkPolyDataNormals.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkStructuredPointsReader.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkContourFilter.h>
#include <vtkDoubleArray.h>
#include <vtkTIFFReader.h>
#include <vtkTIFFWriter.h>
#include <vtkPointData.h>
#include <vtkPoints.h>

#include "_cc_finder.h"
#include "_database.h"
#include "_eblob.h"

#define TIFF_ORIENTATION_READER 1
#define DEBUG
