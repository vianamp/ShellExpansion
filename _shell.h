#ifndef _SHELL_H
#define _SHELL_H

#include "includes.h"

    struct _voxel {
        vtkIdType i;
        unsigned int x;
        unsigned int y;
        unsigned int z;
        double w;
    };

    class _shell {

        private:

            std::vector<_voxel> Voxels;

        public:

            void Clear();
            void Decompose();        
            void AddVoxel(vtkIdType i, unsigned int x, unsigned int y, unsigned int z, double w);
            void GetVoxel(vtkIdType i, _voxel *voxel);
            vtkIdType GetNumberOfVoxels();
            double GetAverageIntensity();

            _shell() {

            }
            ~_shell() {

            }
        
    };

    class _blob {

        private:

            std::vector<_shell> Shells;

        public:
        
            void Clear();
            void AddShell(_shell shell);
            void GetShell(vtkIdType i, _shell *shell);
            void GetMostRecentShell(_shell *shell);
            vtkIdType GetNumberOfShells();

            _blob() {

            }
            ~_blob() {

            }
        
    };

    class _expanding_blob {

        private:

            _blob Blob;
            int nneigh;
            vtkImageData *Volume;
            vtkUnsignedCharArray *State;
            std::vector<int> _dx, _dy, _dz;

        public:

            void PrintNeigh();        
            void SetSeed(unsigned int x, unsigned int y, unsigned int z);
            void SetNeighborhood(unsigned char neighsize);
            int Expand(bool decompose);
            void ExportCurrentShell(std::string fname);
            void SwapImageContent() {
                Volume -> GetPointData() -> GetScalars() -> DeepCopy(State);
            }

            _expanding_blob(unsigned int neightype, vtkImageData *Image) {
                Volume = Image;
                SetNeighborhood(neightype);
                State = vtkUnsignedCharArray::New();
                State -> SetNumberOfComponents(1);
                State -> SetNumberOfTuples(Volume->GetNumberOfPoints());
                State -> FillComponent(0,0);
            }
            ~_expanding_blob() {

            }
        
    };

    void FillHoles(vtkSmartPointer<vtkImageData> ImageData);

    int GetNumberOfComponents(vtkSmartPointer<vtkImageData> ImageData);

    class _cc_finder {

        private:

            int *Dim;
            int *Ext;
            vtkIdType N;
            vtkImageData *Vol;
            std::vector<vtkIdType> CSize;
            double _foreground, _background;

        public:

            vtkIdType Update();
            void SetBackground(double value);
            void SetForeground(double value);

            vtkIdType GetNumberOfComponents();
            vtkIdType GetComponentSize(unsigned short i);

            _cc_finder(vtkImageData *Image) {
                Vol = Image;
                Ext = Image -> GetExtent();
                Dim = Image -> GetDimensions();
                N = Image -> GetNumberOfPoints();
            }
            ~_cc_finder() {

            }
        
    };

    //=============

    struct _point {
        double x, y, z;
    };

    class _Supernova {
        
        private:    

            int _rmax, _nrays, _freq;
            double _xo, _yo, _zo, _scalefactor;
            vtkPolyData *Rays, *Peaks, *Cell, *CellOuter;

        public:

            int GetMaximumRadius() {
                return _rmax;
            }    
            void SetMaximumRadius(int rmax) {
                _rmax = rmax;
            }
            void SetCenter(double xo, double yo, double zo) {
                _xo = xo;
                _yo = yo;
                _zo = zo;
            }
            void SetScaleFactor(double _dxy, double _dz) {
                _scalefactor = _dz/_dxy;
            }
            int GetNumberOfRays() {
                return _nrays;
            }
            void SetNumberOfRays(int nrays) {
                _nrays = nrays;
            }
            void SetFrequency(int freq) {
                _freq = freq;
            }
            void SetPoints(vtkPoints *Points) {
                Rays -> SetPoints(Points);
            }
            void SetLines(vtkCellArray *Array) {
                Rays -> SetLines(Array);
            }
            vtkDataArray *GetIntensities() {
                return Rays -> GetPointData() -> GetScalars();
            } 
            void SetLinesIntensity(vtkUnsignedShortArray *Intensities) {
                Rays -> GetPointData() -> SetScalars(Intensities);
            }
            void Initialize();

            void Probe(vtkImageData *Image);

            void Save(const char FileName[], const char FileName2[]);

            void SaveRays(const char FileName[]);

            void ApplyLimits(const double r1, const bool force);

            void Segmentation();

            void ClipImageData(const char MitoFileName[], vtkImageData *Image, vtkImageData *ClipImage, const int _id);

            void ScalePolyData(double _dxy, double _dz);

            void AjustCoordinates(int Ly, double _dxy, double _dz);

            void EstimateImageMeanAndStd(vtkImageData *Image, double *_mean, double *_std);

            void GetXYZFromRay(const int ray, double *x, double *y, double *z);

            void SaveMassProperties(const char MassFileName[]);

            _Supernova() {
                _freq = 40;
                _rmax = 120;
                _nrays = 5000;
            }

            ~_Supernova();

    };

#endif