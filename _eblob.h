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
            void Decompose(vtkIdType Smin);        
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
            int Expand(double threshold, vtkIdType Smin);
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

#endif