#ifndef _CC_FINDER_H
#define _CC_FINDER_H

#include "includes.h"

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

#endif