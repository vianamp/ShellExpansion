#ifndef _DATABASE_H
#define _DATABASE_H

#include "includes.h"

class _database {
    
    private:

        bool _checkmode;
        double _dxy, _dz;
        std::string _Prefix;
        std::string _RootFolder;
        std::string _MitoFolder;
        std::string _CellFolder;
        std::vector<double> _r1, _r2;
        std::vector<int> _id, _x, _y, _z, _type;

        void Clear() {
              _id.clear();
            _type.clear();
               _x.clear();
               _y.clear();
               _z.clear();
              _r1.clear();
              _r2.clear();
        }

        void AddFromVector(std::vector<double> Vector) {
              _id.push_back((int)Vector[0]);
            _type.push_back((int)Vector[1]);
               _x.push_back((int)Vector[2]);
               _y.push_back((int)Vector[3]);
               _z.push_back((int)Vector[4]);
              _r1.push_back(Vector[5]);
              _r2.push_back(Vector[6]);
        }

    public:
         _database() {

         }
         ~_database();

           int GetId(int i) { return   _id[i]; }
         int GetType(int i) { return _type[i]; }
            int GetX(int i) { return    _x[i]; }
            int GetY(int i) { return    _y[i]; }
            int GetZ(int i) { return    _z[i]; }
        double GetR1(int i) { return   _r1[i]; }
        double GetR2(int i) { return   _r2[i]; }

        double GetDxy() { return _dxy; }
        double  GetDz() { return  _dz; }

        void Print();
        
        std::string GetFullCellName();

        std::string GetFullMitoName();

        std::string MakeGenericFileName(int i, std::string Name, std::string Ext);

        int GetNumberOfCenters();

        void PopulateFromFile(const std::string CentersFileName);

        void SetCheckModeOn() { _checkmode = 1; }

        bool CheckMode() { return _checkmode; }

};

#endif