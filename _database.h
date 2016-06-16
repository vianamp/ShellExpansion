#ifndef _DATABASE_H
#define _DATABASE_H

#include "includes.h"

class _database {
    
    private:

        bool _checkmode;
        double _dxy, _dz;
        std::string _Prefix;
        std::string _RootFolder;
        std::vector<int> _id, _x, _y, _z;

        void Clear() {
              _id.clear();
               _x.clear();
               _y.clear();
               _z.clear();
        }

        void AddFromVector(std::vector<double> Vector) {
              _id.push_back((int)Vector[0]);
               _x.push_back((int)Vector[1]);
               _y.push_back((int)Vector[2]);
               _z.push_back((int)Vector[3]);
        }

    public:
         _database() {

         }
         ~_database();

           int GetId(int i) { return   _id[i]; }
            int GetX(int i) { return    _x[i]; }
            int GetY(int i) { return    _y[i]; }
            int GetZ(int i) { return    _z[i]; }

        double GetDxy() { return _dxy; }
        double  GetDz() { return  _dz; }

        void Print();
        
        std::string GetPrefix();

        std::string GetRootFolder();

        std::string MakeGenericFileName(int i, std::string Name, std::string Ext);

        int GetNumberOfCenters();

        void PopulateFromFile(const std::string CentersFileName);

};

#endif