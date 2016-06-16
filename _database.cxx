#include "_database.h"

    void _database::Print() {
        for (int i = 0; i < _x.size(); i++)
            printf("%d\t%d\t%d\t%d\n",_id[i],_x[i],_y[i],_z[i]);
    }
    
    std::string _database::GetRootFolder() {
        return _RootFolder;
    }

    std::string _database::GetPrefix() {
        return _Prefix;
    }

    std::string _database::MakeGenericFileName(int i, std::string Name, std::string Ext) {
        return _RootFolder + _Prefix + "-" + std::to_string(i) + Name + Ext;
    }

    int _database::GetNumberOfCenters() {
        return (int)_x.size();
    }

    void _database::PopulateFromFile(const std::string CentersFileName) {
        this -> Clear();
        std::string line;
        std::ifstream infile(CentersFileName.c_str());
        #ifdef DEBUG
            printf("Reading .info file...\n");
        #endif
        for( std::string line; getline( infile, line ); ) {
            if ( line == "[RootFolder]" ) {
                getline(infile,line);
                this -> _RootFolder = line;
                printf("> %s\n",this -> _RootFolder.c_str());
            }
            if ( line == "[Prefix]" ) {
                getline(infile,line);
                this -> _Prefix = line;
                printf("> %s\n",this -> _Prefix.c_str());
            }
            if ( line == "[SpacingXY]" ) {
                getline(infile,line);
                this -> _dxy = atof(line.c_str());
                printf("> dxy = %1.2f\n",this -> _dxy);
            }
            if ( line == "[SpacingZ]" ) {
                getline(infile,line);
                this -> _dz = atof(line.c_str());
                printf("> dz = %1.2f\n",this -> _dz);
            }
            if ( line == "[Centers]" ) {
                getline(infile,line);
                std::vector<double> Vector;
                while (line != "[end]") {
                    while (line.find(",") != std::string::npos) {
                        std::string value = line.substr(0,line.find(','));
                        line.erase(0,value.length()+1);
                        Vector.push_back(std::stof(value));
                    }
                    Vector.push_back(std::stof(line));
                    this -> AddFromVector(Vector);
                    getline(infile,line);
                    Vector.clear();
                }
                this -> Print();
            }
        }

    }
