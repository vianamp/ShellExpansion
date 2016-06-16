#include "_eblob.h"

    void _shell::Clear() {
        Voxels.clear();
    }

    void _shell::AddVoxel(vtkIdType i, unsigned int x, unsigned int y, unsigned int z, double w) {
        _voxel voxel = {i,x,y,z,w};
        Voxels.push_back(voxel);
    }
    void _shell::GetVoxel(vtkIdType i, _voxel *voxel) {
        voxel -> i = Voxels[i].i;
        voxel -> x = Voxels[i].x;
        voxel -> y = Voxels[i].y;
        voxel -> z = Voxels[i].z;
        voxel -> w = Voxels[i].w;
    }

    double _shell::GetAverageIntensity() {
        double avg = 0.0;
        for (int i = 0; i < Voxels.size(); i++) {
            avg += Voxels[i].w;
        }
        return avg/Voxels.size();
    }

    void _shell::Decompose(vtkIdType Smin) {
        vtkIdType i;
        int xmin = 1E4, ymin = 1E4, zmin = 1E4;
        int xmax = 0E4, ymax = 0E4, zmax = 0E4;
        for (i = 0; i < Voxels.size(); i++) {
            xmin = (Voxels[i].x < xmin) ? Voxels[i].x : xmin;
            ymin = (Voxels[i].y < ymin) ? Voxels[i].y : ymin;
            zmin = (Voxels[i].z < zmin) ? Voxels[i].z : zmin;
            xmax = (Voxels[i].x > xmax) ? Voxels[i].x : xmax;
            ymax = (Voxels[i].y > ymax) ? Voxels[i].y : ymax;
            zmax = (Voxels[i].z > zmax) ? Voxels[i].z : zmax;
        }

        vtkSmartPointer<vtkImageData> Vol = vtkSmartPointer<vtkImageData>::New();
        Vol -> SetExtent(xmin,xmax,ymin,ymax,zmin,zmax);
        Vol -> SetSpacing(1,1,1);

        vtkSmartPointer<vtkUnsignedCharArray> Scalars = vtkSmartPointer<vtkUnsignedCharArray>::New();
        Scalars -> SetNumberOfComponents(1);
        Scalars -> SetNumberOfTuples((xmax-xmin+1)*(ymax-ymin+1)*(zmax-zmin+1));
        Scalars -> FillComponent(0,0);

        for (i = 0; i < Voxels.size(); i++) {
            Scalars->SetTuple1(Vol->FindPoint(Voxels[i].x,Voxels[i].y,Voxels[i].z),255);
        }

        Vol -> GetPointData() -> SetScalars(Scalars);

        _cc_finder CCFinder(Vol);

        if (CCFinder.Update() < Smin) {
            vtkIdType csize;
            unsigned short cluster;
            for (std::vector<_voxel>::iterator it=Voxels.begin(); it!=Voxels.end();) {
                cluster = (unsigned short)Vol -> GetScalarComponentAsDouble((*it).x,(*it).y,(*it).z,0);
                csize = CCFinder.GetComponentSize(cluster);
                if (csize < Smin) {
                    it = Voxels.erase(it);
                } else {
                  ++it;
                }
            }
        }

    }

    vtkIdType _shell::GetNumberOfVoxels() {
        return (vtkIdType)Voxels.size();
    }

    void _blob::Clear() {
        Shells.clear();
    }

    void _blob::AddShell(_shell shell) {
        _voxel v;
        _shell new_shell;
        for (vtkIdType i = 0; i < shell.GetNumberOfVoxels(); i++) {
            shell.GetVoxel(i,&v);
            new_shell.AddVoxel(v.i,v.x,v.y,v.z,v.w);
        }
        Shells.push_back(new_shell);
    }

    void _blob::GetShell(vtkIdType i, _shell *shell) {
        _voxel v;
        for (vtkIdType p = 0; p < Shells[i].GetNumberOfVoxels(); p++) {
            Shells[i].GetVoxel(p,&v);
            shell -> AddVoxel(v.i,v.x,v.y,v.z,v.w);
        }
    }

    void _blob::GetMostRecentShell(_shell *shell) {
        _voxel v;
        for (vtkIdType p = 0; p < Shells.back().GetNumberOfVoxels(); p++) {
            Shells.back().GetVoxel(p,&v);
            shell -> AddVoxel(v.i,v.x,v.y,v.z,v.w);
        }        
    }

    vtkIdType _blob::GetNumberOfShells() {
        return (vtkIdType)Shells.size();
    }

    void _expanding_blob::PrintNeigh() {
        for (int i = 0; i < _dx.size(); i++) printf("%d ",_dx[i]); printf("\n");
        for (int i = 0; i < _dy.size(); i++) printf("%d ",_dy[i]); printf("\n");
        for (int i = 0; i < _dz.size(); i++) printf("%d ",_dz[i]); printf("\n");
    }

    void _expanding_blob::SetSeed(unsigned int x, unsigned int y, unsigned int z) {
        _shell Shell;
        vtkIdType id = Volume -> FindPoint(x,y,z);
        State -> SetTuple1(id,255);
        Shell.AddVoxel(id,x,y,z,Volume->GetScalarComponentAsDouble(x,y,z,0));
        Blob.AddShell(Shell);
    }

    void _expanding_blob::SetNeighborhood(unsigned char neightype) {
        int r = (int)sqrt(neightype);
        for (int x = -r; x <= r; x++) {
            for (int y = -r; y <= r; y++) {
                for (int z = -r; z <= r; z++) {
                    if ( (x*x+y*y+z*z > 0) & (x*x+y*y+z*z <= neightype) ) {
                        _dx.push_back(x);
                        _dy.push_back(y);
                        _dz.push_back(z);
                    }
                }
            }
        }
        nneigh = (int)_dx.size();
    }

    int _expanding_blob::Expand(double threshold, vtkIdType Smin) {
        _voxel v;
        double w, t;
        int k, x, y, z;
        vtkIdType i, id;
        _shell curr_shell, next_shell;
        Blob.GetMostRecentShell(&curr_shell);
        for (i = 0; i < curr_shell.GetNumberOfVoxels(); i++) {
            curr_shell.GetVoxel(i,&v);
            for (k = 0; k < nneigh; k++) {
                x = v.x + _dx[k];
                y = v.y + _dy[k];
                z = v.z + _dz[k];
                t = threshold * (1-0.2*abs(_dz[k]));
                id = Volume -> FindPoint(x,y,z);
                if (id >= 0) {
                    if (!State->GetTuple1(id)) {
                        w = Volume -> GetScalarComponentAsDouble(x,y,z,0);
                        if (w < t)
                            next_shell.AddVoxel(id,x,y,z,w);
                        State -> SetTuple1(id,255);
                    }
                }
            }
        }
        if (next_shell.GetNumberOfVoxels()) {
            if (Smin) next_shell.Decompose(Smin);
            Blob.AddShell(next_shell);
        } else {
            return 0;
        }
        return 1;
    }

    void _expanding_blob::_expanding_blob::ExportCurrentShell(std::string fname) {
        _voxel v;
        _shell Shell;
        Blob.GetMostRecentShell(&Shell);
        
        vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();
        for (vtkIdType i = 0; i < Shell.GetNumberOfVoxels(); i++) {
            Shell.GetVoxel(i,&v);
            Points -> InsertNextPoint(v.x,v.y,v.z);
        }

        vtkSmartPointer<vtkPolyData> PolyData = vtkSmartPointer<vtkPolyData>::New();
        PolyData -> SetPoints(Points);

        vtkSmartPointer<vtkSphereSource> Sphere = vtkSmartPointer<vtkSphereSource>::New();
        Sphere -> SetThetaResolution(3);
        Sphere -> SetPhiResolution(3);
        Sphere -> SetRadius(1);

        vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
        glyph3D -> SetSourceConnection(Sphere -> GetOutputPort());
        glyph3D -> SetInputData(PolyData);
        glyph3D -> Update();

        vtkSmartPointer<vtkPolyDataWriter> Writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        Writer -> SetInputData(glyph3D->GetOutput());
        Writer -> SetFileName(fname.c_str());
        Writer -> Write();
    }
