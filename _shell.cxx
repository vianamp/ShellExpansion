#include "_shell.h"

    int ssdx_sort[26] = { 0,-1, 0, 1, 0, 0,-1, 0, 1, 0,-1, 1, 1,-1,-1, 0, 1, 0, -1, 1, 1,-1,-1, 1, 1,-1};
    int ssdy_sort[26] = { 0, 0,-1, 0, 1, 0, 0,-1, 0, 1,-1,-1, 1, 1, 0,-1, 0, 1, -1,-1, 1, 1,-1,-1, 1, 1};
    int ssdz_sort[26] = {-1, 0, 0, 0, 0, 1,-1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1, -1,-1,-1,-1, 1, 1, 1, 1};
    double _PI = 3.141592;


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

    void _shell::Decompose() {
        vtkIdType i;
        vtkIdType smin = 50;
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
        //printf("%d, %d, %d, %d, %d, %d\n",xmin,xmax,ymin,ymax,zmin,zmax);

        vtkSmartPointer<vtkImageData> Vol = vtkSmartPointer<vtkImageData>::New();
        Vol -> SetExtent(xmin,xmax,ymin,ymax,zmin,zmax);
        //Vol -> SetDimensions(xmax-xmin+1,ymax-ymin+1,zmax-zmin+1);
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

        if (CCFinder.Update() < smin) {
            vtkIdType csize;
            unsigned short cluster;
            for (std::vector<_voxel>::iterator it=Voxels.begin(); it!=Voxels.end();) {
                cluster = (unsigned short)Vol -> GetScalarComponentAsDouble((*it).x,(*it).y,(*it).z,0);
                csize = CCFinder.GetComponentSize(cluster);
                if (csize < smin) {
                    it = Voxels.erase(it);
                    //printf("erased...\n");
                } else {
                  ++it;
                }
            }
        }

        // vtkSmartPointer<vtkStructuredPointsWriter> Writer = vtkSmartPointer<vtkStructuredPointsWriter>::New();
        // Writer -> SetInputData(Vol);
        // Writer -> SetFileName("Seg.vtk");
        // Writer -> Update();


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

    int _expanding_blob::Expand(bool decompose) {
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
                t = 10000 * (1-0.2*abs(_dz[k]));
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
            if (decompose) next_shell.Decompose();
            Blob.AddShell(next_shell);
            //printf("S = %d	%1.8f\n",(int)next_shell.GetNumberOfVoxels(),next_shell.GetAverageIntensity());
        } else {
            return 0;
        }
        return 1;
    }

    void _expanding_blob::_expanding_blob::ExportCurrentShell(std::string fname) {
        vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();
        _voxel v;
        _shell Shell;
        Blob.GetMostRecentShell(&Shell);
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

    vtkIdType _cc_finder::Update() {
        
        double r[3];
        int x, y, z;
        bool find = 1;
        vtkIdType s, id, ido = N, scluster, label, smin = N;

        bool *State = new bool[N];
        vtkSmartPointer<vtkIdList> Curr = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> Next = vtkSmartPointer<vtkIdList>::New();

        vtkDataArray *Scalars = Vol -> GetPointData() -> GetScalars();

        for (id = N; id--;) State[id] = Scalars->GetTuple1(id) ? 1 : 0;

        label = 0;
        while (find) {
            for (s = 0; s < Curr->GetNumberOfIds(); s++) {
                Vol -> GetPoint(Curr->GetId(s),r);
                for (char i = 0; i < 26; i++) {
                    id = Vol -> FindPoint((int)r[0]+ssdx_sort[i],(int)r[1]+ssdy_sort[i],(int)r[2]+ssdz_sort[i]);
                    if (id >= 0) {
                        if (State[id]==1) {
                            Next -> InsertNextId(id);
                            State[id] = 0;
                            Scalars -> SetTuple1(id,label);
                            scluster++;
                        }
                    }
                }
            }
            if ( Next->GetNumberOfIds() ) {

                Curr -> Reset();
                Curr -> DeepCopy(Next);
                Next -> Reset();

            } else {

                if (label) {
                    CSize.push_back(scluster);
                    smin = (scluster<smin) ? scluster : smin;
                }

                find = 0;
                for (id=ido; id--;) {
                    if (State[id]) {
                        find = 1;
                        ido = id;
                        break;
                    }
                }
                        
                if (find) {
                    label++;
                    scluster = 1;
                    State[id] = 0;
                    Scalars -> SetTuple1(id,label);
                    Curr -> InsertNextId(id);
                }

            }
        }

        //printf("NCC = %d, S = %d\n",(int)CSize.size(),(int)CSize[0]);

        return smin;
    }
    void _cc_finder::SetBackground(double value) {
        _background = value;
    }
    void _cc_finder::SetForeground(double value) {
        _foreground = value;
    }

    vtkIdType _cc_finder::GetNumberOfComponents() {
        return (vtkIdType)CSize.size();
    }
    vtkIdType _cc_finder::GetComponentSize(unsigned short i) {
        return (vtkIdType)CSize[i-1];
    }


    //==========

    double FAux(double a, double t) {
        return -(t/_PI)*exp(-a * (t+_PI));
    }

    void _Supernova::GetXYZFromRay(const int ray, double *x, double *y, double *z) {
        double t = 2.0 * _PI*(1.0*ray)/_nrays - _PI;
        double alpha = 0.05;
        *z = (t<=0) ? FAux(alpha,t) : -FAux(alpha,-t);
        *x = cos(_freq*t);
        *y = sin(_freq*t);
    }

    void _Supernova::Initialize() {

        int m, p;
        double i, j, k;
        double t, x, y, z, n;
        vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> RayArray = vtkSmartPointer<vtkCellArray>::New();

        vtkIdType id = 0;
        for (m = 0; m < _nrays; m++) {

            GetXYZFromRay(m,&x,&y,&z);

            n = sqrt(x*x+y*y+z*z);

            RayArray -> InsertNextCell(_rmax);
            for (p = 0; p < _rmax; p++) {
                i = _xo + p * (x / n);
                j = _yo + p * (y / n);
                k = _zo + p * (z / n);
                Points -> InsertNextPoint(i,j,k);
                RayArray -> InsertCellPoint(id);
                id++;
            }

        }

        Rays = vtkPolyData::New();
        Rays -> SetPoints(Points);
        Rays -> SetLines(RayArray);
        
    }

    void _Supernova::Probe(vtkImageData *Image) {
        double v;
        int i, j, k;
        vtkIdType id;
        int *Dim = Image -> GetDimensions();

        vtkSmartPointer<vtkCellArray> RayArray = vtkSmartPointer<vtkCellArray>::New();

        vtkSmartPointer<vtkUnsignedShortArray> Intensities = vtkSmartPointer<vtkUnsignedShortArray>::New();
        Intensities -> SetNumberOfComponents(1);
        Intensities -> SetNumberOfTuples(_nrays*_rmax);
        Intensities -> FillComponent(0,0);

        double r[3];
        for (id = 0; id < Rays->GetNumberOfPoints(); id++) {
            v = 0;
            Rays -> GetPoint(id,r);
            i = (int)r[0]; j = (int)r[1]; k = (int)r[2];
            if ( (i >= 0 && i < Dim[0]) && (j >= 0 && j < Dim[1]) && (k >= 0 && k < Dim[2]) ) {
                v = Image -> GetScalarComponentAsDouble(i,j,k,0);
            }
            Intensities -> SetTuple1(id,v);
        }

        Rays -> GetPointData() -> SetScalars(Intensities);

        #ifdef DEBUG
            vtkSmartPointer<vtkPolyDataWriter> Writer =  vtkSmartPointer<vtkPolyDataWriter>::New();
            Writer -> SetFileTypeToBinary();
            Writer -> SetInputData(Rays);
            Writer -> SetFileName("Rays.vtk");
            Writer -> Write();
            printf("Rays saved in Supernove folder under the name Rays.vtk!\n");
        #endif

    }

    void _Supernova::ApplyLimits(const double r1, bool force) {
        vtkIdType r, i, id;
        double v, x, y, z, w, fr, gr;
        vtkSmartPointer<vtkDataArray> Scalars = Rays -> GetPointData() -> GetScalars();
        for (r = 0; r < _nrays; r++) {
            GetXYZFromRay(r,&x,&y,&z);
            fr = pow(fabs(z/sqrt(x*x+y*y+z*z)),2.0);
            w = r1 / ( fr * (_scalefactor-1) + 1 );
            for (i = 0; i < _rmax; i++) {
                id = i + r * _rmax;
                if (force) {
                    if (i == (int)w) Scalars -> SetTuple1(id,65535);
                } else {
                    v = Scalars -> GetTuple1(id);
                    v = (i<w)?v:0.0;
                    Scalars -> SetTuple1(id,v);
                }
            }
        }
    }

    void _Supernova::Save(const char FileName[], const char FileName2[]) {

        vtkSmartPointer<vtkPolyDataWriter> SurfaceWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
        SurfaceWriter -> SetFileTypeToBinary();
        SurfaceWriter -> SetInputData(Cell);
        SurfaceWriter -> SetFileName(FileName);
        SurfaceWriter -> Write();

        SurfaceWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
        SurfaceWriter -> SetFileTypeToBinary();
        SurfaceWriter -> SetInputData(CellOuter);
        SurfaceWriter -> SetFileName(FileName2);
        SurfaceWriter -> Write();

    }

    void _Supernova::SaveRays(const char FileName[]) {
        vtkSmartPointer<vtkPolyDataWriter> WriterR =  vtkSmartPointer<vtkPolyDataWriter>::New();
        WriterR -> SetFileTypeToBinary();
        WriterR -> SetInputData(Rays);
        WriterR -> SetFileName(FileName);
        WriterR -> Write();
    }

    void _Supernova::Segmentation() {

       #ifdef DEBUG
            printf("Saving Max projection...\n");
        #endif

        vtkSmartPointer<vtkImageData> Proj = vtkSmartPointer<vtkImageData>::New();
        Proj -> SetDimensions(_rmax,_nrays,1);
        
        vtkSmartPointer<vtkDataArray> Intensities = Rays -> GetPointData() -> GetScalars();
        
        #ifdef DEBUG
            Proj -> GetPointData() -> SetScalars(Intensities);
            vtkSmartPointer<vtkTIFFWriter> TIFFWriter = vtkSmartPointer<vtkTIFFWriter>::New();
            TIFFWriter -> SetFileName("Projection.tif");
            TIFFWriter -> SetFileDimensionality(2);
            TIFFWriter -> SetInputData(Proj);
            TIFFWriter -> Write();
            printf("Projection saved in Supernove folder under the name Projection.tif!\n");
        #endif

        int i, j, m;
        std::vector<double> Prec;
        std::vector<double>::iterator it;

        vtkSmartPointer<vtkDoubleArray> Map = vtkSmartPointer<vtkDoubleArray>::New();
        Map -> DeepCopy(Intensities);

        for (i = 1; i < _nrays; i++) {
            for (j = 0; j < _rmax; j++) {
                for (m = -1; m <=1; m++) {
                    if ( (j+m >= 0) && (j+m < _rmax) ) {
                        Prec.push_back(Map->GetTuple1((j+m)+(i-1)*_rmax));
                    }
                }
                it = std::max_element(Prec.begin(),Prec.end());
                Map -> SetTuple1(j+i*_rmax,Intensities->GetTuple1(j+i*_rmax)+(*it));
                Prec.clear();
            }
        }

        Map -> Modified();

        double v, map_range[2];
        Map -> GetRange(map_range);
        printf("Range = %1.3f / %1.3f\n",map_range[0],map_range[1]);

        for (i = 0; i < Intensities->GetNumberOfTuples(); i++) {
            v = 0;
            if (Intensities->GetTuple1(i) > 0) {
                v = 65535*(Map->GetTuple1(i)-map_range[0])/(map_range[1]-map_range[0]);
            }
            Intensities -> SetTuple1(i,v);
            Map -> SetTuple1(i,v);
        }

        Intensities -> Modified();

        #ifdef DEBUG
            Proj -> GetPointData() -> SetScalars(Intensities);
            TIFFWriter = vtkSmartPointer<vtkTIFFWriter>::New();
            TIFFWriter -> SetFileName("Map.tif");
            TIFFWriter -> SetFileDimensionality(2);
            TIFFWriter -> SetInputData(Proj);
            TIFFWriter -> Write();
            printf("Map saved in Supernove folder under the name Map.tif!\n");
        #endif

        for (j = _rmax; j--;) {
            if (Intensities->GetTuple1(j+(_nrays-1)*_rmax) == 65535)
                break;
        }
    
        #ifdef DEBUG
            printf("\tJo = %d\n",j);
        #endif

        std::vector< std::pair<int,int> > Path;
        Path.push_back(std::make_pair(_nrays-1,j));
        
        #ifdef DEBUG
            printf("\t\t[1]\n");
        #endif

        int delta;
        double mapv, mapmax;
        Intensities -> FillComponent(0,0);

        #ifdef DEBUG
            printf("\t\t[2]\n");
        #endif

        for (i = _nrays-1; i--;) {
            mapmax = 0;
            for (m = -1; m <=1; m++) {
                if ( (j+m >= 0) && (j+m < _rmax) ) {
                    mapv = Map -> GetTuple1((j+m)+(i-1)*_rmax);
                    if (mapv > mapmax) {
                        delta = m;
                        mapmax = mapv;
                    }
                }
            }
            j += delta;
            Path.push_back(std::make_pair(i,j));            
            Intensities -> SetTuple1(j+i*_rmax,65535);

        }

        #ifdef DEBUG
            printf("\tDone.\n");
        #endif


        #ifdef DEBUG
            Proj -> GetPointData() -> SetScalars(Intensities);
            TIFFWriter = vtkSmartPointer<vtkTIFFWriter>::New();
            TIFFWriter -> SetFileName("Path.tif");
            TIFFWriter -> SetFileDimensionality(2);
            TIFFWriter -> SetInputData(Proj);
            TIFFWriter -> Write();
            printf("Path saved in Supernove folder under the name Path.tif!\n");
        #endif

        vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();

        #ifdef DEBUG
            printf("\tScaling z-coordinate for surface reconstruction...\n");
        #endif

        int p, a, b;
        a = Path.front().second;
        b = Path.back().second;
        double t, x, y, z, n, zcoord, zmax, zmin;
        zmax = 0.0; zmin = 1E3;
        do {
            m = Path.back().first;
            p = Path.back().second;
            GetXYZFromRay(m,&x,&y,&z);            
            n = sqrt(x*x+y*y+z*z);
            zcoord = (_zo + p * (z / n)) * _scalefactor;
            zmax = (zcoord>zmax)?zcoord:zmax;
            zmin = (zcoord<zmin)?zcoord:zmin;
            Points -> InsertNextPoint(_xo + p * (x / n),_yo + p * (y / n),zcoord);
            Path.pop_back();
        } while (Path.size() > 0);

        #ifdef DEBUG
            printf("Z = [%1.3f,%1.3f]\n",zmin,zmax);
        #endif

        Points -> InsertNextPoint(_xo,_yo,zmin);
        Points -> InsertNextPoint(_xo,_yo,zmax);

        Peaks = vtkPolyData::New();
        Peaks -> SetPoints(Points);

        vtkSmartPointer<vtkPolyDataWriter> Writer;

        #ifdef DEBUG
            Writer =  vtkSmartPointer<vtkPolyDataWriter>::New();
            Writer -> SetFileTypeToBinary();
            Writer -> SetInputData(Peaks);
            Writer -> SetFileName("Peaks.vtk");
            Writer -> Write();
            printf("Peaks saved in Supernove folder under the name Peaks.vtk!\n");
        #endif

        vtkSmartPointer<vtkSurfaceReconstructionFilter> PoissonRecFilter = vtkSmartPointer<vtkSurfaceReconstructionFilter>::New();
        PoissonRecFilter -> SetNeighborhoodSize(100);
        PoissonRecFilter -> SetSampleSpacing(sqrt(10));
        PoissonRecFilter -> SetInputData(Peaks);
        PoissonRecFilter -> Update();

        vtkSmartPointer<vtkContourFilter> ContourFilter = vtkSmartPointer<vtkContourFilter>::New();
        ContourFilter -> SetInputConnection(PoissonRecFilter->GetOutputPort());
        ContourFilter -> ComputeNormalsOn();
        ContourFilter -> SetValue(0, 0.0);
        ContourFilter -> Update();

        Cell = vtkPolyData::New();
        Cell -> DeepCopy(ContourFilter->GetOutput());

        #ifdef DEBUG
            printf("\tScaling z-coordinate back...\n");
        #endif

        double r[3];
        vtkIdType id;
        for (id = Cell->GetNumberOfPoints(); id--;) {
            Cell -> GetPoint(id,r);
            Cell -> GetPoints() -> SetPoint(id,r[0],r[1],r[2]/_scalefactor);
        }
        Cell -> Modified();

        #ifdef DEBUG
            Writer =  vtkSmartPointer<vtkPolyDataWriter>::New();
            Writer -> SetFileTypeToBinary();
            Writer -> SetInputData(Cell);
            Writer -> SetFileName("Cell.vtk");
            Writer -> Write();
            printf("Cell saved in Supernove folder under the name Cell.vtk!\n");
        #endif

    }

    void FillHoles(vtkSmartPointer<vtkImageData> ImageData) {

        #ifdef DEBUG
            printf("\tSearching for holes in the image...\n");
        #endif

        int *Dim = ImageData -> GetDimensions();
        vtkIdType N = ImageData -> GetNumberOfPoints();

        vtkIdType i, s, ido, id;

        int x, y, z;
        double v, r[3];
        bool find = true;
        long long int ro[3];
        long int scluster, label;
        ro[0] = Dim[0] * Dim[1] * Dim[2];
        ro[1] = Dim[0] * Dim[1];

        vtkSmartPointer<vtkIdList> CurrA = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> NextA = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkLongArray> CSz = vtkSmartPointer<vtkLongArray>::New();
        vtkSmartPointer<vtkLongArray> Volume = vtkSmartPointer<vtkLongArray>::New();
        Volume -> SetNumberOfComponents(1);
        Volume -> SetNumberOfTuples(N);
        Volume -> FillComponent(0,0);

        for (x = 0; x < Dim[0]; x++) {
            for (y = 0; y < Dim[1]; y++) {
                for (z = 0; z < Dim[2]; z++) {
                    id = ImageData -> FindPoint(x,y,z);
                    if (id >= 0) {
                        if ((unsigned short int)ImageData->GetScalarComponentAsDouble(x,y,z,0)) {
                            Volume -> SetTuple1(id,1);
                        } else {
                            Volume -> SetTuple1(id,0);
                        }
                    }
                }
            }
        }

        Volume -> Modified();

        label = 0;
        while (find) {
            for (s = 0; s < CurrA->GetNumberOfIds(); s++) {
                ido = CurrA -> GetId(s);
                ImageData -> GetPoint(ido,r);
                x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
                for (i = 0; i < 6; i++) {
                    id = ImageData -> FindPoint(x+ssdx_sort[i],y+ssdy_sort[i],z+ssdz_sort[i]);
                    if (id >= 0) {
                        v = Volume -> GetTuple1(id);
                        if ((long int)v == 0) {
                            NextA -> InsertNextId(id);
                            Volume -> SetTuple1(id,-label);
                            scluster++;
                        }
                    }
                }
            }
            if (!NextA->GetNumberOfIds()) {
                find = false;
                for (id=ro[0]; id--;) {
                    v = Volume -> GetTuple1(id);
                    if ((long int)v == 0) {
                        find = true;
                        ro[0] = id;
                        break;
                    }
                }
                if (label) {
                    CSz -> InsertNextTuple1(scluster);
                }
                if (find) {
                    label++;
                    scluster = 1;
                    Volume -> SetTuple1(id,-label);
                    CurrA -> InsertNextId(id);
                }
            } else {
                CurrA -> Reset();
                CurrA -> DeepCopy(NextA);
                NextA -> Reset();
            }
        }

        int biggest;
        long int max_size = 0;
        for (id = 0; id < CSz->GetNumberOfTuples(); id++) {
            if (CSz->GetTuple1(id) > max_size) {
                biggest = id;
                max_size = CSz->GetTuple1(id);
            }
        }

        biggest  = -(biggest+1);

        printf("Biggest = %d\n",biggest);

        for (id = N; id--;) {
            if ((long int)Volume->GetTuple1(id) == biggest) {
                ImageData -> GetPointData() -> GetScalars() -> SetTuple1(id,0);
            } else if ((long int)Volume->GetTuple1(id) < biggest) {
                ImageData -> GetPointData() -> GetScalars() -> SetTuple1(id,255);
            }
        }
        ImageData -> GetPointData() -> GetScalars() -> Modified();

        #ifdef DEBUG
            printf("\tNumber of filled holes: %ld\n",(long int)CSz->GetNumberOfTuples()-1);
        #endif

    }


    int GetNumberOfComponents(vtkSmartPointer<vtkImageData> ImageData) {

        #ifdef DEBUG
            printf("\tLooking for Connected Components in the Image...\n");
        #endif

        int *Dim = ImageData -> GetDimensions();
        vtkIdType N = ImageData -> GetNumberOfPoints();

        vtkIdType i, s, ido, id;

        int x, y, z;
        double r[3];
        bool find = true;
        long long int ro[3];
        long int scluster, label;
        ro[0] = Dim[0] * Dim[1] * Dim[2];
        ro[1] = Dim[0] * Dim[1];

        vtkSmartPointer<vtkIdList> CurrA = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> NextA = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkLongArray> CSz = vtkSmartPointer<vtkLongArray>::New();
        vtkSmartPointer<vtkUnsignedCharArray> State = vtkSmartPointer<vtkUnsignedCharArray>::New();
        State -> SetNumberOfComponents(1);
        State -> SetNumberOfTuples(N);
        State -> FillComponent(0,0);

        for (id = N; id--;) State -> SetTuple1(id,ImageData->GetPointData()->GetScalars()->GetTuple1(id)?1:0);

        label = 0;
        while (find) {
            for (s = 0; s < CurrA->GetNumberOfIds(); s++) {
                ido = CurrA -> GetId(s);
                ImageData -> GetPoint(ido,r);
                x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
                for (i = 0; i < 6; i++) {
                    id = ImageData -> FindPoint(x+ssdx_sort[i],y+ssdy_sort[i],z+ssdz_sort[i]);
                    if (id >= 0) {
                        if (State->GetTuple1(id)) {
                            NextA -> InsertNextId(id);
                            State -> SetTuple1(id,0);
                            scluster++;
                        }
                    }
                }
            }
            if (!NextA->GetNumberOfIds()) {
                find = false;
                for (id=ro[0]; id--;) {
                    if (State->GetTuple1(id)) {
                        find = true;
                        ro[0] = id;
                        break;
                    }
                }
                
                if (label) CSz -> InsertNextTuple1(scluster);

                if (find) {
                    label++;
                    scluster = 1;
                    State -> SetTuple1(id,0);
                    CurrA -> InsertNextId(id);
                }
            } else {
                CurrA -> Reset();
                CurrA -> DeepCopy(NextA);
                NextA -> Reset();
            }
        }

        printf("NCC = %d, S = %d\n",(int)CSz->GetNumberOfTuples(),(int)CSz->GetTuple1(0));

        return CSz -> GetNumberOfTuples();

    }

    int LabelComponents(vtkSmartPointer<vtkImageData> ImageData) {

        #ifdef DEBUG
            printf("\tLooking for Connected Components in the Image...\n");
        #endif

        int *Dim = ImageData -> GetDimensions();
        vtkIdType N = ImageData -> GetNumberOfPoints();

        vtkIdType i, s, ido, id;

        int x, y, z;
        double r[3];
        bool find = true;
        long long int ro[3];
        long int scluster, label;
        ro[0] = Dim[0] * Dim[1] * Dim[2];
        ro[1] = Dim[0] * Dim[1];

        vtkSmartPointer<vtkIdList> CurrA = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> NextA = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkLongArray> CSz = vtkSmartPointer<vtkLongArray>::New();
        vtkSmartPointer<vtkFloatArray> State = vtkSmartPointer<vtkFloatArray>::New();
        State -> SetNumberOfComponents(1);
        State -> SetNumberOfTuples(N);
        State -> FillComponent(0,0);

        for (id = N; id--;) State -> SetTuple1(id,ImageData->GetPointData()->GetScalars()->GetTuple1(id)?1:0);

        label = 0;
        while (find) {
            for (s = 0; s < CurrA->GetNumberOfIds(); s++) {
                ido = CurrA -> GetId(s);
                ImageData -> GetPoint(ido,r);
                x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
                for (i = 0; i < 6; i++) {
                    id = ImageData -> FindPoint(x+ssdx_sort[i],y+ssdy_sort[i],z+ssdz_sort[i]);
                    if (id >= 0) {
                        if (State->GetTuple1(id)) {
                            NextA -> InsertNextId(id);
                            State -> SetTuple1(id,-label);
                            scluster++;
                        }
                    }
                }
            }
            if (!NextA->GetNumberOfIds()) {
                find = false;
                for (id=ro[0]; id--;) {
                    if (State->GetTuple1(id)) {
                        find = true;
                        ro[0] = id;
                        break;
                    }
                }
                
                if (label) CSz -> InsertNextTuple1(scluster);

                if (find) {
                    label++;
                    scluster = 1;
                    State -> SetTuple1(id,-label);
                    CurrA -> InsertNextId(id);
                }
            } else {
                CurrA -> Reset();
                CurrA -> DeepCopy(NextA);
                NextA -> Reset();
            }
        }

        printf("NCC = %d, S = %d\n",(int)CSz->GetNumberOfTuples(),(int)CSz->GetTuple1(0));

        return CSz -> GetNumberOfTuples();

    }

    void _Supernova::EstimateImageMeanAndStd(vtkImageData *Image, double *_mean, double *_std) {
        std::vector<unsigned short> v;
        for (vtkIdType id = Image->GetNumberOfPoints(); id--;) {
            v.push_back(Image->GetPointData()->GetScalars()->GetTuple1(id));
        }
        double sum = std::accumulate(v.begin(), v.end(), 0.0);
        double mean = sum / v.size();
        std::vector<double> diff(v.size());
        std::transform(v.begin(), v.end(), diff.begin(),std::bind2nd(std::minus<double>(), mean));
        double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
        double stdev = std::sqrt(sq_sum / v.size());
        *_mean = mean;
        *_std = stdev;
        printf("\t <I> = %1.3f\n",mean);
        printf("\t sd_I = %1.3f\n",stdev);
    }

    void _Supernova::ClipImageData(const char MitoFileName[], vtkImageData *Image, vtkImageData *ClipImage, const int _id) {

        #ifdef DEBUG
            printf("Clipping 3D Volume...\n");
        #endif

        vtkSmartPointer<vtkSelectEnclosedPoints> Filter = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
        Filter -> SetInputData(Image);
        Filter -> SetSurfaceData(Cell);
        Filter -> Update();

        double range[2];
        Filter -> GetOutput() -> GetScalarRange(range);

        printf("Range = [%1.3f,%1.3f]\n",range[0],range[1]);

        vtkSmartPointer<vtkImageShiftScale> Cast = vtkSmartPointer<vtkImageShiftScale>::New();
        Cast -> SetInputData(Filter->GetOutput());
        Cast -> SetShift(0.0);
        Cast -> SetScale(65535.0/range[1]);
        Cast -> SetOutputScalarTypeToUnsignedShort();
        Cast -> Update();

        vtkSmartPointer<vtkImageData> BinaryImage = Cast -> GetOutput();

        unsigned short *qPixel, *pPixel = static_cast<unsigned short*>(BinaryImage->GetScalarPointer());
        for( vtkIdType id = BinaryImage->GetNumberOfPoints(); id--; )
            if (pPixel[id]>0) pPixel[id] = 65535;

        vtkSmartPointer<vtkImageGaussianSmooth> Gauss = vtkSmartPointer<vtkImageGaussianSmooth>::New();
        Gauss -> SetInputData(BinaryImage);
        Gauss -> SetStandardDeviations(3,3,2);
        Gauss -> SetRadiusFactors(3,3,2);
        Gauss -> Update();

        vtkSmartPointer<vtkImageData> SmoothImage = vtkSmartPointer<vtkImageData>::New();
        SmoothImage -> DeepCopy(Gauss->GetOutput());

        #ifdef DEBUG
            printf("Capping smooth image...\n");
        #endif

        int i, j, *Dim = BinaryImage -> GetDimensions();
        for (i = 0; i < Dim[0]; i++) {
            for (j = 0; j < Dim[1]; j++) {
                SmoothImage -> SetScalarComponentFromDouble(i,j,0,0,0);
                SmoothImage -> SetScalarComponentFromDouble(i,j,Dim[2]-1,0,0);
            }
        }

        #ifdef DEBUG
            printf("Saving smooth image...\n");
            vtkSmartPointer<vtkStructuredPointsWriter> Writer = vtkSmartPointer<vtkStructuredPointsWriter>::New();
            Writer -> SetInputData(SmoothImage);
            Writer -> SetFileName("Smooth.vtk");
            Writer -> Update();
        #endif


        double threshold = 10000;//0.5*65535;
        vtkSmartPointer<vtkContourFilter> ContourFilter = vtkSmartPointer<vtkContourFilter>::New();
        ContourFilter -> SetInputData(SmoothImage);
        ContourFilter -> SetValue(0,threshold);
        ContourFilter -> Update();

        Cell -> DeepCopy(ContourFilter->GetOutput());

        ContourFilter -> SetInputData(SmoothImage);
        ContourFilter -> SetValue(0,0.05*65535);
        ContourFilter -> Update();

        CellOuter = vtkPolyData::New();
        CellOuter -> DeepCopy(ContourFilter->GetOutput());

        vtkSmartPointer<vtkIntArray> ID = vtkSmartPointer<vtkIntArray>::New();
        ID -> SetName("ID");
        ID -> SetNumberOfTuples(Cell->GetNumberOfPoints());
        ID -> SetNumberOfComponents(1);
        ID -> FillComponent(0,_id);

        Cell -> GetPointData() -> SetScalars(ID);
        Cell -> Modified();

        #ifdef DEBUG
            printf("Loading mitochondrial image...\n");
        #endif

        vtkSmartPointer<vtkTIFFReader> Reader = vtkSmartPointer<vtkTIFFReader>::New();
        Reader -> SetFileName(MitoFileName);
        Reader -> SetOrientationType(TIFF_ORIENTATION_READER);
        Reader -> Update();

        vtkSmartPointer<vtkImageData> Mito = Reader -> GetOutput();

        double mean, stdev;
        EstimateImageMeanAndStd(Mito,&mean,&stdev);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<> NormalDist(mean,0.5*stdev);

        int voi[6];
        double r[3];
        voi[0] = voi[2] = voi[4] = 1E3;
        voi[1] = voi[3] = voi[5] = 0;
        qPixel = static_cast<unsigned short*>(Mito->GetScalarPointer());
        pPixel = static_cast<unsigned short*>(SmoothImage->GetScalarPointer());
        for( vtkIdType id = SmoothImage->GetNumberOfPoints(); id--; ) {
            if (pPixel[id]>threshold) {
                BinaryImage -> GetPoint(id,r);
                voi[0] = ((int)r[0] < voi[0]) ? (int)r[0] : voi[0];
                voi[1] = ((int)r[0] > voi[1]) ? (int)r[0] : voi[1];
                voi[2] = ((int)r[1] < voi[2]) ? (int)r[1] : voi[2];
                voi[3] = ((int)r[1] > voi[3]) ? (int)r[1] : voi[3];
                voi[4] = ((int)r[2] < voi[4]) ? (int)r[2] : voi[4];
                voi[5] = ((int)r[2] > voi[5]) ? (int)r[2] : voi[5];
            } else {
                qPixel[id] = (unsigned short)NormalDist(gen);
            }
        }

        voi[0] -= (voi[0]-5>0) ? 5 : 0;
        voi[2] -= (voi[2]-5>0) ? 5 : 0;
        voi[1] += (voi[1]+5<Dim[0]) ? 5 : 0;
        voi[3] += (voi[3]+5<Dim[1]) ? 5 : 0;

        printf("VOI: x = [%d,%d], y = [%d,%d], z = [%d,%d]\n",voi[0],voi[1],voi[2],voi[3],0,Dim[2]-1);

        vtkSmartPointer<vtkExtractVOI> ExtractVOI = vtkSmartPointer<vtkExtractVOI>::New();
        ExtractVOI -> SetVOI(voi[0],voi[1],voi[2],voi[3],voi[4],voi[5]);
        ExtractVOI -> SetSampleRate(1,1,1);
        ExtractVOI -> SetInputData(Mito);
        ExtractVOI -> Update();

        ClipImage -> DeepCopy(ExtractVOI->GetOutput());

        #ifdef DEBUG
            printf("\t>Done.\n");
        #endif

    }

    void _Supernova::AjustCoordinates(int Ly, double _dxy, double _dz) {
        vtkSmartPointer<vtkCleanPolyData> Clean = vtkSmartPointer<vtkCleanPolyData>::New();
        Clean -> SetInputData(Cell);
        Clean -> Update();

        Cell -> DeepCopy(Clean -> GetOutput());

        vtkIdType i;
        double r[3], x, y, z;
        double xcm, ycm, zcm;
        xcm = ycm = zcm = 0.0;
        for (i = 0; i < Cell -> GetNumberOfPoints(); i++) {
            Cell -> GetPoint(i,r);
            xcm += r[0];
            ycm += r[1];
            zcm += r[2];
        }
        xcm /= Cell -> GetNumberOfPoints();
        ycm /= Cell -> GetNumberOfPoints();
        zcm /= Cell -> GetNumberOfPoints();

        printf("%1.3f\t%1.3f\t%1.3f\n",xcm,ycm,zcm);

        for (i = 0; i < Cell -> GetNumberOfPoints(); i++) {
            Cell -> GetPoint(i,r);
            x = 1.4 * _dxy * ( r[0] - xcm ) + _dxy * xcm;
            y = 1.4 * _dxy * ( r[1] - ycm ) + _dxy * ycm;
            z = 1.4 * _dz  * ( r[2] - zcm ) + _dz  * zcm;
            Cell -> GetPoints() -> SetPoint(i,x,_dxy*Ly-y,z);
        }
        Cell -> Modified();

        vtkSmartPointer<vtkFillHolesFilter> Holes = vtkSmartPointer<vtkFillHolesFilter>::New();
        Holes -> SetInputData(Cell);
        Holes -> SetHoleSize(100.0);
        Holes -> Update();

        vtkSmartPointer<vtkSmoothPolyDataFilter> Smooth = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
        if (Holes->GetOutput()->GetNumberOfCells() > 0) {
            Smooth -> SetInputData(Holes->GetOutput());
        } else {
            Smooth -> SetInputData(Cell);
        }
        Smooth -> SetNumberOfIterations(20);
        Smooth -> SetRelaxationFactor(0.1);
        Smooth -> Update();

        vtkSmartPointer<vtkPolyDataNormals> Normals = vtkSmartPointer<vtkPolyDataNormals>::New();
        Normals -> SetInputData(Smooth->GetOutput());
        Normals -> ComputeCellNormalsOn();
        Normals -> ComputePointNormalsOn();
        Normals -> AutoOrientNormalsOn();
        Normals -> Update();

        Cell -> DeepCopy(Normals->GetOutput());
        Cell -> Modified();
    }

    void _Supernova::ScalePolyData(double _dxy, double _dz) {
        double r[3];
        vtkSmartPointer<vtkPoints> P = Rays -> GetPoints();
        for (vtkIdType i = 0; i < P->GetNumberOfPoints(); i++) {
            P -> GetPoint(i,r);
            r[0] *= _dxy; r[1] *= _dxy; r[2] *= _dz;
            P -> SetPoint(i,r);
        }
        P = Peaks -> GetPoints();
        for (vtkIdType i = 0; i < P->GetNumberOfPoints(); i++) {
            P -> GetPoint(i,r);
            r[0] *= _dxy; r[1] *= _dxy; r[2] *= _dz;
            P -> SetPoint(i,r);
        }

        vtkSmartPointer<vtkTransform> Scale = vtkSmartPointer<vtkTransform>::New();
        Scale -> Scale(_dxy,_dxy,_dz);
        Scale -> Update();

        vtkSmartPointer<vtkTransformFilter> TFilter = vtkSmartPointer<vtkTransformFilter>::New();
        TFilter -> SetTransform(Scale);
        TFilter -> SetInputData(Cell);
        TFilter -> Update();

        Cell -> DeepCopy(TFilter->GetOutput());

        TFilter = vtkSmartPointer<vtkTransformFilter>::New();
        TFilter -> SetTransform(Scale);
        TFilter -> SetInputData(CellOuter);
        TFilter -> Update();

        CellOuter -> DeepCopy(TFilter->GetOutput());


    }

    void _Supernova::SaveMassProperties(const char MassFileName[]) {

        FILE *f = fopen(MassFileName,"w");
        fprintf(f,"volume\tsurface_area\tshape_factor\n");

        vtkSmartPointer<vtkMassProperties> Mass = vtkSmartPointer<vtkMassProperties>::New();
        Mass -> SetInputData(Cell);
        Mass -> Update();

        fprintf(f,"%1.3f\t%1.3f\t%1.3f\n",Mass->GetVolume(),Mass->GetSurfaceArea(),Mass->GetNormalizedShapeIndex());

        fclose(f);
    }
