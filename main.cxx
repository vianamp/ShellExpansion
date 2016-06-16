/*

    Region Growing Segmentation
    ---------------------------

    File xx.info contains information about the z-stack and a list of seeds coordinates.
    For each seed, we grow a blob shell-by-shell as we monitore the number of connected
    components of the current shell. Only components larger than `Smin` pixels are allowed
    to grow. In addition, voxels of th current shell are allowed to recruit new voxels
    for the next shell only if their intensity is lower than `threshold`.

    Parameters
    ----------

        threshold: maximum internsity allowed for a voxel
                   to be incorporated in the growing blob.

        Smin:      Minimum number of voxels that a component of
                   a shell has to have in order to grow.


    Matheus Palhares Viana, IBM Research | Brazil, 06.16.2016

    Things to do and further improvements
    -------------------------------------

    - The algorithm prevents the growing region from leaking in many circustances, but
      this requires the leakage to occurs as a disconnected component. If this is not
      the case, the algorithm will fail. Something else has to be done in those cases.

    - Monitore the curvature of the growing region voxel-by-voxel and use it as another
      parameter to determine whether a given voxel can grow.

    - Resample the resulting surface according to `SpacingXY` and `SpacingZ` specified
      info file.

    - Implement a new method for generating seeds automatically.

*/

#include "includes.h"

int SaveImageData(const char FileName[], vtkImageData *Image) {
    vtkSmartPointer<vtkStructuredPointsWriter> W = vtkSmartPointer<vtkStructuredPointsWriter>::New();
    W -> SetInputData(Image);
    W -> SetFileName(FileName);
    W -> Write();
    return EXIT_SUCCESS;
}

int ExpandBlobsFromSeed(_database *DataBase, double threshold, vtkIdType Smin) {

    #ifdef DEBUG
        printf("Shell Expansion [DEBUG mode]\n");
        printf("File name: %s\n",DataBase->GetPrefix().c_str());
    #endif

    vtkSmartPointer<vtkTIFFReader> TIFFReader = vtkSmartPointer<vtkTIFFReader>::New();
    TIFFReader -> SetOrientationType(TIFF_ORIENTATION_READER);
    TIFFReader -> SetFileName((DataBase->GetRootFolder()+DataBase->GetPrefix()+".tif").c_str());
    TIFFReader -> Update();

    vtkSmartPointer<vtkImageFlip> Flip = vtkSmartPointer<vtkImageFlip>::New();
    Flip -> SetFilteredAxis(1);
    Flip -> SetInputData(TIFFReader->GetOutput());
    Flip -> Update();

    vtkSmartPointer<vtkImageData> ImageData = Flip -> GetOutput();
    
    vtkSmartPointer<vtkImageData> BackupImage = vtkSmartPointer<vtkImageData>::New();
    BackupImage -> DeepCopy(ImageData);
    BackupImage -> Modified();

    int i, j, *Dim = BackupImage -> GetDimensions();

    #ifdef DEBUG
        printf("Volume dimensions: %dx%dx%d\n",Dim[0],Dim[1],Dim[2]);
    #endif

    int xo, yo, zo, tp, run, MAX_IT = 150;
    for (int id = 0; id < DataBase->GetNumberOfCenters(); id++) {

        xo = DataBase -> GetX(id);
        yo = DataBase -> GetY(id);
        zo = DataBase -> GetZ(id);

        printf("Running Seed [%d]: %d, %d, %d\n",DataBase->GetId(id),xo,yo,zo);

        _expanding_blob EBlob(1,ImageData);
        EBlob.SetSeed(xo,Dim[1]-yo,zo);

        run = 0;
        while(EBlob.Expand(threshold,(run>15)?Smin:0) && run < MAX_IT) {
            //EBlob.ExportCurrentShell((DataBase->GetRootFolder()+DataBase->GetPrefix()+"-"+std::to_string(DataBase->GetId(id))+"-points-"+std::to_string(run)+".vtk").c_str());
            run++;
        }

        (run<MAX_IT-1) ? printf("\tDone.\n") : printf("\tDone (MAX_IT reached).\n");

        EBlob.SwapImageContent();

        vtkSmartPointer<vtkImageGaussianSmooth> Gauss = vtkSmartPointer<vtkImageGaussianSmooth>::New();
        Gauss -> SetInputData(ImageData);
        Gauss -> SetStandardDeviations(3,3,2);
        Gauss -> SetRadiusFactors(3,3,2);
        Gauss -> Update();

        ImageData -> DeepCopy(Gauss->GetOutput());

        for (i = 0; i < Dim[0]; i++) {
            for (j = 0; j < Dim[1]; j++) {
                ImageData -> SetScalarComponentFromDouble(i,j,0,0,0);
                ImageData -> SetScalarComponentFromDouble(i,j,Dim[2]-1,0,0);
            }
        }

        ImageData -> GetPointData() -> GetScalars() -> Modified();

        vtkSmartPointer<vtkContourFilter> ContourFilter = vtkSmartPointer<vtkContourFilter>::New();
        ContourFilter -> SetInputData(ImageData);
        ContourFilter -> SetValue(0,120);
        ContourFilter -> Update();

        vtkSmartPointer<vtkPolyDataWriter> Writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        Writer -> SetInputData(ContourFilter->GetOutput());
        Writer -> SetFileName(DataBase->MakeGenericFileName(DataBase->GetId(id),"-cellsurface",".vtk").c_str());
        Writer -> Update();

        ImageData -> DeepCopy(BackupImage);

    }   

    return EXIT_SUCCESS;

}

int ScanFolderForThisExtension(const char _root[], const char ext[], std::vector<std::string> *List) {
    DIR *dir;
    int ext_p;
    struct dirent *ent;
    std::string _dir_name;
    if ((dir = opendir (_root)) != NULL) {
      while ((ent = readdir (dir)) != NULL) {
        _dir_name = std::string(ent->d_name);
        ext_p = (int)_dir_name.find(std::string(ext));
        if (ext_p > 0) {
            #ifdef DEBUG
                printf("File found: %s\n",_dir_name.c_str());
            #endif
            List -> push_back(std::string(_root)+_dir_name.substr(0,ext_p));
        }
      }
      closedir (dir);
    } else {
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {     
    
    std::string Path;
    vtkIdType Smin = 50;
    double threshold = 9000;
    _database *DataBase = new _database;
    for (int i = 0; i < argc; i++) {
        if (!strcmp(argv[i],"-path")) {
            Path = std::string(argv[i+1]);
            if (Path.back() != '/')
                Path = Path + '/';
        }
        if (!strcmp(argv[i],"-threshold")) {
            threshold = atof(argv[i+1]);
        }
        if (!strcmp(argv[i],"-smin")) {
            Smin = atoi(argv[i+1]);
        }
    }

    std::vector<std::string> Files;
    ScanFolderForThisExtension(Path.c_str(),".info",&Files);

    for (int i = 0; i < Files.size(); i++) {
        DataBase -> PopulateFromFile(Files[i]+".info");
        ExpandBlobsFromSeed(DataBase,threshold,Smin);
    }

    return EXIT_SUCCESS;
}
