
#include "_database.h"
#include "_shell.h"

int SaveImageData(const char FileName[], vtkImageData *Image) {
    vtkSmartPointer<vtkStructuredPointsWriter> W = vtkSmartPointer<vtkStructuredPointsWriter>::New();
    W -> SetInputData(Image);
    W -> SetFileName(FileName);
    W -> Write();
    return EXIT_SUCCESS;
}

int ExpandBlobsFromSeed(_database *DataBase) {

    #ifdef DEBUG
        printf("Shell Expansion [DEBUG mode]\n");
        printf("File name: %s\n",DataBase->GetFullCellName().c_str());
    #endif

    vtkSmartPointer<vtkTIFFReader> TIFFReader = vtkSmartPointer<vtkTIFFReader>::New();
    TIFFReader -> SetOrientationType(TIFF_ORIENTATION_READER);
    TIFFReader -> SetFileName(DataBase->GetFullCellName().c_str());
    TIFFReader -> Update();

    vtkSmartPointer<vtkImageData> ImageData = TIFFReader -> GetOutput();
    ImageData -> Modified();

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
        EBlob.SetSeed(xo,yo,zo);

        run = 0;
        while(EBlob.Expand((run>15)?1:0) && run < MAX_IT) {
            EBlob.ExportCurrentShell(std::string("temp"+std::to_string(run)+".vtk"));
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

/* ================================================================
   MAIN ROUTINE
=================================================================*/

int main(int argc, char *argv[]) {     
    
    std::string Path;
    _database *DataBase = new _database;
    for (int i = 0; i < argc; i++) {
        if (!strcmp(argv[i],"-path")) {
            Path = std::string(argv[i+1]);
            if (Path.back() != '/')
                Path = Path + '/';
        }
        if (!strcmp(argv[i],"-check")) {
            DataBase->SetCheckModeOn();
        }
    }

    std::vector<std::string> Files;
    ScanFolderForThisExtension(Path.c_str(),".centers",&Files);

    for (int i = 0; i < Files.size(); i++) {
        DataBase -> PopulateFromFile(Files[i]+".centers");
        ExpandBlobsFromSeed(DataBase);
        if ( DataBase->CheckMode() ) {
            #ifdef DEBUG
                printf("Running check mode...\n");
            #endif
            break;
        }
    }

    return EXIT_SUCCESS;
}
