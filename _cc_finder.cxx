#include "_cc_finder.h"

    int ssdx_sort[26] = { 0,-1, 0, 1, 0, 0,-1, 0, 1, 0,-1, 1, 1,-1,-1, 0, 1, 0, -1, 1, 1,-1,-1, 1, 1,-1};
    int ssdy_sort[26] = { 0, 0,-1, 0, 1, 0, 0,-1, 0, 1,-1,-1, 1, 1, 0,-1, 0, 1, -1,-1, 1, 1,-1,-1, 1, 1};
    int ssdz_sort[26] = {-1, 0, 0, 0, 0, 1,-1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1, -1,-1,-1,-1, 1, 1, 1, 1};

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
