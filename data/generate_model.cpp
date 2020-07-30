#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main()
{
    int nx=2000;
    int nz=400;
    arma::Mat<float> model(nz,nx,fill::zeros);

    for(int ix=0; ix<nx; ix++)
    {
        for(int iz=0; iz<nz; iz++)
        {
            model(iz,ix)=2000;
        }
        for(int iz=nz/3; iz<nz;iz++)
        {
            model(iz,ix)=2500;
        }
        for(int iz=nz*2.0/3; iz<nz;iz++)
        {
            model(iz,ix)=3000;
        }
    }

    string fn="model_nz"+to_string(nz)+"_nx"+to_string(nx)+".dat";
    //string fn_nx=std::to_string(nx);
    model.save(fn,raw_binary);

    return 0;
}
