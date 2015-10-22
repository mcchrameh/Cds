#ifndef _CAHNHILL_H
#define _CAHNHILL_H
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "VTK.h"
class CahnHill2D
{
  protected:
       int nx, ny;
       double delta_t, M, delta_x, delta_y,u,b,k, C1,C2;
       double tau,A,F,v,D,dxx,B,r;
       std::vector<std::vector<double> >phi;
       double **PHI,**PHI_old,**gamma,**Laplacian2;
       int *down1x, *down2x, *up1x, *up2x;
       //void setBC();
      // void computeBoundaryPoints();
       //void computeInnerPoints();
       void UpdateSolution();
       void initialCondition();
       void setLaplacianBase();
       void setIndex();
       double g(double phi);
       void SetSecondLaplacian2();
       void FiniteDifferenceScheme();
      // double gamma(int i, int j);


  public:
        CahnHill2D();
        CahnHill2D(int Nx,int Ny,double dt);
        virtual ~CahnHill2D();
       void  display();
       std::string make_output_filename(int index)
           {
              std::ostringstream ss;

              //ss << "output_" << index << ".dat";
              ss << "output_" << index << ".vtk";
              return ss.str();
           }

       void result(int count);
       void Solve();
       void Solver2();
       friend void save_vtk(double *U, int Nx, int Ny) ;


};

class Parallel_CahnHill2D: public  CahnHill2D
{
   protected:
            int  nlocalx,nlocaly, remainder, N, Nx,Ny;
            int  rank, size, tag,low;
            int  offset,startrow,endrow;
            double **PHI_p,**PHI_old_p,**gamma_p,**Laplacian2_p,*u_lastrow,*u_firstrow, *Boundary_bottom, *Boundary_top;
            double *Boundary_top1, *Boundary_bottom1, *array_gamma_top,*array_gamma_bottom, *u_store,*bc;
            int  *right1x,*left1x;

   public:

         Parallel_CahnHill2D();
         Parallel_CahnHill2D(int N1x, int N1y, int numProc, double dt);
        ~Parallel_CahnHill2D();
         void  Initialize_parallel();
         void  ExchangeData();
         void  ExchangeData2();
         void  ExchangeGammaData();
         void  ComputeSecondLaplacian();
         void  ComputeLaplacianBase_p();
         void  setIndex_p();
         void  setSecond_laplacian();
         void  FiniteDifferenceScheme_p();
         void  UpdateSolution_p();
         void  WriteToFile_MPI();
         void  SendToMaster(int count); //this  increases communication thus less efficient. 
         void  parallel_solver();
         void  WriteToFile();
         void    ReadFile(std::ifstream& infile);




};


