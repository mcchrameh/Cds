
#include "CahnH.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <time.h>
#include "mpi.h"
#define Random_min  -0.05 //-0.25 // for mim random number generator
#define Random_max  0.05 //0.25
#define PHI_old(i,j)     PHI_old[i][j]
#define gamma(i,j)       gamma[i][j]
#define Laplacian2(i,j)   Laplacian2[i][j]
#define gamma_p(i,j)     gamma_p[i][j]
#define PHI_old_p(i,j)   PHI_old_p[i][j]
#define Laplacian2_p(i,j) Laplacian2_p[i][j]
using namespace std;




CahnHill2D::CahnHill2D()
      {
        nx =0;
        ny=0;
        //phi=0; 
       
       }
CahnHill2D::CahnHill2D(int Nx, int Ny,double dt)
      {
          
          nx=Nx;ny=Ny;delta_x=1.0;delta_y=1.0; delta_t =0.001;M =1.0;b=M;k=M, u=0.5;
          C1=1.0/6.0; C2=1.0/12.0; //C1=1.0/6, C2=1.0/12
          dxx = delta_x;
          B=0.005; D=0.5;A=1.3;v=1.5;tau=0.36;F=0.5;r=0.5;
        for(int i=0;i<nx;i++) 
            {
              phi.push_back(vector<double>(ny));           }
         for(int i=0;i<nx;i++)
            {
             for(int j=0;j<ny;j++)
                 phi[i][j]=0.0; 
            }
          
          PHI=new double*[nx];
          PHI_old = new double*[nx];
          gamma = new double*[nx];
          Laplacian2 = new double*[nx];
         for(int i=0;i<nx; i++)
            {
              PHI[i]=new double[ny];
              PHI_old[i]=new double[ny]; 
              gamma[i]=new double[ny];
              Laplacian2[i]=new double [ny];
            }
          
          for(int i=0;i<nx;i++)
          {
              for(int j=0;j<ny;j++)
              {
                  PHI[i][j]=0.0;
                  PHI_old[i][j]=0.0;
                  gamma[i][j]=0.0;
                  Laplacian2[i][j]=0.0;
                
              }
          }
          
          
          //use these for getting the right indices
          down1x = new int [nx]; memset(down1x, 0, nx*sizeof(int));
          down2x = new int [nx]; memset(down2x, 0, nx*sizeof(int));
          up1x   = new int [nx]; memset(up1x, 0, nx*sizeof(int));
          up2x   = new int [nx]; memset(up2x, 0, nx*sizeof(int));
          cout<<"Constructor called"<<endl;
       }
CahnHill2D::~CahnHill2D()
      {

          for(int i=0;i<nx;i++)
             {
               delete [] PHI[i];
             }
           delete [] PHI;
 
            for(int i=0;i<nx;i++)
             {
               delete [] PHI_old[i];
             }
          delete [] PHI_old;
          
          for(int i=0;i<nx;i++)
          {
              delete [] gamma[i];
          }
          for(int i=0;i<nx;i++)
          {
              delete [] Laplacian2[i];
          }
          
          delete [] gamma;
          delete [] down1x;
          delete [] down2x;
          delete [] up1x;
          delete [] up2x;
          delete [] Laplacian2;
        cout <<"Destructor called"<<endl;

      }
void CahnHill2D:: setIndex()
    {
        for(int s=0; s<nx; s++)
           {
                up1x[s]=s+1;
                up2x[s]=s+2;
                down1x[s]=s-1;
                down2x[s]=s-2;
           }
        
        
        
            down1x[0]=nx-1;
            down2x[0]=nx-2;
            down2x[1]=nx-1;
            up1x[nx-1]=0;
            up2x[nx-1]=1;
            up2x[nx-2]=0;
        
        
        

    }

double CahnHill2D::g(double phi)
{
    double q=0.0;
   // q=(1.0 + tau - A*pow((1.0-2.0*F),2))*PHI_old(i,j)-v*(1.0-2.0*F)*pow(PHI_old(i,j),2)-u*pow(PHI_old(i,j),3);
    //q=(1.0 + tau - A*pow((1.0-2.0*F),2))*phi-v*(1.0-2.0*F)*pow(phi,2)-u*pow(phi,3);
    q = A*phi- (A/3.0)*pow(phi,3); // map from paper
   // cout<<"g(i,j)="<<q<<"\t";
    return q;
    
}
void CahnHill2D::setLaplacianBase()
     {
      //   cout<<"gamma(0,0)="<<gamma(0,0)<<"\t";
         double AP=0.0, BP=0.0,  ATP=0.0;//gtemp=0.0;
         setIndex();
         
         for(int i=0;i<nx;i++)
         {
             for(int j=0;j<ny;j++)
             {
                 
                 AP=C1*(PHI_old(up1x[i],j) + PHI_old(down1x[i],j)
                        + PHI_old(i,up1x[j]) + PHI_old(i,down1x[j]));
                 BP=C2*(PHI_old(down1x[i],up1x[j]) + PHI_old(down1x[i],down1x[j])
                        +PHI_old(up1x[i],up1x[j]) + PHI_old(up1x[i],down1x[j]));
                 ATP = AP + BP;
                // gtemp= (1.0 + tau - A*pow((1.0-2.0*F),2))*PHI_old(i,j)-v*(1.0-2.0*F)*pow(PHI_old(i,j),2)-u*pow(PHI_old(i,j),3);
                // cout<<"g(i,j)m="<<gtemp<<"\t";
                 //g_temp =g(i,j);
                 
                 gamma(i,j)=g(PHI_old(i,j))+ D*(ATP-PHI_old(i,j))-PHI_old(i,j);
        //         cout<<"gamma="<<gamma(i,j)<<"\t";
                 
                 
              }
             cout<<""<<endl;
         }
         
     }
void CahnHill2D::SetSecondLaplacian2()
{
    double AP=0.0, BP=0.0;
    setIndex();
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            AP=C1*(gamma(up1x[i],j) + gamma(down1x[i],j)
                   + gamma(i,up1x[j]) + gamma(i,down1x[j]));
            BP=C2*(gamma(down1x[i],up1x[j]) + gamma(down1x[i],down1x[j])
                   +gamma(up1x[i],up1x[j]) + gamma(up1x[i],down1x[j]));
            Laplacian2(i,j) = AP + BP;

            
            
        }
        
    }
    
    
}
void CahnHill2D::FiniteDifferenceScheme()
{
 //   double r=0.5;
    cout<<"in finite difference scheme"<<endl;
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
          // PHI[i][j]= PHI_old(i,j) + delta_t*(Laplacian2[i][j]-gamma(i,j) + B*(1.0)*PHI_old(i,j));
            PHI[i][j]= PHI_old(i,j) - B*(PHI_old(i,j)-1.0 + 2*r)+ gamma(i,j)-Laplacian2(i,j);
            
        }
        
    }
    
}

void CahnHill2D:: display()
        {

           for(int i=0; i<nx;i++)
              {
                for(int j=0;j<ny;j++)
                  {
                    cout<<""<<PHI_old[i][j]<<"\t";
                  }
                cout <<""<<endl;
              }


        }


void CahnHill2D::initialCondition()
     {
       
         for(int i=0; i<nx;i++)
            {
              for(int j=0;j<ny;j++)
                 {
                    double range =Random_max-Random_min;
                    double div =RAND_MAX / range;
                    PHI_old[i][j]=Random_min + (rand()/div);
                     
                 }

            }
      

     }


void CahnHill2D::result(int count)
{
    
        FILE *file = fopen(make_output_filename(count).c_str(), "w");
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            
            fprintf(file,"%15.11f \t",PHI[i][j]);
            //fprintf(file,"%15.11f \t",PHI_old[i][j]);
            
        }
        fprintf(file,"\n");
    }
    fclose(file);
    
    
}


    


void CahnHill2D::UpdateSolution()
{
    
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            PHI_old[i][j]=PHI[i][j];
        }
        
    }
    
    
}

void CahnHill2D::Solver2()
{

   initialCondition();
    int count =0;
    double t=0.0;
    while(t<1002)
    {
        std::cout<<"in while loop solver"<<std::endl;
        setLaplacianBase();
        SetSecondLaplacian2();
        FiniteDifferenceScheme();
        UpdateSolution();
        if (count==1000|| count==0)
            result( count);
       
        t+=1.0;
        count++;
        cout<<"time="<<t<<"\t";
        
    }

  
  
    cout<<""<<endl;
  
    
    
}

void CahnHill2D::Solve()
{
    initialCondition();
    
    //setBC();
    int count =0;
    double t=0.0;
    while(t<13)
    {
        //display();
       // computeBoundaryPoints();

      //  computeInnerPoints();
      //  computeBoundaryPoints();
        UpdateSolution();
       if (count==10 || count==0)
          result( count);
        //t+=delta_t;
        t+=1.0;
        count++;
     }
}


Parallel_CahnHill2D::Parallel_CahnHill2D(int N1x,int N1y, int numProc,double dt):CahnHill2D( N1x, N1y, dt)
{

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
 
  
    std::cout<<"size="<<size<<std::endl;
  
      //remainder=N1x % size;
     N =N1y;
     Nx=N1x;
     Ny = N1y;
    nlocalx = (int)Nx/size;

   startrow = rank*nlocalx;
   if(Nx%size!=0)
   {
    if(rank<size-1)
     {
      

    
       endrow = startrow + nlocalx-1;
     }

    else 
      {
      
        endrow=Nx-1;
        nlocalx=endrow-startrow + 1;
      }
    }
   else
     {
       nlocalx = (int)Nx/size;
     }
     nlocaly = nlocalx;

     tag=201;
   /*--------------------------------------------------------------------
    create local array (matrices) on each process
  ---------------------------------------------------------------------*/
    
       
     PHI_p=new double*[nlocalx ];
     PHI_old_p=new double*[nlocalx ];
     gamma_p = new double*[nlocalx ];
     Laplacian2_p= new double*[nlocalx];
    //for transfering boundary points for each process
    u_lastrow=new double[Ny];
    u_firstrow=new double[Ny];
    Boundary_top = new double[Ny];
    Boundary_bottom =new double[Ny];
    Boundary_top1 = new double[Ny];
    Boundary_bottom1=new double[Ny];
    u_store = new double[nlocalx*Ny]; // use by each process to store its solution and latter send it to the master for result viewing. 
    right1x = new int[Ny];
    left1x = new int [Ny];
    array_gamma_top = new double[N];
    array_gamma_bottom = new double[N];
    bc = new double[N];   
    std::memset(array_gamma_top, 0, N*sizeof(double));
    std::memset(array_gamma_bottom, 0, N*sizeof(double));
    std::memset(u_store, 0, nlocalx*Ny*sizeof(double));
    std::memset(bc, 0, N*sizeof(double));   
    std::memset(left1x, 0, Ny*sizeof(int));
    std::memset(u_firstrow, 0, Ny*sizeof(double));
    std::memset(u_lastrow, 0, Ny*sizeof(double));
    std::memset(Boundary_top, 0, Ny*sizeof(double));
    std::memset(Boundary_bottom, 0, Ny*sizeof(double));
    std::memset(Boundary_bottom1, 0, Ny*sizeof(double));
    std::memset(Boundary_top1, 0, Ny*sizeof(double));


         for(int i=0;i<nlocalx; i++)
            {
              PHI_p[i]=new double[Ny];
              PHI_old_p[i]=new double[Ny];
              gamma_p[i]=new double[Ny];
              Laplacian2_p[i]=new double [Ny];
            }

          for(int i=0;i<nlocalx  ;i++)
          { 
             for(int j=0;j<Ny;j++)
              {
                  PHI_p[i][j]=0.0;
                  PHI_old_p[i][j]=0.0;
                  gamma_p[i][j]=0.0;
                  Laplacian2_p[i][j]=0.0;

              }
          }
   
  std::cout<<"Parallel constructor called"<<std::endl;
}


Parallel_CahnHill2D::~Parallel_CahnHill2D()
{
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //MPI_Comm_size(MPI_COMM_WORLD, &size);

  for(int i =0; i<nlocalx; i++)
    {
      delete [] PHI_p[i];
      delete [] PHI_old_p[i];
      delete [] gamma_p[i];
      delete [] Laplacian2_p[i];
    }
     delete []PHI_p;
     delete []PHI_old_p;
     delete []gamma_p;
     delete []Laplacian2_p;
     delete []u_lastrow;
     delete []u_firstrow;
     delete []Boundary_top;
     delete []Boundary_bottom;
     delete []left1x;
     delete []right1x;
     delete []array_gamma_top;
     delete []array_gamma_bottom;
     delete []u_store;
     delete []bc;
    std::cout<<"destroyed parallel matrices"<<std::endl;

}


void Parallel_CahnHill2D::Initialize_parallel()
{

         MPI_Comm_rank(MPI_COMM_WORLD, &rank);
         MPI_Comm_size(MPI_COMM_WORLD, &size);
         srand(time(NULL) + rank);
         std::cout<<"I am rank="<<rank<<std::endl;
         for(int i=0; i<nlocalx ;i++)
            {
              for(int j=0;j<N;j++)
                 {
                    double range =Random_max-Random_min;
                    double div =RAND_MAX / range;
                    PHI_old_p[i][j]=Random_min + (rand()/div);
                  //  std::cout<<PHI_old_p[i][j]<<"\t";

                 }
                  // std::cout<<""<<std::endl;

            }
 std::cout<<" intialize data structure"<<std::endl;   
 
}

void Parallel_CahnHill2D:: setIndex_p()
    {   
        for(int s=0; s<N; s++)
           {    
                right1x[s]=s+1;
                
                left1x[s]=s-1;
                
           }


            
            left1x[0]=N-1;
            right1x[N-1]=0;
            
            



    
    }
void Parallel_CahnHill2D::ExchangeData2()
    { 
       MPI_Comm_rank(MPI_COMM_WORLD, &rank);
       MPI_Comm_size(MPI_COMM_WORLD, &size);
       tag = 202;
       MPI_Status status;
       //on each process copy all the data points to be sent 
       for(int k=0;k<N;k++)
          {
             u_lastrow[k]=PHI_old_p[nlocalx-1][k];
             u_firstrow[k]=PHI_old_p[0][k];
      
           }
            
    //This is for periodic boundary conditions used by process 0 and process size-1
    if (rank==0)
     {
         for(int j=0;j<N;j++)
              {
                 Boundary_top[j]= PHI_old_p[0][j];
                
              }
     }
     if (rank==size-1)
      {
         for(int j=0;j<N;j++)
              {
                Boundary_bottom[j]= PHI_old_p[nlocalx -1][j];
              }


      } 
    MPI_Barrier(MPI_COMM_WORLD);
    int  up_nbr =0, down_nbr=0; 
         up_nbr = rank + 1;
    if (up_nbr >= size) up_nbr = MPI_PROC_NULL;
    down_nbr = rank - 1;
    if (down_nbr < 0) down_nbr = MPI_PROC_NULL;
 //--------------exchange periodic boundary condition-----------------------------
    if (rank==0)
       {
          //  MPI_Send(u_firstrow,N,MPI_DOUBLE,size-1,tag,MPI_COMM_WORLD);
          //  MPI_Recv(u_lastrow,N,MPI_DOUBLE,size-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
          MPI_Sendrecv(Boundary_top,N,MPI_DOUBLE,size-1, tag,
                       Boundary_bottom, N, MPI_DOUBLE,size-1, tag,MPI_COMM_WORLD, &status);
          for(int j=0;j<N;j++)
             {
                 bc[j] = Boundary_bottom[j];
             }
                   

       }
    if (rank==size-1)
       {
            MPI_Sendrecv(Boundary_bottom,N,MPI_DOUBLE,0, tag,
                         Boundary_top, N, MPI_DOUBLE,0, tag,MPI_COMM_WORLD, &status);
             for(int j=0;j<N;j++)
             {
                 bc[j] = Boundary_top[j];
             }
 
       }
 //----------------end of exchange periodic boundary condition------------------------
    if((rank%2)==0)
      {
         MPI_Sendrecv(u_lastrow, N, MPI_DOUBLE, up_nbr, tag,
                      u_firstrow,N ,MPI_DOUBLE, up_nbr, tag, MPI_COMM_WORLD,  &status);

      }
    else
       {
        // std::cout<<"rank="<<rank<<std::endl;
          MPI_Sendrecv(u_lastrow, N, MPI_DOUBLE, up_nbr, tag,
                       u_firstrow,N ,MPI_DOUBLE, up_nbr, tag, MPI_COMM_WORLD,  &status);
       }
    if((rank%2)==1)
      {
    
          MPI_Sendrecv(u_firstrow, N, MPI_DOUBLE, down_nbr, tag,
                       u_lastrow,N ,MPI_DOUBLE, down_nbr, tag, MPI_COMM_WORLD,  &status);
      }
     else
      {
            MPI_Sendrecv(u_firstrow, N, MPI_DOUBLE, down_nbr, tag,
                         u_lastrow,N ,MPI_DOUBLE, down_nbr, tag, MPI_COMM_WORLD,  &status);

      }
 
}

void Parallel_CahnHill2D::ComputeLaplacianBase_p()
    {
       MPI_Comm_rank(MPI_COMM_WORLD, &rank);
       MPI_Comm_size(MPI_COMM_WORLD, &size);
       setIndex_p();

       if (rank==0)
          {
           // setIndex_p();
            double AP=0.0, BP=0.0, SUM_p0=0.0;
            for(int j=0;j<N;j++)
               {
                 const int i=0;
                 AP = C1*(PHI_old_p(i+1,j) + bc[j] + PHI_old_p(i, right1x[j]) + PHI_old_p(i,left1x[j]) );
                 BP = C2*(PHI_old_p(i+1,right1x[j]) + PHI_old_p(i+1,left1x[j]) + bc[left1x[j]] + bc[right1x[j]]);
                 SUM_p0 = AP + BP;
                 gamma_p(i,j)=CahnHill2D::g(PHI_old_p(i,j))+ D*(SUM_p0-PHI_old_p(i,j))-PHI_old_p(i,j);
               }
             AP =0.0;BP=0.0;SUM_p0=0.0;
            for(int i=1;i<nlocalx -1;i++)
               {
                for(int j=0;j<N;j++)
                  {
                     AP = C1*(PHI_old_p(i+1,j)+PHI_old_p(i-1,j) + PHI_old_p(i,right1x[j]) + PHI_old_p(i,left1x[j]));
                     BP = C2*(PHI_old_p(i+1,right1x[j]) + PHI_old_p(i+1,left1x[j]) + PHI_old_p(i-1,right1x[j]) + PHI_old_p(i-1,left1x[j]));
                     SUM_p0 = AP + BP;
                     gamma_p(i,j)=CahnHill2D::g(PHI_old_p(i,j))+ D*(SUM_p0-PHI_old_p(i,j))-PHI_old_p(i,j);
                   }
               }
             AP =0.0;BP=0.0;SUM_p0=0.0;
             for(int j=0;j<N;j++)
                 {
                   const int i =nlocalx-1;
                   AP = C1*(u_firstrow[j] + PHI_old_p(i-1,j) + PHI_old_p(i,right1x[j]) + PHI_old_p(i,left1x[j]));
                   BP = C2*(u_firstrow[right1x[j]] + u_firstrow[left1x[j]] + PHI_old_p(i-1,right1x[j]) + PHI_old_p(i-1,left1x[j]));
                   SUM_p0 = AP + BP;
                   gamma_p(i,j)=CahnHill2D::g(PHI_old_p(i,j))+ D*(SUM_p0-PHI_old_p(i,j))-PHI_old_p(i,j);

                 }
          }

      else if(rank==size-1)
        {
           // setIndex_p();
            double AP=0.0, BP=0.0, SUM_p0=0.0;
            for(int j=0;j<N;j++)
               {
                 const int i=0;
                 AP = C1*(PHI_old_p(i+1,j) + u_lastrow[j] + PHI_old_p(i, right1x[j]) + PHI_old_p(i,left1x[j]) );
                 BP = C2*(PHI_old_p(i+1,right1x[j]) + PHI_old_p(i+1,left1x[j]) + u_lastrow[left1x[j]] + u_lastrow[right1x[j]]);
                 SUM_p0 = AP + BP;
                 gamma_p(i,j)=CahnHill2D::g(PHI_old_p(i,j))+ D*(SUM_p0-PHI_old_p(i,j))-PHI_old_p(i,j);
               }
            AP =0.0;BP=0.0;SUM_p0=0.0;
            for(int i=1;i<nlocalx-1;i++)
               {
                for(int j=0;j<N;j++)
                  {
                     AP = C1*(PHI_old_p(i+1,j)+PHI_old_p(i-1,j) + PHI_old_p(i,right1x[j]) + PHI_old_p(i,left1x[j]));
                     BP = C2*(PHI_old_p(i+1,right1x[j]) + PHI_old_p(i+1,left1x[j]) + PHI_old_p(i-1,right1x[j]) + PHI_old_p(i-1,left1x[j]));
                     SUM_p0 = AP + BP;
                     gamma_p(i,j)=CahnHill2D::g(PHI_old_p(i,j))+ D*(SUM_p0-PHI_old_p(i,j))-PHI_old_p(i,j);
                   }
               }
             AP =0.0;BP=0.0;SUM_p0=0.0;
             for(int j=0;j<N;j++)
                 {
                     const int i =nlocalx-1;
                     AP = C1*(bc[j] + PHI_old_p(i-1,j) + PHI_old_p(i,right1x[j]) + PHI_old_p(i,left1x[j]) );
                     BP = C2*(bc[right1x[j]] + bc[left1x[j]] + PHI_old_p(i-1,right1x[j]) + PHI_old_p(i-1,left1x[j]));
                     SUM_p0 = AP + BP;
                     gamma_p(i,j)=CahnHill2D::g(PHI_old_p(i,j))+ D*(SUM_p0-PHI_old_p(i,j))-PHI_old_p(i,j);
                 }    

 
        }      
   else
       {
          double AP=0.0, BP=0.0, SUM_p0=0.0;
            for(int j=0;j<N;j++)
               {
                 const int i=0;
                 AP = C1*(PHI_old_p(i+1,j) + u_lastrow[j] + PHI_old_p(i, right1x[j]) + PHI_old_p(i,left1x[j]) );
                 BP = C2*(PHI_old_p(i+1,right1x[j]) + PHI_old_p(i+1,left1x[j]) + u_lastrow[left1x[j]] + u_lastrow[right1x[j]]);
                 SUM_p0 = AP + BP;
                 gamma_p(i,j)=CahnHill2D::g(PHI_old_p(i,j))+ D*(SUM_p0-PHI_old_p(i,j))-PHI_old_p(i,j);
               }
             AP =0.0;BP=0.0;SUM_p0=0.0;
              for(int i=1;i<nlocalx-1;i++)
               {
                for(int j=0;j<N;j++)
                  {
                     AP = C1*(PHI_old_p(i+1,j)+PHI_old_p(i-1,j) + PHI_old_p(i,right1x[j]) + PHI_old_p(i,left1x[j]));
                     BP = C2*(PHI_old_p(i+1,right1x[j]) + PHI_old_p(i+1,left1x[j]) + PHI_old_p(i-1,right1x[j]) + PHI_old_p(i-1,left1x[j]));
                     SUM_p0 = AP + BP;
                     gamma_p(i,j)=CahnHill2D::g(PHI_old_p(i,j))+ D*(SUM_p0-PHI_old_p(i,j))-PHI_old_p(i,j);
                   }
               }
             AP =0.0;BP=0.0;SUM_p0=0.0;
   
              for(int j=0;j<N;j++)
                 {
                   const int i =nlocalx -1;
                   AP = C1*(u_firstrow[j] + PHI_old_p(i-1,j) + PHI_old_p(i,right1x[j]) + PHI_old_p(i,left1x[j]));
                   BP = C2*(u_firstrow[right1x[j]] + u_firstrow[left1x[j]] + PHI_old_p(i-1,right1x[j]) + PHI_old_p(i-1,left1x[j]));
                   SUM_p0 = AP + BP;
                   gamma_p(i,j)=CahnHill2D::g(PHI_old_p(i,j))+ D*(SUM_p0-PHI_old_p(i,j))-PHI_old_p(i,j);

                 }


                                             


      }

    }


void Parallel_CahnHill2D::ExchangeGammaData()
    {
       MPI_Comm_rank(MPI_COMM_WORLD, &rank);
       MPI_Comm_size(MPI_COMM_WORLD, &size);
       tag = 202;
       MPI_Status status;
       //on each process copy all the data points to be sent 
       for(int k=0;k<N;k++)
          {
             u_lastrow[k]=gamma_p[nlocalx -1][k];
             u_firstrow[k]=gamma_p[0][k];
      //     std::cout<<u_lastrow[k]<<"\t";
           }
            std::cout<<""<<std::endl;
    //This is for periodic boundary conditions used by process 0 and process size-1
    if (rank==0)
     {
         for(int j=0;j<N;j++)
              {
                 Boundary_top[j]= gamma_p[0][j];
                // std::cout<<"Boundary_top="<< Boundary_top[j]<<std::endl;
              }
     }
     if (rank==size-1)
      {
         for(int j=0;j<N;j++)
              {
                Boundary_bottom[j]= gamma_p[nlocalx -1][j];
              }


      }
                                                                                                                                                                                       MPI_Barrier(MPI_COMM_WORLD);
    int  up_nbr =0, down_nbr=0; 
         up_nbr = rank + 1;
    if (up_nbr >= size) up_nbr = MPI_PROC_NULL;
    down_nbr = rank - 1;
    if (down_nbr < 0) down_nbr = MPI_PROC_NULL;
 //--------------exchange periodic boundary condition-----------------------------
    if (rank==0)
       {
          //  MPI_Send(u_firstrow,N,MPI_DOUBLE,size-1,tag,MPI_COMM_WORLD);
          //  MPI_Recv(u_lastrow,N,MPI_DOUBLE,size-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
          MPI_Sendrecv(Boundary_top,N,MPI_DOUBLE,size-1, tag,
                       Boundary_bottom, N, MPI_DOUBLE,size-1, tag,MPI_COMM_WORLD, &status);
          for(int j=0;j<N;j++)
             {
                 bc[j] = Boundary_bottom[j];
             }
                   

       }
    if (rank==size-1)
       {
            MPI_Sendrecv(Boundary_bottom,N,MPI_DOUBLE,0, tag,
                         Boundary_top, N, MPI_DOUBLE,0, tag,MPI_COMM_WORLD, &status);
             for(int j=0;j<N;j++)
             {
                 bc[j] = Boundary_top[j];
             }
 
       }
 //----------------end of exchange periodic boundary condition------------------------
    if((rank%2)==0)
      {
         MPI_Sendrecv(u_lastrow, N, MPI_DOUBLE, up_nbr, tag,
                      u_firstrow,N ,MPI_DOUBLE, up_nbr, tag, MPI_COMM_WORLD,  &status);

      }
    else
       {
        // std::cout<<"rank="<<rank<<std::endl;
          MPI_Sendrecv(u_lastrow, N, MPI_DOUBLE, up_nbr, tag,
                       u_firstrow,N ,MPI_DOUBLE, up_nbr, tag, MPI_COMM_WORLD,  &status);
       }
    if((rank%2)==1)
      {
    
          MPI_Sendrecv(u_firstrow, N, MPI_DOUBLE, down_nbr, tag,
                       u_lastrow,N ,MPI_DOUBLE, down_nbr, tag, MPI_COMM_WORLD,  &status);
      }
     else
      {
            MPI_Sendrecv(u_firstrow, N, MPI_DOUBLE, down_nbr, tag,
                         u_lastrow,N ,MPI_DOUBLE, down_nbr, tag, MPI_COMM_WORLD,  &status);

      }
       




    }
 
void Parallel_CahnHill2D::ComputeSecondLaplacian()
    {
       MPI_Comm_rank(MPI_COMM_WORLD, &rank);
       MPI_Comm_size(MPI_COMM_WORLD, &size);
       setIndex_p();

       if (rank==0)
          {
           // setIndex_p();
            double AP=0.0, BP=0.0;
            for(int j=0;j<N;j++)
               {
                 const int i=0;
                 AP = C1*(gamma_p(i+1,j) + bc[j] + gamma_p(i, right1x[j]) + gamma_p(i,left1x[j]) );
                 BP = C2*(gamma_p(i+1,right1x[j]) + gamma_p(i+1,left1x[j]) + bc[left1x[j]] + bc[right1x[j]]);
                 Laplacian2_p(i,j) = AP + BP;
                 
               }
             AP =0.0;BP=0.0;
            for(int i=1;i<nlocalx -1;i++)
               {
                for(int j=0;j<N;j++)
                  {
                     AP = C1*(gamma_p(i+1,j) + gamma_p(i-1,j) + gamma_p(i,right1x[j]) + gamma_p(i,left1x[j]));
                     BP = C2*(gamma_p(i+1,right1x[j]) + gamma_p(i+1,left1x[j]) + gamma_p(i-1,right1x[j]) + gamma_p(i-1,left1x[j]));
                     Laplacian2_p(i,j) = AP + BP;
                     
                   }
               }
             AP =0.0;BP=0.0;
             for(int j=0;j<N;j++)
                 {
                   const int i =nlocalx -1;
                   AP = C1*(u_firstrow[j] + gamma_p(i-1,j) + gamma_p(i,right1x[j]) + gamma_p(i,left1x[j]));
                   BP = C2*(u_firstrow[right1x[j]] + u_firstrow[left1x[j]] + gamma_p(i-1,right1x[j]) + gamma_p(i-1,left1x[j]));
                   Laplacian2_p(i,j) = AP + BP;
                   

                 }
          }

      else if(rank==size-1)
        {
           // setIndex_p();
            double AP=0.0, BP=0.0;
            for(int j=0;j<N;j++)
               {
                 const int i=0;
                 AP = C1*(gamma_p(i+1,j) + u_lastrow[j] + gamma_p(i, right1x[j]) + gamma_p(i,left1x[j]) );
                 BP = C2*(gamma_p(i+1,right1x[j]) + gamma_p(i+1,left1x[j]) + u_lastrow[left1x[j]] + u_lastrow[right1x[j]]);
                 Laplacian2_p(i,j) = AP + BP;
                 
               }
            AP =0.0;BP=0.0;
            for(int i=1;i<nlocalx -1;i++)
               {
                for(int j=0;j<N;j++)
                  {
                     AP = C1*(gamma_p(i+1,j) + gamma_p(i-1,j) + gamma_p(i,right1x[j]) + gamma_p(i,left1x[j]));
                     BP = C2*(gamma_p(i+1,right1x[j]) + gamma_p(i+1,left1x[j]) + gamma_p(i-1,right1x[j]) + gamma_p(i-1,left1x[j]));
                     Laplacian2_p(i,j) = AP + BP;
                     
                   }
               }
             AP =0.0;BP=0.0;
             for(int j=0;j<N;j++)
                 {
                     const int i =nlocalx-1;
                     AP = C1*(bc[j] + gamma_p(i-1,j) + gamma_p(i,right1x[j]) + gamma_p(i,left1x[j]) );
                     BP = C2*(bc[right1x[j]] + bc[left1x[j]] + gamma_p(i-1,right1x[j]) + gamma_p(i-1,left1x[j]));
                     Laplacian2_p(i,j) = AP + BP;
                     
                 }    

 
        }      
   else
       {
          double AP=0.0, BP=0.0;
            for(int j=0;j<N;j++)
               {
                 const int i=0;
                 AP = C1*(gamma_p(i+1,j) + u_lastrow[j] + gamma_p(i, right1x[j]) + gamma_p(i,left1x[j]) );
                 BP = C2*(gamma_p(i+1,right1x[j]) + gamma_p(i+1,left1x[j]) + u_lastrow[left1x[j]] + u_lastrow[right1x[j]]);
                 Laplacian2_p(i,j) = AP + BP;
                 
               }
             AP =0.0;BP=0.0;
              for(int i=1;i<nlocalx-1;i++)
               {
                for(int j=0;j<N;j++)
                  {
                     AP = C1*(gamma_p(i+1,j) + gamma_p(i-1,j) + gamma_p(i,right1x[j]) + gamma_p(i,left1x[j]));
                     BP = C2*(gamma_p(i+1,right1x[j]) + gamma_p(i+1,left1x[j]) + gamma_p(i-1,right1x[j]) + gamma_p(i-1,left1x[j]));
                     Laplacian2_p(i,j) = AP + BP;
                     
                   }
               }
             AP =0.0;BP=0.0;
   
              for(int j=0;j<N;j++)
                 {
                   const int i =nlocalx -1;
                   AP = C1*(u_firstrow[j] + gamma_p(i-1,j) + gamma_p(i,right1x[j]) + gamma_p(i,left1x[j]));
                   BP = C2*(u_firstrow[right1x[j]] + u_firstrow[left1x[j]] + gamma_p(i-1,right1x[j]) + gamma_p(i-1,left1x[j]));
                   Laplacian2_p(i,j) = AP + BP;
                   

                 }


                                             


      }



    }

void Parallel_CahnHill2D::ExchangeData()
 {
       MPI_Comm_rank(MPI_COMM_WORLD, &rank);
       MPI_Comm_size(MPI_COMM_WORLD, &size);
      // tag = 202;
        MPI_Status status;
    // on each process copy all the data points to be sent 
                  for(int k=0;k<N;k++)
                      {
                          u_lastrow[k]=PHI_old_p[nlocalx-1][k];
                          u_firstrow[k]=PHI_old_p[0][k];
                         // std::cout<<u_lastrow[k]<<"\t";
                        }
      //      std::cout<<""<<std::endl;
                        
            MPI_Barrier(MPI_COMM_WORLD);
     /*--- Start the exchange of data process -------*/ 
  if(rank<size-1)
    {

                if(rank==0)
                {

                 MPI_Send(u_lastrow,N,MPI_DOUBLE,rank+1,tag,MPI_COMM_WORLD);

                 MPI_Recv(u_firstrow,N,MPI_DOUBLE,rank+1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
                 }
                if(rank>0)
                {
                // std::cout<<"if statement rank="<<rank<<std::endl;
                // std::cout<<""<<std::endl;
                   MPI_Send(u_firstrow,N,MPI_DOUBLE,rank-1,10,MPI_COMM_WORLD);
                   MPI_Send(u_lastrow,N,MPI_DOUBLE,rank+1,10,MPI_COMM_WORLD);
                 }
                
      }
   if(rank==size-1)
      {
                MPI_Send(u_firstrow,N,MPI_DOUBLE,rank-1,tag,MPI_COMM_WORLD);

                MPI_Recv(u_lastrow,N,MPI_DOUBLE,rank-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      }
   if(rank<size-1)
        {
             if(rank>0)
                {
                  MPI_Recv(u_lastrow,N,MPI_DOUBLE,rank-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);

                  MPI_Recv(u_firstrow,N,MPI_DOUBLE,rank+1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);

                }
        }
   MPI_Barrier(MPI_COMM_WORLD);

 /*------------Periodic boundary condition on the y-axis (vertically), only rank=0 and rank =size-1 participate--------------*/
  if(rank==0)
    {
        
        
           for(int j=0;j<N;j++)
              {
                 Boundary_top[j]= PHI_old_p[0][j];
                 std::cout<<"Boundary_top="<< Boundary_top[j]<<std::endl;
              }
        

       // MPI_Sendrecv(u_firstrow,N,MPI_DOUBLE,size-1, 203,
       //               u_lastrow, N, MPI_DOUBLE,
       //              size-1, 203,MPI_COMM_WORLD, &status);
       //
 
       // MPI_Send(u_firstrow,N,MPI_DOUBLE,size-1,tag,MPI_COMM_WORLD);
       // MPI_Recv(u_lastrow,N,MPI_DOUBLE,size-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
       
        MPI_Send(Boundary_top,N,MPI_DOUBLE,size-1,tag,MPI_COMM_WORLD);
        MPI_Recv(Boundary_bottom,N,MPI_DOUBLE,size-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
   
    }
  if(rank==size-1)
    { 
     // copy the last  row to be sent into a long array
      
        
           for(int j=0;j<N;j++)
              {
                Boundary_bottom[j]= PHI_old_p[nlocalx-1][j];
              }
     

       // MPI_Sendrecv(u_lastrow,N,MPI_DOUBLE,0, 2,
       //               u_firstrow, N, MPI_DOUBLE,
       //              0, 2,MPI_COMM_WORLD, &status); 
       
            
     //  MPI_Recv(u_firstrow,N,MPI_DOUBLE,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
     //  MPI_Send(u_lastrow,N,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
   //
       
       
       MPI_Recv(Boundary_top,N,MPI_DOUBLE,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
       MPI_Send(Boundary_bottom,N,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
    }

/*---------end of vertical periodic boundary exchange --------------*/                                                                                                              
   MPI_Barrier(MPI_COMM_WORLD);
//calculate partial part of the laplacian  and gamma function with exhange data.
if (rank==0)
   {
     // int i;
      setIndex_p();
      double AP1=0.0, BP1=0.0;
      double  SUM_L=0.0;
      std::cout<<"C1="<<C1<<std::endl;
      for(int j=0;j<N;j++)
         {
               const int i =0; //first row on process 0 
                std::cout<<"D in parallel="<<D<<std::endl;              
                AP1 = C1*(PHI_old_p[i+1][j] + Boundary_bottom[j] + PHI_old_p[i][right1x[j]] + PHI_old_p[i][left1x[j]] );
                BP1 = C2*(PHI_old_p[i+1][right1x[j]] + PHI_old_p[i+1][left1x[j]] + Boundary_bottom[right1x[j]] + Boundary_bottom[left1x[j]]);
                SUM_L = AP1 + BP1;
                gamma_p(i,j)=CahnHill2D::g(PHI_old_p(i,j))+ D*(SUM_L-PHI_old_p(i,j))-PHI_old_p(i,j);

         }
             
        
      // std::cout<<"(rank, ATP1)="<<std::setw(4)<<rank<<std::setw(12)<< SUM_L<<std::endl;
      }
if (rank==size-1)
   {
             setIndex_p();
             double AP1=0.0, BP1=0.0;  
             double  SUM_L=0.0;

       for(int j=0; j<N; j++)
          {
             const int i=nlocalx-1; //last row on process size-1
             AP1 = C1*(Boundary_top[j] + PHI_old_p[i-1][j] + PHI_old_p[i][right1x[j]] + PHI_old_p[i][left1x[j]] );
             BP1 = C2*(Boundary_top[right1x[j]] + Boundary_top[left1x[j]] + PHI_old_p[i-1][right1x[j]] + PHI_old_p[i-1][left1x[j]]);
             SUM_L = AP1 + BP1;
             gamma_p(i,j)=CahnHill2D::g(PHI_old_p(i,j))+ D*(SUM_L-PHI_old_p(i,j))-PHI_old_p(i,j);


         }


   }
  MPI_Barrier(MPI_COMM_WORLD);

    //compute the laplacian and gamma function for  last row on each process except process size-1
 if(rank<size-1)
   { 
          setIndex_p();

         double AP1=0.0, BP1=0.0;
         double  SUM_L=0.0;
         const int i=nlocalx-1;
       for(int j =0;j<N;j++)
        { 
         AP1 = C1*(u_firstrow[j] + PHI_old_p(i-1,j) + PHI_old_p(i,right1x[j]) + PHI_old_p(i,left1x[j]) );
         BP1 = C2*(u_firstrow[right1x[j]] + u_firstrow[left1x[j]] + PHI_old_p(i-1,right1x[j]) + PHI_old_p[i-1][left1x[j]]);
         SUM_L = AP1 + BP1;
         gamma_p(i,j)=CahnHill2D::g(PHI_old_p(i,j))+ D*(SUM_L-PHI_old_p(i,j))-PHI_old_p(i,j);
     
       }
       std::cout<<"I am rank, computing lastrow"<<std::setw(8)<<rank<<std::endl;   
   }
//compute the laplacian and gamma function for first row on each process except process 0
 if (rank>0)
    {

      setIndex_p();
      double AP1=0.0, BP1=0.0;
      double  SUM_L=0.0;
//      std::cout<<"C1="<<C1<<std::endl;
      for(int j=0;j<Ny;j++)
         {
               const int i =0; //first row on process 0 

                AP1 = C1*(PHI_old_p[i+1][j] + u_lastrow[j] + PHI_old_p[i][right1x[j]] + PHI_old_p[i][left1x[j]] );
                BP1 = C2*(PHI_old_p[i+1][right1x[j]] + PHI_old_p[i+1][left1x[j]] + u_lastrow[right1x[j]] + u_lastrow[left1x[j]]);
                SUM_L = AP1 + BP1;
                gamma_p(i,j)=CahnHill2D::g(PHI_old_p(i,j))+ D*(SUM_L-PHI_old_p(i,j))-PHI_old_p(i,j);

         }
   


    }
//compute the laplacian and gamma function on all processes except, leave out first row and last row since they have been calculated
 MPI_Barrier(MPI_COMM_WORLD);
 double AP_p=0.0, BP_p=0.0, ATP_p=0.0;
 for(int i=1;i<nlocalx-1;i++)
   {
    for(int j=0;j<N;j++)
       {

         
                 AP_p=C1*(PHI_old_p(i+1,j) + PHI_old_p(i-1,j)
                        + PHI_old_p(i,right1x[j]) + PHI_old_p(i,left1x[j]));
                 BP_p=C2*(PHI_old_p(i+1,right1x[j]) + PHI_old_p(i+1,left1x[j])
                        +PHI_old_p(i-1,right1x[j]) + PHI_old_p(i-1,left1x[j]));
                 ATP_p = AP_p + BP_p;
                 gamma_p(i,j)=CahnHill2D::g(PHI_old_p(i,j))+ D*(ATP_p-PHI_old_p(i,j))-PHI_old_p(i,j);


       }

   }
           
       
    }

// set second laplacian SetSecondLaplacian2()
void Parallel_CahnHill2D:: setSecond_laplacian()
    {

        //each process copies gamma values for the first and last row. These values have to sent to other proecess
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Status status;
        
        for(int k=0;k<N;k++)
            {
               array_gamma_top[k]=gamma_p[0][k];
               array_gamma_bottom[k]=gamma_p[nlocalx-1][k];
                         // std::cout<<u_lastrow[k]<<"\t";
             }
  
       
             /*--- Start the exchange of data process -------*/ 
   if(rank<size-1)
    {

                if(rank==0)
                {

                 MPI_Send(array_gamma_bottom,N,MPI_DOUBLE,rank+1,tag,MPI_COMM_WORLD);

                 MPI_Recv(array_gamma_top,N,MPI_DOUBLE,rank+1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
                 }
                if(rank>0)
                {
                // std::cout<<"if statement rank="<<rank<<std::endl;
                // std::cout<<""<<std::endl;
                   MPI_Send(array_gamma_top,N,MPI_DOUBLE,rank-1,10,MPI_COMM_WORLD);
                   MPI_Send(array_gamma_bottom,N,MPI_DOUBLE,rank+1,10,MPI_COMM_WORLD);
                 }
                
      }
   if(rank==size-1)
      {
                MPI_Send(array_gamma_top,N,MPI_DOUBLE,rank-1,tag,MPI_COMM_WORLD);

                MPI_Recv(array_gamma_bottom,N,MPI_DOUBLE,rank-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      }
   if(rank<size-1)
        {
             if(rank>0)
                {
                  MPI_Recv(array_gamma_bottom,N,MPI_DOUBLE,rank-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);

                  MPI_Recv(array_gamma_top,N,MPI_DOUBLE,rank+1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);

                }
        }
   MPI_Barrier(MPI_COMM_WORLD);

 /*------------Periodic boundary condition for gamma values, only rank=0 and rank =size-1 participate--------------*/


  if(rank==0)
    {
        
        
           for(int j=0;j<N;j++)
              {
                 Boundary_top1[j]= gamma_p[0][j];
                 //std::cout<<"Boundary_top="<< Boundary_top[j]<<std::endl;
              }
        



      MPI_Send(&Boundary_top1[0],N,MPI_DOUBLE,size-1,tag,MPI_COMM_WORLD);
      MPI_Recv(&Boundary_bottom1[0],N,MPI_DOUBLE,size-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
   
    }
  if(rank==size-1)
    { 
     // copy the last  rows to be sent into a long array
      
        
           for(int j=0;j<N;j++)
              {
                Boundary_bottom1[j]= gamma_p[nlocalx-1][j];
              }
        
     MPI_Recv(&Boundary_top1[0],N,MPI_DOUBLE,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
     MPI_Send(&Boundary_bottom1[0],N,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
    }
MPI_Barrier(MPI_COMM_WORLD);

    if (rank==0)
   {
     // int i;
      setIndex_p();
      double AP=0.0, BP=0.0;
      //double  SUM_L=0.0;
     // std::cout<<"C1="<<C1<<std::endl;
      for(int j=0;j<N;j++)
         {
               const int i =0; //first row on process 0 
              
            AP=C1*(gamma_p(i+1,j) + Boundary_bottom1[j]+ gamma_p(i,right1x[j]) + gamma_p(i,left1x[j]));
            BP=C2*(Boundary_bottom1[right1x[j]] + Boundary_bottom1[left1x[j]] + gamma_p(i+1,right1x[j]) + gamma_p(i+1,left1x[j]));
            Laplacian2_p(i,j) = AP + BP;

         }
     }
   if (rank==size-1)
      {
       setIndex_p();
       double AP=0.0, BP=0.0;
      //double  SUM_L=0.0;
      //std::cout<<"C1="<<C1<<std::endl;
      for(int j=0;j<N;j++)
         {
               const int i =nlocalx-1; //last row on process size-1
              
            AP=C1*(Boundary_top1[j] + gamma_p(i-1,j)+ gamma_p(i,right1x[j]) + gamma_p(i,left1x[j]));
            BP=C2*(gamma_p(i-1,right1x[j]) + gamma_p(i-1,left1x[j]) + Boundary_top1[right1x[j]] + Boundary_top1[left1x[j]]);
            Laplacian2_p(i,j) = AP + BP;

         }

      }
  //compute the gamma function and second laplacian using the last row on each process except process size-1
  if (rank<size-1)
     {
        setIndex_p();
       double AP=0.0, BP=0.0;
      //double  SUM_L=0.0;
     // std::cout<<"C1="<<C1<<std::endl;
      for(int j=0;j<N;j++)
         {
               const int i =nlocalx-1; //last row on process size-1
              
            AP=C1*(array_gamma_top[j] + gamma_p(i-1,j)+ gamma_p(i,right1x[j]) + gamma_p(i,left1x[j]));
            BP=C2*(gamma_p(i-1,right1x[j]) + gamma_p(i-1,left1x[j]) + array_gamma_top[right1x[j]] + array_gamma_top[left1x[j]]);
            Laplacian2_p(i,j) = AP + BP;

         }

     }
 //compute the gamma function and the second laplacian using first row on each process except rank 0
 if (rank>0)
    {
      setIndex_p();
      double AP=0.0, BP=0.0;
      //double  SUM_L=0.0;
      std::cout<<"C1="<<C1<<std::endl;
      for(int j=0;j<N;j++)
         {
               const int i =0; //first row on process 0 
              
            AP=C1*(gamma_p(i+1,j) + array_gamma_bottom[j]+ gamma_p(i,right1x[j]) + gamma_p(i,left1x[j]));
            BP=C2*(array_gamma_bottom[right1x[j]] + array_gamma_bottom[left1x[j]] + gamma_p(i+1,right1x[j]) + gamma_p(i+1,left1x[j]));
            Laplacian2_p(i,j) = AP + BP;

         } 

    }
//compute Gamma function and second Laplacian for all the points except first and last rows
double AP=0.0, BP=0.0;
 for(int i=1;i<nlocalx-1;i++)
   {
    for(int j=0;j<N;j++)
       {
          AP=C1*(gamma_p(i+1,j) + gamma_p(i-1,j)+ gamma_p(i,right1x[j]) + gamma_p(i,left1x[j]));
          BP=C2*(gamma_p(i-1,right1x[j]) + gamma_p(i-1,left1x[j]) + gamma_p(i+1,right1x[j]) + gamma_p(i+1,left1x[j]));
          Laplacian2_p(i,j) = AP + BP;

       }
    }
}


void Parallel_CahnHill2D::FiniteDifferenceScheme_p()
{
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    //double r=0.5;
    cout<<"in finite difference scheme parallel"<<endl;
    for(int i=0;i<nlocalx;i++)
    {
        for(int j=0;j<N;j++)
        {
          // PHI[i][j]= PHI_old(i,j) + delta_t*(Laplacian2[i][j]-gamma(i,j) + B*(1.0)*PHI_old(i,j));
            PHI_p[i][j]= PHI_old_p(i,j) - B*(PHI_old_p(i,j)-1.0 + 2*r)+ gamma_p(i,j)-Laplacian2_p(i,j);
            
        }
        
    }
    
}

void Parallel_CahnHill2D::UpdateSolution_p()
    {
       
       for(int i=0;i<nlocalx;i++)
          {
           for(int j=0;j<N;j++)
              {
                PHI_old_p[i][j]=PHI_p[i][j];
             }

          }

    }

void Parallel_CahnHill2D::WriteToFile_MPI()
    {
       
       char const *fmt ="%20.8f";
       char const *endfmt="%20.8f\n";
       MPI_Status status;
       MPI_File fileh;
       MPI_Offset disp;
       MPI_Datatype num_as_string;
       MPI_Datatype localarray;
       const int charspernum=16;
       MPI_Comm_rank(MPI_COMM_WORLD, &rank);
       MPI_Comm_size(MPI_COMM_WORLD, &size);
      
       MPI_Type_contiguous(charspernum, MPI_CHAR,&num_as_string);
       MPI_Type_commit(&num_as_string);
       /*----convert data into txt-----*/
       char *data_as_txt=(char *)malloc(nlocalx*Ny*charspernum*sizeof(char));
       int count =0;
       for(int i=0;i<nlocalx;i++)
          {
           for(int j=0;j<Ny-1;j++)
              {
                sprintf(&data_as_txt[count*charspernum],fmt,PHI_p[i][j]);
                count++;
              }
             sprintf(&data_as_txt[count*charspernum],endfmt,PHI_p[i][Ny-1]);
             count++;

          }
 
       int globalsizes[2]={Nx, Ny};
       int localsizes[2]={nlocalx,Ny};
       int starts[2]={rank*nlocalx,0};
       int order = MPI_ORDER_C;
       MPI_Type_create_subarray(2,globalsizes,localsizes,starts,order,MPI_DOUBLE,&localarray);
       MPI_Type_commit(&localarray);
 
       MPI_File_open(MPI_COMM_WORLD, "mpifile.dat", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fileh);
       //MPI_File_set_view(fileh, 0, MPI_CHAR, localarray, "native", MPI_INFO_NULL);
      // MPI_File_write_all(fileh,data_as_txt,nlocalx*Ny,num_as_string,&status);
       MPI_File_set_view(fileh, rank*nlocalx*sizeof(double), MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
       MPI_File_write_all(fileh,PHI_p,nlocalx*Ny,MPI_DOUBLE,&status);
       MPI_File_close(&fileh);
       MPI_Type_free(&localarray);
   //    MPI_Type_free(&num_as_string);


     


   }

void Parallel_CahnHill2D::SendToMaster(int count)
    {
       
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        double *U;
        for(int i=0;i<nlocalx;i++)
           {
                for(int j=0;j<Ny;j++)
                   {
                           u_store[i*Ny+j]=PHI_p[i][j];
                   }
                }

        if(rank==0)
          {
                  U = new double[Nx*Ny];
                  std::memset(U, 0, Nx*Ny*sizeof(double));
          }
       MPI_Gather( u_store, nlocalx*Ny, MPI_DOUBLE, U, nlocalx*Ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (rank==0)
           {
                  FILE *file = fopen(make_output_filename(count).c_str(), "w");
                   for(int i=0;i<Nx;i++)
                      {
                            for(int j=0;j<Ny;j++)
                              {

                                 fprintf(file,"%15.11f \t",U[i*Ny+j]);

                              }
                              fprintf(file,"\n");
                     }
                fclose(file);

                 delete []U;
           }

   //  save_vtk(U, Nx,Ny);
    }
void Parallel_CahnHill2D::parallel_solver()
    {
       Initialize_parallel();
     //   ReadFile(output_0.dat);

       int count =0;
       double t=0.0;
       while(t<1000)
           {
            // std::cout<<"in while loop solver"<<std::endl;
             ExchangeData();
             setSecond_laplacian();
           //  ExchangeData2();
           //  ComputeLaplacianBase_p();
           //  ExchangeGammaData();
           //  ComputeSecondLaplacian();
             FiniteDifferenceScheme_p();
             UpdateSolution_p();
            

          //  if (count==100|| count==0)
           //    {
            //     SendToMaster( count);
                 //result( count);
            //   }
        //t+=delta_t;
                t+=1.0;
                count++;
        //                        cout<<"#"<<"\t"
           }
           //SendToMaster( count);
           WriteToFile();
         //  WriteToFile_MPI();


            cout<<""<<endl;
    
     

    }


void Parallel_CahnHill2D::WriteToFile()
    {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);

      char fname[200];
      sprintf(fname,"solution%d.dat",rank);
      FILE* fp;
      fp = fopen(fname,"w");
      for(int i=0;i<nlocalx;i++)
         {
           for(int j =0; j<N;j++)
             {
                fprintf(fp,"%20.16f", PHI_p[i][j]);
               
             }
           fprintf(fp,"\n");
   
        }
     fclose(fp);  


    }

void Parallel_CahnHill2D::ReadFile(std::ifstream& infile)
    {
     /*//  std::ifstream& output_0.dat;
       MPI_Comm_rank(MPI_COMM_WORLD, &rank);
       MPI_Comm_size(MPI_COMM_WORLD, &size);
       for(int i=0;i<nlocalx;i++)
          {
           for(int j=0;j<Ny;j++)
              {

                 infile>>PHI_old_p[i][j];
              }
           }
 
       
          */


         
            MPI_Status status;
         
         
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &size);
          
          
            tag = 201;
            low=rank*nlocalx;
          //high=low+nlocalx;
          double **u1,*u2;
           u1=new double*[Nx];
           u2 =new double[Nx*Ny];
           std::memset(u2, 0, Nx*Ny*sizeof(double));
           for(int i=0;i<Nx;i++)
              {       
                 u1[i]=new double[Ny];
              }
 
               for(int i=0;i<Nx;i++)
                  {
                   for(int j=0;j<Ny;j++)
                      {
                        u1[i][j]=0.0;
                       }
                   }
                   
          


              for(int i=0;i<Nx;i++)
                 {
                  for(int j=0;j<Ny;j++)
                     {
              
                         infile>>u2[i*Ny+j];
                       //  u2[i*Ny+j]=u1[i][j];
                      //  printf("rank=%d, u2[%d]=%lf\n",rank,i*Ny+j,u2[i*Ny+j]);
                     }
                 }
          
       
           
          if(rank==0)
          {

                       
             
              for(int i=0;i<Nx;i++)
                 {
                  for(int j=0;j<Ny;j++)
                     {

                         infile>>u1[i][j];
                        // u2[i*Ny+j]=u1[i][j];
                       // printf("rank=%d, u2[%d]=%lf\n",rank,i*Ny+j,u2[i*Ny+j]);
                     }
                 }
        

                                                                                                                                                                                                                                                                                                                                                                                          
         for(int i=1;i<size;i++)
            {
              int offset1=i*nlocalx;
              MPI_Send(&u1[offset1][0],nlocalx*Ny,MPI_DOUBLE,i,tag,MPI_COMM_WORLD);
           }
                 
  
          }
       else
          {
              MPI_Recv(&u1[rank*nlocalx][0],nlocalx*Ny,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&status);
                int k=0;
             for(int i=rank*nlocalx;i<rank*nlocalx+nlocalx;i++)
                {
                  
                  for(int j=0;j<Ny;j++)
                     {
                     //   PHI_old_p[k][j]=u2[i*Ny+j];
                      //  printf("rank=%d, u2[%d][%d]=%lf\n",rank,i,j,u2[i*Ny+j]);
                     //   printf("rank=%d, phi_old_p[%d][%d]=%lf\n",rank,k,j,PHI_old_p[k][j]);
                    }
                //  k++;
               }
         }
       int k=0;
      for(int i =rank*nlocalx; i<rank*nlocalx+nlocalx; i++)
        
          {
            for(int j=0;j<Ny; j++)
               {
                 
                  PHI_old_p[k][j]=u2[i*Ny+j];
                  printf("rank=%d, phi_old[%d][%d]=%lf\n",rank,k,j,PHI_old_p[k][j]);
                }
             k++;
           }
  
 
                                                                                                                                                 for(int i =0;i<Nx;i++)                                                                                                                                               {     
                      delete []u1[i];
                                                                                                                                                                                          }
               
               
            delete []u1;
            delete []u2;
                                                                                                                                                                                                                                                                                                                                              

     }
