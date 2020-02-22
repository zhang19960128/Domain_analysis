#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <string>
#include <cmath>
#include <fstream>
#include <list>
double average(std::list<double>& listA){
  double sum=0.0;
  int count=0;
  for(std::list<double>::iterator a=listA.begin();a!=listA.end();a++){
    sum=sum+*a;
    count=count+1;
  }
  if(count==0){
    return 0.0;
  }
  else{
    return sum/count;
  }
}
int next(int i,int direction,int period){
  if(direction==0){
    return (i+1)%period;
  }
  else if(direction==1){
    return (i+period)%(period*period);
  }
  else{
    return (i+period*period)%(period*period*period);
  }
}
/*count the domain from X direction,polarization along k direction*/
size_t findstartx(int cell,int ny,int nz,int k,std::list<double>& px,double* reducepolar,double decay_rate){
	  size_t nindex;
		size_t record_start;
		px.clear();
    for(size_t j=0;j<cell;j++){
     nindex=j+ny*cell+nz*cell*cell;
     if(px.size()==0){
      px.push_back(reducepolar[k+nindex*3]);
     }
     else{
      if(std::fabs(average(px)-reducepolar[k+nindex*3])/std::fabs(average(px))>decay_rate){
        record_start=j;
        break;
      }
			else{
				px.push_back(reducepolar[k+nindex*3]);
			}
     }
    }
		return record_start;
}
/*count the the domain from Y direction,polar along the k direciton*/
size_t findstarty(int cell,int nx,int nz,int k,std::list<double>& py,double* reducepolar,double decay_rate){
	size_t nindex;
	size_t record_start;
	py.clear();
    for(size_t j=0;j<cell;j++){
     nindex=nx+j*cell+nz*cell*cell;
     if(py.size()==0){
      py.push_back(reducepolar[k+nindex*3]);
     }
     else{
      if(std::fabs(average(py)-reducepolar[k+nindex*3])>decay_rate*std::fabs(average(py))){
        record_start=j;
        break;
      }
			else{
				py.push_back(reducepolar[k+nindex*3]);
			}
     }
		}
		return record_start;
}
/*count the the domain from Z direction,polar along the k direciton*/
size_t findstartz(int cell,int nx,int ny,int k,std::list<double>& pz,double* reducepolar,double decay_rate){
	size_t nindex;
	size_t record_start;
	pz.clear();
   for(size_t j=0;j<cell;j++){
   nindex=(nx+ny*cell+j*cell*cell);
   if(pz.size()==0){
    pz.push_back(reducepolar[k+nindex*3]);
   }
   else{
    if(std::fabs(average(pz)-reducepolar[k+nindex*3])>decay_rate*std::fabs(average(pz))){
      record_start=j;
      break;
    }
		else{
			pz.push_back(reducepolar[k+nindex*3]);
		}
   }
  }
	return record_start;
}
/*count along the x direciton, polar along the k direction*/
void countdomainX(int cell,int ny,int nz,int k,size_t record_start,std::list<double>& px,double* reducepolar,int* domainsizecount,double decay_rate){
	px.clear();
	size_t nindex;
  for(size_t j=record_start;j<cell+record_start;j++){
      nindex=j%cell+ny*cell+nz*cell*cell;
      if(px.size()==0){
        px.push_back(reducepolar[k+nindex*3]);
      }
      else{
        if(std::fabs(average(px)-reducepolar[k+nindex*3])>decay_rate*std::fabs(average(px))){
          domainsizecount[px.size()]=domainsizecount[px.size()]+1;
          px.clear();
        }
        else{
          px.push_back(reducepolar[k+nindex*3]);
        }
      }
    }
}
/*Count the domain along the Y direction, polar along the k direction*/
void countdomainY(int cell,int nx,int nz,int k,size_t record_start,std::list<double>& py,double* reducepolar,int* domainsizecount,double decay_rate){
	py.clear();
	size_t nindex;
    for(size_t j=record_start;j<cell+record_start;j++){
      nindex=(nx+(j%cell)*cell+nz*cell*cell);
      if(py.size()==0){
        py.push_back(reducepolar[k+nindex*3]);
      }
      else{
        if(std::fabs(average(py)-reducepolar[k+nindex*3])>decay_rate*std::fabs(average(py))){
          domainsizecount[py.size()]=domainsizecount[py.size()]+1;
          py.clear();
        }
        else{
          py.push_back(reducepolar[k+nindex*3]);
        }
      }
    }
}
/*Count the domain along the Z direction, polar along the k direction*/
void countdomainZ(int cell,int nx,int ny,int k,size_t record_start,std::list<double>& pz,double* reducepolar,int* domainsizecount,double decay_rate){
	pz.clear();
	size_t nindex;
       for(size_t j=record_start;j<cell+record_start;j++){
      nindex=(nx+ny*cell+(j%cell)*cell*cell);
      if(pz.size()==0){
        pz.push_back(reducepolar[k+nindex*3]);
      }
      else{
        if(std::fabs(average(pz)-reducepolar[k+nindex*3])>decay_rate*std::fabs(average(pz))){
          domainsizecount[pz.size()]=domainsizecount[pz.size()]+1;
          pz.clear();
        }
        else{
          pz.push_back(reducepolar[k+nindex*3]);
        }
      }
    }
}
int main(){
  MPI_Init(NULL,NULL);
  int simulationtime=1800000/200;
	int initial_offset=0;
  int cell=48;
  int world_rank,world_size;
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  MPI_File mpifile;
  MPI_Status status;
  MPI_Offset offset;
	std::list<double> px;
	std::list<double> py;
	std::list<double> pz;
  double* polaraverage=new double [3*cell*cell*cell];
  double* reducepolar=new double [3*cell*cell*cell];
	double* tempsum=new double [3*cell*cell*cell];
  double* localpolar=new double [3*cell*cell*cell];
  for(size_t i=0;i<3*cell*cell*cell;i++){
    polaraverage[i]=0.0;
    reducepolar[i]=0.0;
		tempsum[i]=0.0;
 }
  std::string filename="local_polar.bin";
  MPI_File_open(MPI_COMM_WORLD,filename.c_str(),MPI_MODE_RDONLY,MPI_INFO_NULL,&mpifile);
  for(size_t frame=world_rank+initial_offset;frame<simulationtime;frame=frame+world_size){
    offset=(3*cell*cell*cell*frame)*sizeof(double);
    MPI_File_read_at(mpifile,offset,localpolar,3*cell*cell*cell,MPI::DOUBLE,&status);
		for(size_t i=0;i<3*cell*cell*cell;i++){
			tempsum[i]=tempsum[i]+localpolar[i];
		}
  }
	MPI_File_close(&mpifile);
  int* domainsizecount=new int [cell];
  int* domainsizecount_reduce=new int [cell];
  MPI_Allreduce(tempsum,reducepolar,3*cell*cell*cell,MPI::DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	for(size_t i=0;i<3*cell*cell*cell;i++){
		reducepolar[i]=reducepolar[i]/simulationtime;
	}
  /*starting from X direction*/
  int nx,ny,nz,nindex;
  nx=0;
  size_t record_start;
	size_t record_start_direction[3]={0,0,0};
  double decay_rate=0.7;
  for(size_t i=0;i<cell;i++){
    domainsizecount[i]=0;
    domainsizecount_reduce[i]=0;
  }
  for(size_t i=world_rank;i<cell*cell;i=i+world_size){
    ny=i%cell;
    nz=(i-ny)/cell;
		if(world_rank==0){
		std::cout<<i<<std::endl;
		}
		for(size_t k=0;k<3;k++){
    px.clear();
    record_start_direction[k]=findstartx(cell,ny,nz,k,px,reducepolar,decay_rate);
		px.clear();
		countdomainX(cell,ny,nz,k,record_start_direction[k],px,reducepolar,domainsizecount,decay_rate);
		}
  }
  //starting from y direction
  for(size_t i=world_rank;i<cell*cell;i=i+world_size){
    nx=i%cell;
    nz=(i-nx)/cell;
		if(world_rank==0){
		std::cout<<i<<std::endl;
		}
    for(size_t k=0;k<3;k++){
    py.clear();
		record_start_direction[k]=findstarty(cell,nx,nz,k,py,reducepolar,decay_rate);
		py.clear();
		countdomainY(cell,nx,nz,k,record_start_direction[k],py,reducepolar,domainsizecount,decay_rate);
		}
    }
  //z direction
  for(size_t i=world_rank;i<cell*cell;i=i+world_size){
    nx=i%cell;
    ny=(i-nx)/cell;
		if(world_rank==0){
		std::cout<<i<<std::endl;
		}
		for(size_t k=0;k<3;k++){
			pz.clear();
		record_start_direction[k]=findstartz(cell,nx,ny,k,pz,reducepolar,decay_rate);
		pz.clear();
		countdomainZ(cell,nx,ny,k,record_start_direction[k],pz,reducepolar,domainsizecount,decay_rate);
		}
  }
  MPI_Allreduce(domainsizecount,domainsizecount_reduce,cell,MPI::INT,MPI_SUM,MPI_COMM_WORLD);
  if(world_rank==0){
		std::fstream fs;
		fs.open("domainsize.txt",std::fstream::out);
    for(size_t i=0;i<cell;i++){
      fs<<i<<" "<<domainsizecount_reduce[i]<<std::endl;
    }
		fs.close();
  }
  delete [] localpolar;
  delete [] polaraverage;
	delete [] tempsum;
  delete [] reducepolar;
  px.clear();
  py.clear();
  pz.clear();
  MPI_Finalize();
}
