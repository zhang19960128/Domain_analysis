#include <mpi.h>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <list>
#include <cstdlib>
#include <iostream>
#include <time.h>
#include <stdio.h>
#include <cstring>
#include <vector>
#include "interface.h"
#include "constant.h"
#define PI 3.141592653
int main(int argc,char* argv[]){
  MPI_Init(NULL,NULL);
	int world_rank,world_size;
	MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&world_size);
	clock_t start=clock();
	double distance=0.0;
	double temp_double=0.0;
	int cell=48;
	int inteval=150;
	int frame=1500000/200/inteval;
	int* pointA;
	int* pointB;
	updatefileinfo(argv[1],cell);
	MPI_File fh;
	MPI_File_open(MPI_COMM_WORLD,filesys::polar_direction_file,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
	//pair pattern (distance,timedelay correlation)
	MPI_Offset offsetone,offsettwo;
	MPI_Status status;
	double directA[3]={0.0,0.0,0.0};
	double directB[3]={0.0,0.0,0.0};
	double* timedelay=new double [frame];
	MPI_File correlation;
	MPI_File_open(MPI_COMM_WORLD,filesys::polar_correlation_file,MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&correlation);
	double* frameA=new double [cell*cell*cell*3];
	double* frameB=new double [cell*cell*cell*3];
	double* sum=new double [filesys::pairsize];
	double* angle_pair_local=new double [filesys::pairsize];
	double* angle_delay_local=new double [filesys::pairsize];
	double* angle_pair_global=new double [filesys::pairsize];
	for(size_t i=0;i<filesys::pairsize;i++){
		angle_pair_local[i]=0.0;
		angle_pair_global[i]=0.0;
		angle_delay_local[i]=0.0;
	}
	//parallel over the time delay
	for(size_t j=world_rank;j<frame;j=j+world_size){
		for(size_t k=0;k<filesys::pairsize;k++){
			sum[k]=0.0;
		}
		for(size_t k=0;k<frame;k++){
			//(k,k+j)
			if(k+j>=frame){
			}
			else{
				offsetone=(k*inteval*cell*cell*cell)*3*sizeof(double);
				offsettwo=((k+j)*inteval*cell*cell*cell)*3*sizeof(double);
				MPI_File_read_at(fh,offsetone,frameA,3*cell*cell*cell,MPI::DOUBLE,&status);
				MPI_File_read_at(fh,offsettwo,frameB,3*cell*cell*cell,MPI::DOUBLE,&status);
				for(size_t m=0;m<filesys::pairsize;m++){
					for(size_t n=0;n<3;n++){
						sum[m]=sum[m]+frameA[filesys::pairA[m]*3+n]*frameB[filesys::pairB[m]*3+n];
					}
				}
			}
		}
			for(size_t m=0;m<filesys::pairsize;m++){
				sum[m]=sum[m]/frame;
			}
		for(size_t j=0;j<filesys::pairsize;j++){
			angle_delay_local[j]=180/PI*acos(sum[j]);
		}
		for(size_t j=0;j<filesys::pairsize;j++){
			angle_pair_local[j]=angle_pair_local[j]+angle_delay_local[j];
		}
	offsetone=(j*filesys::pairsize)*sizeof(double);
	MPI_File_write_at(correlation,offsetone,sum,filesys::pairsize,MPI::DOUBLE,&status);
	}
	MPI_File_close(&fh);
	MPI_File_close(&correlation);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(angle_pair_local,angle_pair_global,filesys::pairsize,MPI::DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for(size_t i=0;i<filesys::pairsize;i=i+1){
		angle_pair_global[i]=angle_pair_global[i]/frame;
	}
	if(world_rank==0){
		std::fstream Corre_angle;
	  Corre_angle.open(filesys::correlation_angle_file,std::fstream::out);
	for(size_t i=0;i<filesys::pairsize;i=i+1){
		Corre_angle<<angle_pair_global[i]<<std::endl;
	}
	Corre_angle.close();
	std::cout<<"I finished calculating correlation angles"<<std::endl;
	}
	else{
	}
	MPI_Barrier(MPI_COMM_WORLD);
  if(world_rank==0){
	std::fstream correlation_distance;
	correlation_distance.open(filesys::correlation_distance_file,std::fstream::out);
	for(size_t j=0;j<filesys::pairsize;j=j+1){
		distance=0.0;
		pointA=changeindex(filesys::pairA[j],cell);
		pointB=changeindex(filesys::pairB[j],cell);
		for(size_t t=0;t<3;t++){
			temp_double=pointA[t]-pointB[t]-round((pointA[t]-pointB[t]+0.0)/cell)*cell;
			distance=distance+temp_double*temp_double;
		}
		distance=sqrt(distance);
		distance=distance*filesys::celllength;
		correlation_distance<<distance<<std::endl;
	}
	correlation_distance.close();
	clock_t duration=clock()-start;
	std::cout<<((double)duration)/CLOCKS_PER_SEC<<" seconds"<<std::endl;
	}
	else{
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	std::cout<<"I finished"<<std::endl;
	exit(0);
}
