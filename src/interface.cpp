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
int* changeindex(int index,int cell){
	int* re=new int[3];
	re[2]=floor(index/(cell*cell));
	index=index-re[2]*cell*cell;
	re[1]=floor(index/cell);
	re[0]=index-cell*re[1];
	return re;
}
std::vector<std::string> split(std::string s, std::string delimiter) {
		    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
				std::string token;
				std::vector<std::string> res;
				while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos) {
							token = s.substr(pos_start, pos_end - pos_start);
							pos_start = pos_end + delim_len;
							res.push_back(token);
				}
				res.push_back(s.substr(pos_start));
			  return res;
}
void updatefileinfo(char* fileinput,int cell){
	std::string temp;
	std::stringstream ss;
	int atomid;
	std::list<int> pairlist;
	std::list<int> atomlist;
	std::fstream fs;
	int* pointA;
	int* pointB;
	clock_t start=clock();
	double distance;
	double temp_double;
	int world_rank,world_size;
	MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  if(world_rank==0){
		std::vector<std::string> tempstringvector;
	fs.open(fileinput,std::fstream::in);
	while(getline(fs,temp)){
		if(temp.find("filesys::cutin")!=std::string::npos){
			tempstringvector=split(temp,"=");
			filesys::cutin=std::atof(tempstringvector[1].c_str());
		}
		else if(temp.find("cutout")!=std::string::npos){
			tempstringvector=split(temp,"=");
			filesys::cutoff=std::atof(tempstringvector[1].c_str());
		}
		else if(temp.find("listfile")!=std::string::npos){
			tempstringvector=split(temp,"=");
			filesys::atomlistfile=tempstringvector[1];
		}
		else if(temp.find("polar_direction_file")!=std::string::npos){
			tempstringvector=split(temp,"=");
			 filesys::polar_direction_file=new char [tempstringvector[1].length()+1];
			std::strcpy(filesys::polar_direction_file,tempstringvector[1].c_str());
		}
		else if(temp.find("polar_correlation_file")!=std::string::npos){
			tempstringvector=split(temp,"=");
			filesys::polar_correlation_file=new char [tempstringvector[1].length()+1];
			std::strcpy(filesys::polar_correlation_file,tempstringvector[1].c_str());
		}
		else if(temp.find("correlation_angle_file")!=std::string::npos){
			tempstringvector=split(temp,"=");
			filesys::correlation_angle_file=new char [tempstringvector[1].length()+1];
			std::strcpy(filesys::correlation_angle_file,tempstringvector[1].c_str());
		}
		else if(temp.find("correlation_distance_file")!=std::string::npos){
			tempstringvector=split(temp,"=");
			filesys::correlation_distance_file=new char [tempstringvector[1].length()+1];
			std::strcpy(filesys::correlation_distance_file,tempstringvector[1].c_str());
		}
	}
	fs.close();
	}
	else{
	}
	MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&filesys::cutin,1,MPI::DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&filesys::cutoff,1,MPI::DOUBLE,0,MPI_COMM_WORLD);
	int stringlength=0;
	stringlength=strlen(filesys::correlation_distance_file);
	MPI_Bcast(&stringlength,1,MPI::INT,0,MPI_COMM_WORLD);
	if(world_rank!=0){
		filesys::correlation_distance_file=new char [stringlength+1];
	}
	MPI_Bcast(filesys::correlation_distance_file,stringlength+1,MPI::CHAR,0,MPI_COMM_WORLD);
	stringlength=strlen(filesys::polar_direction_file);
	MPI_Bcast(&stringlength,1,MPI::INT,0,MPI_COMM_WORLD);
	if(world_rank!=0){
		filesys::polar_direction_file=new char [stringlength+1];
	}
	MPI_Bcast(filesys::polar_direction_file,stringlength+1,MPI::CHAR,0,MPI_COMM_WORLD);
	stringlength=strlen(filesys::polar_correlation_file);
	MPI_Bcast(&stringlength,1,MPI::INT,0,MPI_COMM_WORLD);
	if(world_rank!=0){
		filesys::polar_correlation_file=new char [stringlength+1];
	}
	MPI_Bcast(filesys::polar_correlation_file,stringlength+1,MPI::CHAR,0,MPI_COMM_WORLD);
	stringlength=strlen(filesys::correlation_angle_file);
	MPI_Bcast(&stringlength,1,MPI::INT,0,MPI_COMM_WORLD);
	if(world_rank!=0){
		filesys::correlation_angle_file=new char [stringlength+1];
	}
	MPI_Bcast(filesys::correlation_angle_file,stringlength+1,MPI::CHAR,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if(world_rank==0){
	fs.open(filesys::atomlistfile.c_str(),std::fstream::in);
	while(getline(fs,temp)){
		ss.clear();
		ss.str(temp);
		ss>>atomid;
		atomlist.push_back(atomid);
	}
	fs.close();
	atomlist.sort([](int a,int b)->bool{ return a<b;});
	/*searching the pairs that satifying the filesys::cutoff criterial*/
	for(std::list<int>::iterator a=atomlist.begin();a!=atomlist.end();a++)
		for(std::list<int>::iterator b=atomlist.begin();b!=atomlist.end();b++){
			pointA=changeindex(*a,cell);
			pointB=changeindex(*b,cell);
			distance=0.0;
			for(size_t i=0;i<3;i++){
				temp_double=pointA[i]-pointB[i]-round((pointA[i]-pointB[i]+0.0)/cell)*cell;
				distance=temp_double*temp_double+distance;
			}
			distance=sqrt(distance);
			if((*a < *b)&&( distance < filesys::cutoff /filesys::celllength )&&( distance > filesys::cutin /filesys::celllength )){
				pairlist.push_back(*a);
				pairlist.push_back(*b);
			}
			delete [] pointA;
			delete [] pointB;
		}
	}
	else{
	/*doing nothing and waiting for instructions*/
	}
	filesys::pairsize=pairlist.size()/2;
	MPI_Bcast(&filesys::pairsize,1,MPI::INT,0,MPI_COMM_WORLD);
	filesys::pairA=new int [filesys::pairsize];
	filesys::pairB=new int [filesys::pairsize];
	int count=0;
	for(std::list<int>::iterator a=pairlist.begin();a!=pairlist.end();a++){
		if(count%2==0){
			filesys::pairA[count/2]=*a;
		}
		else{
			filesys::pairB[(count-1)/2]=*a;
	  }
   count=count+1;
	}
	MPI_Bcast(filesys::pairA,filesys::pairsize,MPI::INT,0,MPI_COMM_WORLD);
	MPI_Bcast(filesys::pairB,filesys::pairsize,MPI::INT,0,MPI_COMM_WORLD);
	if(world_rank==1){
	std::cout<<"I have "<<filesys::pairsize<<" pairs"<<std::endl;
	std::cout<<filesys::cutin<<" "<<filesys::cutoff<<std::endl;
	printf("direction file %s\n",filesys::polar_direction_file);
	printf("correlation file %s\n",filesys::polar_correlation_file);
	printf("angle_file %s\n",filesys::correlation_angle_file);
	printf("distance_file %s\n",filesys::correlation_distance_file);
	}
}
