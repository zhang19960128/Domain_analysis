#ifndef constant_h
#define constant_h
#include <iostream>
#include <string>
namespace filesys{
extern double celllength;/*cell length*/
extern double cutin;
extern double cutoff;/*angstrom*/
extern char* polar_direction_file;
extern char* polar_correlation_file;
extern char* correlation_angle_file;
extern char* correlation_distance_file;
extern std::string atomlistfile;
extern int* pairA;
extern int* pairB;
extern int pairsize;
}
#endif
