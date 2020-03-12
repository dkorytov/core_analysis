
#include <stdlib.h>
#include <iostream>
#include <unistd.h>
#include <numeric>

#include <assert.h>
#include <mpi.h>
#include <omp.h>
#include <H5Cpp.h>

#include <pthread.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "GenericIO.h"
#include "dtk/all.hpp"
#include "n2_merger.hpp"
#include "dtk/hdf5_util.hpp"
#include "chainingmesh.hpp"
#include "fof.hpp"

std::vector<float> x, y, z;
std::vector<int> w;
std::vector<int64_t> colors;
int size;

void load_cluster(char* fname){
  dtk::AutoTimer t;
  dtk::read_hdf5(fname, "x", x);
  dtk::read_hdf5(fname, "y", y);
  dtk::read_hdf5(fname, "z", z);
  w.resize(x.size());
  colors.resize(x.size());
  for(int i=0;i<w.size();++i){
    w[i] = i;
    colors[i] = i;
  }
  size = x.size();
  // int* srt = dtk::arg_sort(x.data(), size);
  // dtk::reorder(x.data(),size,srt);
  // dtk::reorder(y.data(),size,srt);
  // dtk::reorder(z.data(),size,srt);
  // dtk::reorder(w,size,srt);
  // std::cout<<"Load time: "<<t<<std::endl;
}
void merge(){
  std::cout<<"\n\n"<<std::endl;
  load_cluster("tmp_hdf5/single_cluster.hdf5");
  std::cout<<"N2 raw"<<std::endl;
  float distance = 0.1;
  dtk::AutoTimer t;
  n2_merger3d<float>(x.data(), y.data(), z.data(), w.data(), &size, distance, colors.data(), 40000000);
  std::cout<<"Merging time: "<<t<<std::endl;
  std::cout<<"size: "<<size<<std::endl;
  std::cout<<"\n\nN2 x-sort"<<std::endl;
  load_cluster("tmp_hdf5/single_cluster.hdf5");
  t.start();
  n2_merger3d<float>(x.data(), y.data(), z.data(), w.data(), &size, distance, colors.data(), 400);
  std::cout<<"Merging time: "<<t<<std::endl;
  std::cout<<"size: "<<size<<std::endl;

}
void merge_scaling(){
  load_cluster("tmp_hdf5/single_cluster.hdf5");
  for(int i =0;i<size;++i){
    std::cout<<i<<" ";
    dtk::Timer t;t.start();
    int current_size = i;
    float r = 0.02;
    n2_merger3d<float>(x.data(), y.data(), z.data(), w.data(), &current_size, r, colors.data(), i+1);
    t.stop();
    double time = t.get_useconds();
    std::cout<<time<<" ";
    t.start();
    n2_merger3d<float>(x.data(), y.data(), z.data(), w.data(), &current_size, r, colors.data(), 0);
    t.stop();
    time = t.get_useconds();
    std::cout<<time<<std::endl;
  }
}
int main(int argc, char** argv){
  // std::cout<<"Yo"<<std::endl;


  // merge();
  merge_scaling();
}
