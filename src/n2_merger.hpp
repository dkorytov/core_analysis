#include <set>
#include <vector>
#include <stddef.h>
#include <iostream>
#include "chainingmesh.hpp"
extern "C" void n2_merger_float(float* x, float* y, int* w, int* size,float merg_len,int64_t* colors_out);
extern "C" void n2_merger_double(double* x, double* y, int* w, int* size,double merg_len,int64_t* colors_out);
extern "C" void n2_merger_float3d(float* x, float* y, float* z, int* w, int* size,float merg_len,int64_t* colors_out);
extern "C" void n2_merger_double3d(double* x, double* y, double* z,int* w, int* size,double merg_len,int64_t* colors_out);
				    


void merg_colors(int i , int j, std::vector<int64_t>& colors);

template<class T>
T average(std::vector<T> data){
  size_t size = data.size();
  T result = 0;
  for(size_t i =0;i<size;++i){
    result += data[i];
  }
  return result/static_cast<T>(size);
}
template<class T>
T median(std::vector<T> data){
  //std::
  return 0;
}

template<typename T>
void n2_merger(T* x,T* y, int* w, int* size,T merg_len,int64_t* colors_out){
  //std::cout<<"We are in the function!!"<<std::endl;
  //  std::cout<<"size: "<<size[0]<<" merg_len: "<<merg_len<<std::endl;
  //for(int i = 0;i<size[0];++i){
  //  std::cout<<x[i]<<" "<<y[i]<<std::endl;
  //}

  float merger_len;
  if(merg_len > 0)
    merger_len = merg_len*merg_len;
  else
    merger_len = 0;
  std::vector<int64_t> colors(size[0]);
  for(int i =0;i<size[0];++i){
    colors[i] =i;
  }
  //
  for(int i =0;i<size[0];++i){
      for(int j =0;j<i;++j){
	float dist_x = x[i]-x[j];
	float dist_y = y[i]-y[j];
	float dist = dist_x*dist_x+dist_y*dist_y;
	if(dist <= merger_len)
	  merg_colors(i,j,colors);
      }
  }
  //write out colors to the output if requested
  if(colors_out != NULL){
    for(int64_t i =0;i<colors.size();++i){
      colors_out[i]=colors[i];
    }
  }

  std::vector<float> avg_x;
  std::vector<float> avg_y;
  std::vector<float> avg_z;
  std::vector<int> avg_weight;

  std::set<int> set(colors.begin(),colors.end());
  std::vector<float> clr_x;
  std::vector<float> clr_y;
  std::vector<float> clr_z;

  //average all position of each color group
  for(std::set<int>::iterator it = set.begin();it!=set.end();++it){
    int c1 = *it;
    clr_x.clear();
    clr_y.clear();
    for(int i =0;i<size[0];++i){
      if(colors[i] ==c1){
	clr_x.push_back(x[i]);
	clr_y.push_back(y[i]);
      }
    }
    avg_x.push_back(average(clr_x));
    avg_y.push_back(average(clr_y));
    avg_weight.push_back(clr_x.size());
  }
  //overwrite the result to x,y,z,w& size
  for( int i =0;i<avg_x.size();++i){
    x[i]=avg_x[i];
    y[i]=avg_y[i];
    w[i]=avg_weight[i];
    //std::cout<<"xy: "<<x[i]<<" "<<y[i]<<std::endl;
  } 
  for(int i=avg_x.size();i<size[0];++i){//for debug
    x[i]=0.0;
    y[i]=0.0;
    w[i]=0.0;
  }
  size[0] = avg_x.size();
  return;
}

template<typename T>
void n2_merger3d(T* x,T* y, T* z, int* w, int* size,T merg_len,int64_t* colors_out){
  //std::cout<<"We are c n2_merger3d"<<std::endl;
  //std::cout<<"size: "<<size[0]<<" merg_len: "<<merg_len<<std::endl;
  //  for(int i = 0;i<10;++i){
  //  std::cout<<x[i]<<" "<<y[i]<<std::endl;
  // }
  float merger_len;
  if(merg_len > 0)
    merger_len = merg_len*merg_len;
  else
    merger_len = 0;
  std::vector<int64_t> colors(size[0]);
  for(int i =0;i<size[0];++i){
    colors[i] =i;
  }
  if(size[0]<=400000){
    //std::cout<<"n2 method"<<std::endl;
    for(int i =1;i<size[0];++i){
      for(int j =0;j<i;++j){
	float dist_x = x[i]-x[j];
	float dist_y = y[i]-y[j];
	float dist_z = z[i]-z[j];
	float dist = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
	if(dist <= merger_len)
	  merg_colors(i,j,colors);
      }
    }
  }
  else if(true){
    // dtk::Timer t1;t1.start();
    int* srt = dtk::arg_sort(x, size[0]);
    dtk::reorder(x,size[0],srt);
    dtk::reorder(y,size[0],srt);
    dtk::reorder(z,size[0],srt);
    dtk::reorder(w,size[0],srt);
    delete [] srt;
    float orginal_merger_len = std::sqrt(merger_len);
    for(int i =1;i<size[0];++i){
      for(int j =0;j<i;++j){
	float dist_x = x[i]-x[j];
	if(dist_x>orginal_merger_len)
	  break;
	float dist_y = y[i]-y[j];
	float dist_z = z[i]-z[j];
	float dist = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
	if(dist <= merger_len)
	  merg_colors(i,j,colors);
      }
    }
    // t1.stop();
    // dtk::Timer t2;t2.start();
    // for(int i =1;i<size[0];++i){
    //   for(int j =0;j<i;++j){
    // 	float dist_x = x[i]-x[j];
    // 	float dist_y = y[i]-y[j];
    // 	float dist_z = z[i]-z[j];
    // 	float dist = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
    // 	if(dist <= merger_len)
    // 	  merg_colors(i,j,colors);
    //   }
    // }
    // t2.stop();
    // std::cout<<"merging timing: n2: "<<t1<<"   n2 smart: "<<t2<<std::endl;
  }
  else{
    //if we have too many points a naive n2 would take too long. 
    //so we use a chaining mesh. 
    std::cout<<"cm method"<<std::endl;
    ChainingMesh<T> cmesh(256.0,256.0,256.0,16,16,16);
    cmesh.set_data(x,y,z,size[0]);
    size_t cell_num = cmesh.get_num_cells();
    //check the distance between all data points in each cell and 
    //neighbor cells. 
    for(size_t cell_i =0;cell_i<cell_num;++cell_i){
      std::vector<size_t> data_i = cmesh.get_points(cell_i,1,1,1);
      //std::cout<<"result size: "<<data_i.size()<<std::endl;
      for(int i=1;i<data_i.size();++i){
	for(int j =0;j<i;++j){
	  size_t i2 = data_i[i];
	  size_t j2 = data_i[j];
	  float dist_x = x[i2]-x[j2];
	  float dist_y = y[i2]-y[j2];
	  float dist_z = z[i2]-z[j2];
	  float dist = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
	  if(dist <= merger_len)
	    merg_colors(i2,j2,colors);
      }
    }
      
    }
  }
  //if an output is specified, write out fof color groups. 
  if(colors_out != NULL){
    for(int64_t i =0;i<colors.size();++i){
      colors_out[i]=colors[i];
    }
  }
  std::vector<float> avg_x;
  std::vector<float> avg_y;
  std::vector<float> avg_z;
  std::vector<int> avg_weight;

  std::set<int> set(colors.begin(),colors.end());
  std::vector<float> clr_x;
  std::vector<float> clr_y;
  std::vector<float> clr_z;
  //std::cout<<"calculating color mergers"<<std::endl;
  //average all position of each color group

  for(std::set<int>::iterator it = set.begin();it!=set.end();++it){
    int c1 = *it;
    clr_x.clear();
    clr_y.clear();
    clr_z.clear();
    for(int i =0;i<size[0];++i){
      if(colors[i] ==c1){
	clr_x.push_back(x[i]);
	clr_y.push_back(y[i]);
	clr_z.push_back(z[i]);

      }
    }
    //std::cout<<"\ti:"<<avg_x.size()<<"/"<<size[0]<<std::endl;
    avg_x.push_back(average(clr_x));
    avg_y.push_back(average(clr_y));
    avg_z.push_back(average(clr_z));
    avg_weight.push_back(clr_x.size());
  }
  //overwrite the result to x,y,z,w& size
  //  std::cout<<"getting pos"<<std::endl;

  //  std::cout<<"0"<<"->"<<avg_x.size()<<std::endl;
  for( int i =0;i<avg_x.size();++i){
    //std::cout<<"\ti:"<<i<<std::endl;
    x[i]=avg_x[i];
    y[i]=avg_y[i];
    z[i]=avg_z[i];
    w[i]=avg_weight[i];
  } 
  //std::cout<<"fill with zeros"<<std::endl;
  for(int i=avg_x.size();i<size[0];++i){//for debug
    //std::cout<<"\ti:"<<i<<std::endl;
    x[i]=0.0;
    y[i]=0.0;
    z[i]=0.0;
    w[i]=0.0;
  }
  size[0] = avg_x.size();
  return;
}



extern "C" void n2_merger3d_write_out(int size, 
				      float* x,  float* y,  float* z, 
				      float* vx, float* vy, float* vz,
				      int64_t* core_tag,
				      float* infall_time,
				      float* infall_mass,
				      float merg_len,
				      const char* file_loc);
