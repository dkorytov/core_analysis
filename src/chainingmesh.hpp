#ifndef CHAININGMESH_HPP
#define CHAININGMESH_HPP
#include <vector>
#include <cstdlib>
#include "dtk/all.hpp"
#include <cmath>

template<typename T>
class ChainingMesh{
  T min_x,max_x;
  T min_y,max_y;
  T min_z,max_z;
  int grid_i;
  int grid_j;
  int grid_k;
  T grid_idx;
  T grid_jdy;
  T grid_kdz;
  T* data_x,*data_y,*data_z;
  int* cell_assignment;
  size_t* data_srt;
  size_t num_data;
  std::vector<size_t> cell_offset;
  std::vector<size_t> cell_size;
  bool periodic;
public:

  ChainingMesh():data_srt(NULL){
  }
  ChainingMesh(T max_x, T max_y, T mat_z,
	       int  grid_i, int  grid_j, int  grid_k):max_x(max_x),
						      max_y(max_y),
						      max_z(max_z),
						      grid_i(grid_i),
						      grid_j(grid_j),
						      grid_k(grid_k),
						      grid_idx(max_x/grid_i),
						      grid_jdy(max_y/grid_j),
						      grid_kdz(max_z/grid_k),
						      data_srt(NULL){
    cell_offset.resize(grid_i*grid_j*grid_k,0);
    cell_size.resize(grid_i*grid_j*grid_k,0);
  }
  ChainingMesh(T* x, T* y, T* z, size_t size,int grid_i, int grid_j, int grid_k,bool periodic=false):
    grid_i(grid_i), grid_j(grid_j), grid_k(grid_k),periodic(periodic),data_x(x),data_y(y),data_z(z),
    data_srt(NULL){
    max_x = dtk::max(x,size);
    min_x = dtk::min(x,size);
    max_y = dtk::max(y,size);
    min_y = dtk::min(y,size);
    max_z = dtk::max(z,size);
    min_z = dtk::min(z,size);
    //int* data_cell = new int[size];
  }
  
  ~ChainingMesh(){
    //if(cell_assignment != NULL)
      //delete [] cell_assignment;
  }
  void clear(){
    if(data_srt != NULL)
      delete [] data_srt;
    if(cell_assignment != NULL)
      delete [] cell_assignment;

  }
  int index1d_from_index3d(int i, int j, int k){
    if(periodic){
      handle_periodic(i,grid_i);
      handle_periodic(j,grid_j);
      handle_periodic(k,grid_k);
    }
    else if(i<0 || i>=grid_i || 
	    j<0 || j>=grid_j ||
	    k<0 || k>=grid_k){
      return -1; //outside the box
    }
    int res = i*grid_j*grid_k + j*grid_k + k;
    return res;
  }
  int index1d_from_pos(T x, T y, T z){

    if(periodic){
      handle_periodic(x,max_x);
      handle_periodic(y,max_y);
      handle_periodic(z,max_z);

    }
    int i = index_i(x);
    int j = index_j(y);
    int k = index_k(z);

    return index1d_from_index3d(i,j,k);
  }
  std::vector<int> index3d_from_pos(T x, T y, T z){
    std::vector<int> result(3);
    result[0]=index_i(x);
    result[1]=index_j(y);
    result[2]=index_k(z);
    return result;
  }
  std::vector<int> index3d_from_index1d(int index1d){
    std::vector<int> ijk(3);
    ijk[0]=index_i_from_index1d(index1d);
    ijk[1]=index_j_from_index1d(index1d);
    ijk[2]=index_k_from_index1d(index1d);
    return ijk;
  }
  int index_i_from_index1d(int index1d){
    return index1d/(grid_j*grid_k);
  }
  int index_j_from_index1d(int index1d){
    return (index1d/grid_k)%grid_j;
  }
  int index_k_from_index1d(int index1d){
    return index1d%grid_k;
  }
  
  void handle_periodic(int& i,int max_i){
    //i = std::mod(std::mod(i,max_i)+max_i,max_i);
    i = ((i%max_i)+max_i)%max_i;
  }
  void handle_periodic(T& x,T max_x){
    x = std::fmod(std::fmod(x,max_x)+max_x,max_x);
  }
  void set_limits(T* x, T* y, T* z, size_t size,int ni, int nj, int nk,
		  bool periodic=false,T buffer=.01){
    set_grids(ni,nj,nk);
    this->periodic = periodic;
    max_x = dtk::max(x,size);
    min_x = dtk::min(x,size);
    max_y = dtk::max(y,size);
    min_y = dtk::min(y,size);
    max_z = dtk::max(z,size);
    min_z = dtk::min(z,size);
    
    //make the boundaries a tiny bit outside the data limits
    if(!periodic){
      T del_x = max_x - min_x;
      T del_y = max_y - min_y;
      T del_z = max_z - min_z;

      max_x = max_x + buffer*del_x;
      min_x = min_x - buffer*del_x;
      max_y = max_y + buffer*del_y;
      min_y = min_y - buffer*del_y;
      max_z = max_z + buffer*del_z;
      min_z = min_z - buffer*del_z;
    }
    grid_idx = (max_x-min_x)/grid_i;
    grid_jdy = (max_y-min_y)/grid_j;
    grid_kdz = (max_z-min_z)/grid_k;
    cell_offset.resize(grid_i*grid_j*grid_k,0);
    cell_size.resize(grid_i*grid_j*grid_k,0);
  }
  void set_limits(T min_x,T max_x,T min_y,T max_y,T min_z,T max_z,
		  int grid_ni, int grid_nj, int grid_nk, bool periodic){
    this->min_x = min_x;
    this->max_x = max_x;
    this->min_y = min_y;
    this->max_y = max_y;
    this->min_z = min_z;
    this->max_z = max_z;
    set_grids(grid_ni,grid_nj, grid_nk);
    this->periodic = periodic;
  }
  void set_grids(int ni,int nj, int nk){
    grid_i = ni;
    grid_j = nj;
    grid_k = nk;
    grid_idx = (max_x-min_x)/grid_i;
    grid_jdy = (max_y-min_y)/grid_j;
    grid_kdz = (max_z-min_z)/grid_k;
    cell_offset.resize(ni*nj*nk);
    cell_size.resize(ni*nj*nk);
  }
  void set_data(T* x, T* y, T* z,size_t size){
    //first of all, make sure all the data is inside the period boundary
    num_data = size;
    int size_i = size;
    cell_assignment = new int[size];
    data_srt = dtk::new_sort_array<size_t>(size);
    for(size_t i =0;i<size;++i){
      cell_assignment[i] = index1d_from_pos(x[i],y[i],z[i]); //assign which cell each datapoint should go
    }
    //reoder the data so data points in the same cell are adjacent;
    std::sort(data_srt,data_srt+size,dtk::Sorter<int>(cell_assignment));
    std::vector<size_t> start;
    std::vector<size_t> len;
    //copy the cell assignment and reorder them into ascending order
    std::vector<int> ordered_assign(cell_assignment,cell_assignment+size); 
    dtk::reorder(ordered_assign,data_srt);
    dtk::group_by_id(ordered_assign,start,len); //find the size of each group of cells
    size_t running_offset=0;
    size_t last_group=0;
    //we must assign each cell an offset and size;
    //this offset will point where the cell starts in the srt array. 
    for(size_t i =0;i<cell_offset.size();++i){
      if(last_group >= start.size() || ordered_assign[start[last_group]] < i){
	cell_offset[i]= running_offset;
	cell_size[i]  = 0; //no data points in this cell;
      }
      else if(ordered_assign[start[last_group]] == i){
	cell_offset[i]=start[last_group];
	cell_size[i]  =len[last_group];
	running_offset = cell_offset[i]; //empty cells will use this value 
	++last_group;                //to point in between filled cells. 
      }
    }
  }
  void print_grids(){
    std::cout<<grid_i<<" "<<grid_j<<" "<<grid_k<<std::endl;
    std::cout<<grid_idx<<" "<<grid_jdy<<" "<<grid_kdz<<std::endl;
  }
  std::vector<size_t> get_points(T x, T y, T z, T r){
    int buffer_i = r/grid_i;
    int buffer_j = r/grid_j;
    int buffer_k = r/grid_k;
    int cell_indx = index(x,y,z);
    return get_points(cell_indx,buffer_i,buffer_j,buffer_k);
  }
  
  int index_i(T x){
    return (int)((x-min_x)/grid_idx);
  }
  int index_j(T y){
    return (int)((y-min_y)/grid_jdy);
  }
  int index_k(T z){
    return (int)((z-min_z)/grid_kdz);
  }
  size_t get_num_cells(){
    return cell_offset.size();
  }
  const size_t* get_srt(){
    return data_srt;
  }
  const int* get_cell_assignment(){
    return cell_assignment;
  }
  std::vector<size_t> get_points(int cell_indx, int buffer_i,int buffer_j, int buffer_k){
    std::vector<size_t> result;
    std::vector<size_t> cells;
    int i = index_i(cell_indx);
    int j = index_j(cell_indx);
    int k = index_k(cell_indx);
    //std::cout<<"cell_indx"<<cell_indx<<std::endl;
    //std::cout<<"i:"<<i<<" j:"<<j<<" k:"<<k<<std::endl;
    //std::cout<<"i:"<<grid_i<<"j:"<<grid_j<<"k:"<<grid_k<<std::endl;
    for(int di =-buffer_i;di<=buffer_i;++di){
      for(int dj =-buffer_j;dj<=buffer_j;++dj){
	for(int dk =-buffer_k;dk<=buffer_k;++dk){
	  int ii = i+di;
	  int jj = j+dj;
	  int kk = k+dk;
	  //std::cout<<"*"<<ii<<" "<<jj<<" "<<kk<<std::endl;
	  int cell_indx = index1d_from_index3d(ii,jj,kk);
	  //std::cout<<"\t"<<cell_indx<<std::endl;
	  //std::cout<<"cell size:"<<cell_size[cell_indx]<<std::endl;
	  for(size_t cell_i=0; cell_i<cell_size[cell_indx]; ++cell_i){
	    //std::cout<<"i:"<<cell_i<<"/"<<cell_size[cell_indx]<<std::endl;
	    size_t tmp_result = data_srt[cell_offset[cell_indx]+cell_i];
	    //std::cout<<"\t"<<tmp_result<<std::endl;
	    result.push_back(tmp_result);
	  }
	}
      }
    }
    //std::cout<<"Done"<<std::endl;
    return result;
  }
  
  void get_cell(int index1d,int& offset,int& size){
    if(index1d == -1){ //an outside bounds index;
      offset = 0; //size ==0 
      size = 0;
    }
    else{
      offset = cell_offset.at(index1d);
      size   = cell_size.at(index1d);
    }
  }
  void get_cell(int i,int j, int k, int& offset, int& size){
    get_cell(index1d_from_index3d(i,j,k),offset,size);
  }
  int get_grid_i(){
    return grid_i;
  }
  int get_grid_j(){
    return grid_j;
  }
  int get_grid_k(){
    return grid_k;
  }

};
#endif// CHAININGMESH_HPP
