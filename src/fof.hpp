#ifndef _DTK_FOF_HPP_
#define _DTK_FOF_HPP_
#include <iostream>
#include "dtk/all.hpp"
#include "chainingmesh.hpp"

inline int find_last_link(int i, std::map<int,int>& recolor_map){
  while(recolor_map.count(i)!=0){
    i = recolor_map[i];
  }
  return i;
}
inline int find_last_link(int i, int* recolor_map){
  while(recolor_map[i] != i ){
    i = recolor_map[i];
  }
  return i;
}
inline void link(int i1,int i2,std::map<int,int>& recolor_map){
  
  if(i1>i2){
    std::swap(i1,i2); //make i1 always the smaller one;
  }
  //std::cout<<"linking "<<i2<<"->"<<i1;
  i1 = find_last_link(i1,recolor_map);
  //std::cout<<"->"<<i1<<std::endl;
  if(recolor_map.count(i2) != 0){
    int old_link = find_last_link(recolor_map[i2],recolor_map);
    //std::cout<<"\tnote: "<<i2<<"->"<<old_link<<std::endl;
    if(old_link > i1){
      //std::cout<<"\tlinked "<<old_link<<"->"<<i1<<std::endl;
      recolor_map[old_link] = i1;
    }
    else if(old_link < i1){
      //std::cout<<"\tlinked "<<i1<<"->"<<old_link<<std::endl;
      recolor_map[i1] = old_link;
      i1 = old_link;
    }
  }
  /*while(recolor_map.count(i1)!=0){//find the final color from ia->ib->ic-id to ia->id and update ib->id;
    i1 = recolor_map[i1]; //TODO update
    } */
  //std::cout<<"\tlinked "<<i2<<"->"<<i1<<std::endl;
  recolor_map[i2]=i1;
}
inline void link(int i1, int i2,int* color_map){
   if(i1>i2){
    std::swap(i1,i2); //make i1 always the smaller one;
  }
  i1 = find_last_link(i1,color_map);
  if(color_map[i2] != i2){
    int old_link = find_last_link(color_map[i2],color_map);
    if(old_link > i1){
      color_map[old_link] = i1;
    }
    else if(old_link < i1){
      color_map[i1] = old_link;
      i1 = old_link;
    }
  }
  color_map[i2]=i1;
}

template<typename T>
inline T dist(int i1, int i2, T* x, T*y, T*z){
  T dx = x[i1]-x[i2];
  T dy = y[i1]-y[i2];
  T dz = z[i1]-z[i2];
  return dx*dx + dy*dy + dz*dz;
}

template<typename T>
void interact_cells(T ll_sq, T* x, T* y, T* z, bool* use, 
		  int offset1, int size1, int offset2, int size2,
		  std::map<int,int>& recolor_map){
  for(int i1 = offset1;i1<(offset1+size1);++i1){
    if(use[i1])
      for(int i2=offset2;i2<(offset2+size2);++i2){
	if(use[i2] && i1!=i2){
	  if(dist(i1,i2,x,y,z) < ll_sq)
	    link(i1,i2,recolor_map);
	}
      }
  }
}
template<typename T>
void interact_cells(T ll_sq, T* x, T* y, T* z, bool* use, int* color,
		    int offset1, int size1, int offset2, int size2){
  for(int i1 = offset1;i1<(offset1+size1);++i1){
    if(use[i1])
      for(int i2=offset2;i2<(offset2+size2);++i2){
	if(use[i2] && i1!=i2){
	  if(dist(i1,i2,x,y,z) < ll_sq)
	    link(i1,i2,color);
	}
      }
  }
}

void assign_color(int* color, size_t size,std::map<int,int> recolor_map){
 //now the assigning a color(int) to each particle. Particles with the same color are in the
  //same fof group.
  for(int i =0;i<size;++i){
    //std::cout<<x[i]<<" "<<y[i]<<" "<<z[i]<<std::endl;
    if(recolor_map.count(i) ==0){
      color[i]=i; //no member;
    }
    else{
      //DEBUG
      if(recolor_map[i] == i){
	//std::cout<<"recolor issue"<<__LINE__<<std::endl;
	throw;
      }
      color[i] = color[recolor_map[i]]; //take the same color of the
      //particle it was linked. Since recolor_map[i] is always less than i, 
      //color[i] is correctly filled. 
    }
    //std::cout<<color[i]<<std::endl;
  }
}
void assign_color(int* color, size_t size){
 for(int i =0;i<size;++i){
      color[i] = color[color[i]]; //take the same color of the
  }
}


//chaining mesh method
template<typename T>
void find_fof_cm(T ll,T* x, T* y, T* z, bool* use,int* color, size_t size,ChainingMesh<T>& cm){
  find_fof_cm(ll,x, y,z, use, color, size,cm,1);
}
template<typename T>
void find_fof_cm(T ll,T* x, T* y, T* z, bool* use,int* color, size_t size,ChainingMesh<T>& cm,int grid_len){
  find_fof_cm(ll,x, y,z, use, color, size,cm,grid_len,grid_len,grid_len);
}
template<typename T>
void find_fof_cm(T ll,T* x, T* y, T* z, bool* use,int* color, size_t size,ChainingMesh<T>& cm,
	      int grid_len_i,int grid_len_j, int grid_len_k){
  //std::map<int,int> recolor_map;
  for(size_t i=0;i<size;++i){
    if(use[i])
      color[i] = i;
    else
      color[i]=-1; //group all the none used guys together. 
  }
  T ll_sq = ll*ll;
  int grid_i = cm.get_grid_i();
  int grid_j = cm.get_grid_j();
  int grid_k = cm.get_grid_k();
  for(int i=0;i<grid_i;++i){
    for(int j=0;j<grid_j;++j){
      for(int k=0;k<grid_k;++k){ //for each cell, interact it with itself & neighbors to one side.
	int base_offset,base_size; //the one side is so we don't get duplicate interactions
	cm.get_cell(i,j,k,base_offset,base_size);
	for(int i2=0;i2<grid_len_i;++i2){
	  for(int j2=0;j2<grid_len_j;++j2){
	    for(int k2=0;k2<grid_len_k;++k2){
	      int offset,size;
	      cm.get_cell(i+i2,j+j2,k+k2,offset,size);
	      // interact_cells(ll_sq,x,y,z,use,base_offset,base_size,offset,size,recolor_map);
	      interact_cells(ll_sq,x,y,z,use,color,base_offset,base_size,offset,size);
	    }
	  }
	}
      }
    }
  }
  assign_color(color,size);
}

//n2 with 
template<typename T>
void find_fof_n2(T ll,T* x, T* y, T* z, bool* use,int* color, size_t size){
  std::map<int,int> recolor_map;
  for(size_t i=0;i<size;++i)
    color[i] = i;
  T ll_sq = ll*ll;
  interact_cells(ll_sq,x,y,z,use,0,size,0,size,recolor_map);
  assign_color(color,size,recolor_map);
}


template<typename T>
void find_fof_srt(T ll,T* x, T* y, T* z,int* color, size_t size){
  for(int i=0;i<size;++i){
    color[i]=i;
  }
  float ll_sq = ll*ll;
  int* srt = dtk::new_sort_array<int>(size);
  std::sort(srt,srt+size,dtk::Sorter<T>(x));
  dtk::reorder(x,size,srt);
  dtk::reorder(y,size,srt);
  dtk::reorder(z,size,srt);
  delete [] srt;
  for(int i =0;i<size;++i){
    for(int j = i+1;(x[j]-x[i])<ll && j<size;++j){
      if(dist(i,j,x,y,z) < ll_sq){
	link(i,j,color);
      }
    }
  }
}

template<typename T>
void fof_group_colors(T* x, T* y, T* z, int* color, size_t size, 
		  size_t min_fof_size,
		  std::vector<std::vector<T> >& fof_x, 
		  std::vector<std::vector<T> >& fof_y, 
		  std::vector<std::vector<T> >& fof_z){
  std::vector<int> reordered_colors(color,color+size);
  int* srt = dtk::new_sort_array<int>(size);
  std::sort(srt,srt+size,dtk::Sorter<int>(color));
  dtk::reorder(reordered_colors,srt);
  std::vector<size_t> group_start,group_size;
  dtk::group_by_id(reordered_colors,group_start,group_size);
  for(int i =0;i<group_start.size();++i){
    int group_id = color[srt[group_start[i]]];
    if(group_id == -1) //unused cores
      continue;
    fof_x.push_back(std::vector<T>());
    fof_y.push_back(std::vector<T>());
    fof_z.push_back(std::vector<T>());
    for(int j =0;j<group_size[i];++j){
      fof_x[i].push_back(x[srt[group_start[i]+j]]);
      fof_y[i].push_back(y[srt[group_start[i]+j]]);
      fof_z[i].push_back(z[srt[group_start[i]+j]]);
    }
  }
}
#endif// _DTK_FOF_HPP_
