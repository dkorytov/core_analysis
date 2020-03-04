
#include "n2_merger.hpp"

void merg_colors(int i, int j,std::vector<int64_t>& colors){
  int c1 = colors[i];
  int c2 = colors[j];
  if(c1 == c2) //no need to merg two of the same colors
    return;
  for(int k=0;k<i;++k){ //take any c2 color and make it c1
    if(colors[k]==c2)
      colors[k]=c1;
  }
}
void merg_colors_fast(int i, int j, std::vector<int64_t>& colors){
  if(colors[i]==i){
    colors[i] = j;
  }
  throw; //not implemented
}
void merg_colors_fast_post_processing(std::vector<int64_t>& colors){
  for(int i =1;i<colors.size();++i){
    colors[i]=colors[colors[i]];
  }
  throw; // not implemented
}
void n2_merger_float(float* x, float* y, int* w, int* size,float merg_len,int64_t* colors_out){
  n2_merger<float>(x,y,w,size,merg_len,colors_out);
}
void n2_merger_double(double* x, double* y, int* w, int* size,double merg_len,int64_t* colors_out){
  n2_merger<double>(x,y,w,size,merg_len,colors_out);
}
void n2_merger_float3d(float* x, float* y, float* z, int* w, int* size,float merg_len,int64_t* colors_out){
  n2_merger3d<float>(x,y,z,w,size,merg_len,colors_out);
}
void n2_merger_double3d(double* x, double* y, double* z, int* w, int* size,double merg_len,int64_t* colors_out){
  n2_merger3d<double>(x,y,z,w,size,merg_len,colors_out);
}
void n2_merger3d_write_out(int size, 
			   float* x,  float* y,  float* z, 
			   float* vx, float* vy, float* vz,
			   int64_t* core_tag,
			   float* infall_time,
			   float* infall_mass,
			   float merg_len,
			   const char* file_loc){
  
  std::cout<<"starting with size:"<<size<<std::endl;
  float merger_len;
  if(merg_len > 0)
    merger_len = merg_len*merg_len;
  else
    merger_len = 0;
  std::vector<int64_t> colors(size);
  for(int i =0;i<size;++i){
    colors[i] =i;
  }
  if(size<=400000){
    std::cout<<"n2 method"<<std::endl;
    for(int i =1;i<size;++i){
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
  else{
    //if we have too many points a naive n2 would take too long. 
    //so we use a chaining mesh. 
    std::cout<<"cm method"<<std::endl;
    ChainingMesh<float> cmesh(256.0,256.0,256.0,16,16,16);
    cmesh.set_data(x,y,z,size);
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
  int* colors_out = new int[size];
  int* w = new int[size];
  for(int64_t i =0;i<colors.size();++i){
    colors_out[i]=colors[i];
  }

  std::vector<float> avg_x;
  std::vector<float> avg_y;
  std::vector<float> avg_z;
  std::vector<int> avg_weight;


  std::vector<float> avg_vx;
  std::vector<float> avg_vy;
  std::vector<float> avg_vz;
  std::vector<float> med_vx;
  std::vector<float> med_vy;
  std::vector<float> med_vz;
  std::vector<int64_t> min_core_tag;
  std::vector<float> max_infall_mass;
  std::vector<float> min_infall_mass;
  std::vector<float> avg_infall_mass;
  std::vector<float> med_infall_mass;
  
  std::set<int> set(colors.begin(),colors.end());
  std::vector<float> clr_x;
  std::vector<float> clr_y;
  std::vector<float> clr_z;

  std::vector<float> clr_vx;
  std::vector<float> clr_vy;
  std::vector<float> clr_vz;
  std::vector<int64_t> clr_core_tags;
  std::vector<float> clr_infall_mass;
  std::vector<float> clr_infall_time;

  //std::cout<<"calculating color mergers"<<std::endl;
  //average all position of each color group

  for(std::set<int>::iterator it = set.begin();it!=set.end();++it){
    int c1 = *it;
    clr_x.clear();
    clr_y.clear();
    clr_z.clear();
    clr_vx.clear();
    clr_vy.clear();
    clr_vz.clear();
    clr_core_tags.clear();
    clr_infall_mass.clear();
    clr_infall_time.clear();
    
    for(int i =0;i<size;++i){
      if(colors[i] ==c1){
	clr_x.push_back(x[i]);
	clr_y.push_back(y[i]);
	clr_z.push_back(z[i]);
	clr_vx.push_back(vx[i]);
	clr_vy.push_back(vy[i]);
	clr_vz.push_back(vz[i]);
  	clr_core_tags.push_back(core_tag[i]);
	clr_infall_mass.push_back(infall_mass[i]);
	clr_infall_time.push_back(infall_time[i]);
      }
    }
    //std::cout<<"\ti:"<<avg_x.size()<<"/"<<size[0]<<std::endl;
    avg_x.push_back(average(clr_x));
    avg_y.push_back(average(clr_y));
    avg_z.push_back(average(clr_z));
    avg_vx.push_back(average(clr_vx));
    avg_vy.push_back(average(clr_vy));
    avg_vz.push_back(average(clr_vz));
    med_vx.push_back(median(clr_vx));
    med_vy.push_back(median(clr_vy));
    med_vz.push_back(median(clr_vz));
    min_core_tag.push_back(dtk::min(clr_core_tags));
    max_infall_mass.push_back(dtk::max(clr_infall_mass));
    min_infall_mass.push_back(dtk::min(clr_infall_mass));
    med_infall_mass.push_back(median(clr_infall_mass));
    avg_weight.push_back(clr_x.size());
  }
  //overwrite the result to x,y,z,w& size
  //  std::cout<<"getting pos"<<std::endl;

  //  std::cout<<"0"<<"->"<<avg_x.size()<<std::endl;

  dtk::write_binary("processed_catalog/avg_x.bin",avg_x);
  dtk::write_binary("processed_catalog/avg_y.bin",avg_y);
  dtk::write_binary("processed_catalog/avg_z.bin",avg_z);
  dtk::write_binary("processed_catalog/avg_vx.bin",avg_vx);
  dtk::write_binary("processed_catalog/avg_vy.bin",avg_vy);
  dtk::write_binary("processed_catalog/avg_vz.bin",avg_vz);
  dtk::write_binary("processed_catalog/med_vx.bin",med_vx);
  dtk::write_binary("processed_catalog/med_vy.bin",med_vy);
  dtk::write_binary("processed_catalog/med_vz.bin",med_vz);
  dtk::write_binary("processed_catalog/min_core_tag.bin",min_core_tag);
  dtk::write_binary("processed_catalog/min_infall_mass.bin",min_infall_mass);
  dtk::write_binary("processed_catalog/median_infall_mass.bin",med_infall_mass);
  dtk::write_binary("processed_catalog/max_infall_mass.bin",max_infall_mass);
  std::cout<<"writing to text file"<<std::endl;
  std::ofstream file(file_loc);

  file<<"avg_x\tavg_y\tavg_z\tavg_vx\tavg_vy\tavg_vz\tmed_vx\tmed_vy\tmed_vz\tcore_tag\tmin_mass\tmax_mass\tmedian_mass"<<std::endl;
  for(size_t i =0;i<avg_x.size();++i){
    file<<avg_x.at(i)<<"\t"
	<<avg_y.at(i)<<"\t"
	<<avg_z.at(i)<<"\t"
	<<avg_vx.at(i)<<"\t"
	<<avg_vy.at(i)<<"\t"
	<<avg_vz.at(i)<<"\t"
	<<med_vx.at(i)<<"\t"
	<<med_vy.at(i)<<"\t"
	<<med_vz.at(i)<<"\t"
	<<min_core_tag.at(i)<<"\t"
	<<min_infall_mass.at(i)<<"\t"
	<<max_infall_mass.at(i)<<"\t"
	<<med_infall_mass.at(i)<<"\t"
	<<std::endl;
  }
  file.close();
  std::cout<<"done writing"<<std::endl;
  delete [] w;
  delete [] colors_out;
  return;
  

}
