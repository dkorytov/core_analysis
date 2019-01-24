
#include <stdlib.h>
#include <iostream>
#include <unistd.h>
#include <numeric>

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



struct ZMR{
  //bin edges
  std::vector<float> z_bins; 
  std::vector<float> m_bins;
  std::vector<float> r_bins;
  //bin averages
  std::vector<float> z_bin_avg;
  std::vector<float> m_bin_avg;
  std::vector<float> r_bin_avg;
  
  std::vector<float> r_bin_area;
  std::vector<float> r_bin_length;
  int z_size;
  int m_size;
  int r_size;
  int zm_size;   // z_bins.size()*m_bins.size()
  int mr_size;
  int zmr_size;  // z_bins.size()*m_bins.size()*r_bins.size()
  //Ngal within r200
  std::vector<float> zm_Ngal;
  std::vector<float> zm_Ngal_err;
  std::vector<float> zm_Ngal_var;
  
  std::vector<float> zmr_gal_counts; //the raw count of how many galaxies in each bin
  //Galaxy Denisty per unit r200^2 in each radial bin
  std::vector<float> zmr_gal_density;
  std::vector<float> zmr_gal_density_err;
  std::vector<float> zmr_gal_density_var;
  //Galaxy integrate count (last radial bin = Ngal)
  std::vector<float> zmr_gal_accum;
  std::vector<float> zmr_gal_accum_err;
  std::vector<float> zmr_gal_accum_var;
  //Galaxy dn/dr 
  std::vector<float> zmr_dgal_dr;
  std::vector<float> zmr_dgal_dr_err;
  std::vector<float> zmr_dgal_dr_var;
  //Number of clusters by z,m bins
  std::vector<int64_t> zm_counts;
  //Number of clusters by z,m,r bins
  //not really used since there aren't radial cuts
  std::vector<int64_t> zmr_counts;
  int centrals_found;
  int centrals_expected;
  void zero_out(){
    for(int i =0;i<zm_size;++i){
      zm_Ngal_var[i]=0;
      zm_Ngal_err[i]=0;
      zm_Ngal[i]=0;
      zm_counts[i] = 0;
    }
    for(int i =0;i<zmr_size;++i){
      int j =0;
      zmr_gal_counts[i]=0;
      zmr_gal_density[i]=0.0;    
      zmr_gal_density_err[i]=0.0;
      zmr_gal_density_var[i]=0.0;
      zmr_gal_accum[i]=0.0;    
      zmr_gal_accum_err[i]=0.0;
      zmr_gal_accum_var[i]=0.0;
      zmr_dgal_dr[i]=0.0;
      zmr_dgal_dr_err[i]=0.0;
      zmr_dgal_dr_var[i]=0.0;
      zmr_counts[i]=0.0;      
    }
    centrals_found = 0;
    centrals_expected = 0;
  }
  void clear(){
    if(z_size !=0){
      zm_size = 0;
      zmr_size = 0;
      /*delete [] zm_Ngal;
      delete [] zm_Ngal_err;
      delete [] zm_Ngal_var;
      delete [] zmr_gal_counts;
      delete [] zmr_gal_density;
      delete [] zmr_gal_density_err;
      delete [] zmr_gal_density_var;
      delete [] zmr_gal_accum;
      delete [] zmr_gal_accum_err;
      delete [] zmr_gal_accum_var;
      delete [] zmr_dgal_dr;
      delete [] zmr_dgal_dr_err;
      delete [] zmr_dgal_dr_var;
      delete [] zm_counts;
      delete [] zmr_counts;*/
    }
  }
  void alloc(int z_size,int m_size,int r_size){
    this->z_size = z_size;
    this->m_size = m_size;
    this->r_size = r_size;
    zm_size = z_size*m_size;
    mr_size = m_size*r_size;
    zmr_size = z_size*m_size*r_size;
    zm_Ngal.resize(zm_size);
    zm_Ngal_err .resize(zm_size);
    zm_Ngal_var .resize(zm_size);
    zmr_gal_counts     .resize(zmr_size);
    zmr_gal_density    .resize(zmr_size);
    zmr_gal_density_err.resize(zmr_size);
    zmr_gal_density_var.resize(zmr_size);
    zmr_gal_accum      .resize(zmr_size);
    zmr_gal_accum_err  .resize(zmr_size);
    zmr_gal_accum_var  .resize(zmr_size);
    zmr_dgal_dr        .resize(zmr_size);
    zmr_dgal_dr_err    .resize(zmr_size);
    zmr_dgal_dr_var    .resize(zmr_size);
    zm_counts.resize(zm_size);
    zmr_counts.resize(zmr_size);
    zero_out();
  }
  int index(int z_i,int m_i) const{
    return z_i*m_size + m_i;
  }
  int index(int z_i,int m_i,int r_i) const{
    return z_i*mr_size + m_i*r_size + r_i;
  }
  ZMR(){}
  virtual ~ZMR(){
    clear();
  }
  void load(std::string zmr_file_name){
    H5::H5File zmr_file(zmr_file_name,H5F_ACC_RDONLY);
    dtk::read_hdf5(zmr_file,"z_bins",z_bins);    
    dtk::read_hdf5(zmr_file,"m_bins",m_bins);    
    dtk::read_hdf5(zmr_file,"r_bins",r_bins);    
    alloc(z_bins.size()-1,m_bins.size()-1,r_bins.size()-1);    
    dtk::read_hdf5(zmr_file,"zm_Ngal"    ,zm_Ngal);    
    dtk::read_hdf5(zmr_file,"zm_Ngal_var",zm_Ngal_var);
    dtk::read_hdf5(zmr_file,"zm_Ngal_err",zm_Ngal_err);
    dtk::read_hdf5(zmr_file,"zmr_gal_density",zmr_gal_density);
    dtk::read_hdf5(zmr_file,"zmr_gal_density_var",zmr_gal_density_var);
    dtk::read_hdf5(zmr_file,"zmr_gal_density_err",zmr_gal_density_err);
    dtk::read_hdf5(zmr_file,"zm_counts", zm_counts);
    dtk::read_hdf5(zmr_file,"zmr_counts",zmr_counts);
    z_bin_avg.reserve(z_size-1);
    m_bin_avg.reserve(m_size-1);
    r_bin_avg.reserve(r_size-1);
    r_bin_area.reserve(r_size-1);
    for(int i=0;i<z_size;++i)
      z_bin_avg.push_back((z_bins[i]+z_bins[i+1])/2.0);
    for(int i=0;i<m_size;++i)
      m_bin_avg.push_back((m_bins[i]+m_bins[i+1])/2.0);
    for(int i=0;i<r_size;++i){
      r_bin_avg.push_back((r_bins[i]+r_bins[i+1])/2.0);
      r_bin_area.push_back(M_PI*(r_bins[i+1]*r_bins[i+1] - r_bins[i]*r_bins[i]));
      r_bin_length.push_back(r_bins[i+1]-r_bins[i]);

    }
    //making the dgal/dr profiles
    for(int zi=0;zi<z_size;++zi){
      for(int mi=0;mi<m_size;++mi){
	float gal_accum=0;
	std::vector<float> gal_accum_err;
	std::vector<float> gal_accum_var;
	for(int ri=0;ri<r_size;++ri){
	  int indx = index(zi,mi,ri);
	  zmr_dgal_dr[indx] = zmr_gal_density[indx]*r_bin_area[ri]/r_bin_length[ri];
	  zmr_dgal_dr_err[indx] = zmr_gal_density_err[indx]*r_bin_area[ri]/r_bin_length[ri];
	  zmr_dgal_dr_var[indx] = zmr_gal_density_var[indx]*r_bin_area[ri]/r_bin_length[ri];
	  gal_accum += zmr_gal_density[indx]*r_bin_area[ri];
	  zmr_gal_accum[indx]=gal_accum;
	  zmr_gal_accum_err[indx]=sqrt(gal_accum*zmr_counts[indx])/zmr_counts[indx];
	  zmr_gal_accum_var[indx]=0.0;
	  zmr_gal_counts[indx] = zmr_gal_density[indx]*zmr_counts[indx]*r_bin_area[ri];
	  //gal_accum_err.push_back(zmr_dgal_dr_err[indx]);
	  //gal_accum_var.push_back(zmr_dgal_dr_var[indx]);
	  //zmr_gal_accum_err[indx]=dtk::root_mean_squared(gal_accum_err);
	  //zmr_gal_accum_var[indx]=dtk::root_mean_squared(gal_accum_var);
	}
      }
    }
  }
  void copy_bins(ZMR& ref){
    z_bins = ref.z_bins;
    m_bins = ref.m_bins;
    r_bins = ref.r_bins;
    z_bin_avg = ref.z_bin_avg;
    m_bin_avg = ref.m_bin_avg;
    r_bin_avg = ref.r_bin_avg;
    //r_bin_area = ref.r_bin_area;
    //r_bin_length  = ref.r_bin_length;
    z_size = ref.z_size;
    m_size = ref.m_size;
    r_size = ref.r_size;
    r_bin_area.resize(r_size);
    r_bin_length.resize(r_size);
    for(int i=0;i<r_size;++i){
      r_bin_area.at(i) = M_PI*(r_bins[i+1]*r_bins[i+1]-r_bins[i]*r_bins[i]);
      r_bin_length.at(i) = r_bins[i+1]-r_bins[i];
    }
    alloc(z_size,m_size,r_size);
  }
  void add_cluster(int z_i,int m_i,std::vector<float>& r_cnts,int Ngal){
    // std::cout<<"adding cluster: "<<z_i<<" "<<m_i<<std::endl;
    if(z_i == -1 || m_i == -1)
      return;
    int indx = index(z_i,m_i);
    zm_Ngal.at(indx)+=Ngal;
    zm_counts.at(indx)+=1;
    for(int r_i =0;r_i<r_size;++r_i){
      int indx = index(z_i,m_i,r_i);
      zmr_gal_counts.at(indx) +=r_cnts.at(r_i);
      // std::cout<<r_i<<"/"<<r_size<<" "<<r_cnts.at(r_i)<<std::endl;
      zmr_counts.at(indx)+=1;
    }
    //counting expected centrals
    if(r_bins[0] == 0 ){//if we have the central cluster area in the analysis
      ++centrals_expected;
      if(r_cnts[0]>0)
	++centrals_found;
    }
  }
  void finalize_data(){
    for(int i =0;i<zm_size;++i){
      zm_Ngal_var[i] = 0.0;
      zm_Ngal_err[i] = sqrt(zm_Ngal[i])/zm_counts[i];
      zm_Ngal[i]     = zm_Ngal[i]/zm_counts[i];
    }
    for(int z_i=0;z_i<z_size;++z_i){
      for(int m_i=0;m_i<m_size;++m_i){
	float gal_accum =0;
	for(int r_i=0;r_i<r_size;++r_i){
	  int indx = index(z_i,m_i,r_i);
	  zmr_gal_density[indx] = zmr_gal_counts[indx]/zmr_counts[indx]/r_bin_area.at(r_i);
	  zmr_gal_density_var[indx] = 0; 
	  zmr_gal_density_err[indx] = sqrt(zmr_gal_counts[indx])/zmr_counts[indx]/r_bin_area.at(r_i);
	  zmr_dgal_dr[indx] = zmr_gal_counts[indx]/zmr_counts[indx]/r_bin_length.at(r_i);
	  zmr_dgal_dr_err[indx] = sqrt(zmr_gal_counts[indx])/zmr_counts[indx]/r_bin_length.at(r_i);
	  zmr_dgal_dr_var[indx] = 0;
	  gal_accum += zmr_gal_counts[indx];
	  zmr_gal_accum[indx] = gal_accum/zmr_counts[indx];
	  zmr_gal_accum_err[indx] = sqrt(gal_accum)/zmr_counts[indx];
	  zmr_gal_accum_var[indx] =0;
	  /*std::cout<<z_i<<" "<<m_i<<" "<<r_i<<": "<<zmr_gal_counts[indx]
		   <<" "<<zmr_gal_density[indx]
		   <<" "<<zmr_dgal_dr[indx]
		   <<" "<<zmr_gal_accum[indx]<<std::endl;*/
	}
      }
    }
    //dtk::pause();
  }
  bool check_same_bins(const ZMR& ref)const{
    if(ref.z_bins.size() != z_bins.size()){
      std::cout<<"diff z bins"<<std::endl;
      return false;
    }
    if(ref.m_bins.size() != m_bins.size()){
      std::cout<<"diff m bins"<<std::endl;
      return false;
    }
    if(ref.r_bins.size() != r_bins.size()){
      std::cout<<"diff r bins"<<std::endl;

      return false;
    }
    for(int i =0;i<z_bins.size();++i)
      if(ref.z_bins[i]!=z_bins[i]){
	std::cout<<"diff z["<<i<<"] bins"<<std::endl;
	std::cout<<ref.z_bins[i]<<" "<<z_bins[i]<<std::endl;
	return false;
      }
    for(int i =0;i<m_bins.size();++i)
      if(ref.m_bins[i]!=m_bins[i]){
	std::cout<<"diff m["<<i<<"] bins"<<std::endl;
	return false;
      }
    for(int i =0;i<r_bins.size();++i)
      if(ref.r_bins[i]!=r_bins[i]){
	std::cout<<"diff r["<<i<<"] bins"<<std::endl;
	for(int j =0;j<r_bins.size();++j){
	  std::cout<<ref.r_bins[j]<<" "<<r_bins[j]<<std::endl;
	}
	for(int j =0;j<z_bins.size();++j){
	  std::cout<<ref.z_bins[j]<<" "<<z_bins[j]<<std::endl;
	}

	for(int j =0;j<m_bins.size();++j){
	  std::cout<<ref.m_bins[j]<<" "<<m_bins[j]<<std::endl;
	}
	throw;
	return false;
      }
    //if the bins are the same, we are all good :)
    return true;
  }
  void write_txt(std::string);
  void print_non_zero(){
    std::cout<<"Non-zero gal_density"<<std::endl;
    for(int zi=0;zi<z_size;++zi){
      for(int mi=0;mi<m_size;++mi){
	for(int ri=0;ri<r_size;++ri){
	  int indx = index(zi,mi,ri);
	  float val = zmr_gal_density[indx];
	  if(val != 0)
	    std::cout<<"["<<zi<<","<<mi<<","<<ri<<"]: "<<val<<std::endl;
	  }
	}
      }
  }

};


struct CoreParam{ 
  //all variables are NOT log
  //values
  float m_infall;
  float r_disrupt;
  float r_fof;
  float r_merger;
};

struct HaloCat{
  int64_t size;
  int64_t* htag;
  //float*   fof_mass;
  float*   sod_mass;
  float*   x;
  float*   y;
  float*   z;
  float*   sod_radius;
};
struct Cores{
  int64_t  size;
  int64_t* ctag;
  int64_t* host_htag; //current host htag
  int64_t* infall_htag; //infall htag
  int *    infall_step; //infall step
  float*   infall_mass;
  int64_t* central_htag;
  int*     central_step;
  float*   central_mass;
  int*     is_central;
  float* radius;
  float* x;
  float* y;
  float* z;
  int* compact;
};
struct CorePrtcls{
  int64_t size;
  int64_t* htag;
  int64_t* pid;
};
struct AccumPrtcls{
  int64_t size;
  int64_t* pid;
  float* x;
  float* y;
  float* z;
};

struct Galaxies{
  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> z;
  std::vector<float> w; //the number of mergered galaxies.
  std::vector<float> r; //distance from the center of the cluster.
  std::vector<int>   type;
  void reserve(size_t size){
    x.reserve(size);
    y.reserve(size);
    z.reserve(size);
    w.reserve(size);
    r.reserve(size);
    type.reserve(size);
  }
  void print(int num){
    std::cout<<"Galaxies: "<<x.size()<<std::endl;
    for(int i=0;i<num;++i){
      std::cout<<"["<<i<<"]"<<x[i]<<" "<<y[i]<<" "<<z[i]<<" "<<r[i]<<std::endl;
    }
  }
};
struct Cluster{
  int64_t htag;
  float sod_mass;
  float sod_radius;
  float x,y,z;
  float redshift;
  int step;
  int z_i,m_i; //the redshift and mass bin this cluster belongs in. 
  //the core catalog that belong to this cluster
  int64_t core_size;
  int64_t* core_id,*core_htag;
  float*  core_x,*core_y,*core_z,*core_r,*core_m;
  int*    core_is_central;
  int*    core_step;
  int*    core_cp_offset;
  int*    core_cp_size;
  //the core particles that belong to this cluster
  int64_t cp_size;
  float*  cp_x,*cp_y,*cp_z;
  bool*   cp_usage;
  int*    cp_color;
  ChainingMesh<float> cp_cm;
  
  std::vector<float> bcg_x,bcg_y,bcg_z,bcg_n;
  //tmp debug arrays
  std::vector<float> dis_cp_x,dis_cp_y,dis_cp_z;
  Cluster(){}
  Cluster(int64_t htag,float sod_mass,float rad,float x,float y, float z,int step);
  void get_compact_galaxies(float m_infall, float r_disrupt, Galaxies& gal) const;
  void get_compact_merged_galaxies(float m_infall,float r_disrupt, float r_merger, Galaxies& gal, bool verbose = false)const;
  void get_fof_galaxies(float m_infall, float r_disrupt, float r_fof,Galaxies& gal) const;
  void get_fof_galaxies2(float m_infall,float r_disrupt, float r_fof, Galaxies& gal)const;
  void get_fof_only_galaxies(float r_disrupt,Galaxies& gal);

  void get_radial_bins(float m_infall,float r_disrupt,float r_fof, 
		       std::vector<float>& r2_bins, float Ngal_r2_lim,
		       std::vector<float>& r_cnt, float& Ngal) const;
  void get_radial_bins(Galaxies& gal,std::vector<float>& r_bins,
		       std::vector<float>& r2_bins, float Ngal_r2_lim,
		       std::vector<float>& r_cnt, float& Ngal) const;
  
  void write_out(std::string file,CoreParam cp,Galaxies& gal);
  void alloc_cores(int size);
  void alloc_cp(int size);
  void write_out_lb_cores(std::string file);
  void write_out_lb_prop(std::string file_loc);
  void print_prop() const{
    std::cout<<"htag: "<<htag<<std::endl;
    std::cout<<"mass: "<<sod_mass<<std::endl;
    std::cout<<"radius: "<<sod_radius<<std::endl;
    std::cout<<"pos: "<<x<<" "<<y<<" "<<z<<std::endl;
    std::cout<<"redshift: "<<redshift<<std::endl;
    std::cout<<"step: "<<step<<std::endl;
    std::cout<<"z_i,m_i: "<<z_i<<" "<<m_i<<std::endl;
    std::cout<<"core_size: "<<core_size<<std::endl;
    std::cout<<"cp_size: "<<cp_size<<std::endl;
  }
  void print_cores() const{
    std::cout<<"core size: "<<core_size<<std::endl;
    for(int i =0;i<10;++i){
      std::cout<<i<<": "<<core_id[i]<<" "<<core_htag[i]<<" "<<core_m[i]<<" "
	       <<"\t"<<core_x[i]<<" "<<core_y[i]<<" "<<core_z[i]<<" "<<core_r[i]<<std::endl;

    }
    // for(int i =0;i<10;++i){
    //   std::cout<<i<<": "<<core_id[i]<<" "<<core_htag[i]<<" "<<core_m[i]<<" "<<core_cp_offset[i]<<" "<<core_cp_size[i]
    // 	       <<"\t"<<core_x[i]<<" "<<core_y[i]<<" "<<core_z[i]<<" "<<core_r[i]<<std::endl;
    // }
    // for(int i=core_size-10;i<core_size;++i){
    //   std::cout<<i<<": "<<core_id[i]<<" "<<core_htag[i]<<" "<<core_m[i]<<" "<<core_cp_offset[i]<<" "<<core_cp_size[i]
    // 	       <<"\t"<<core_x[i]<<" "<<core_y[i]<<" "<<core_z[i]<<" "<<core_r[i]<<std::endl;
    // }
  }
  void print_cp(){
    std::cout<<"cp_size: "<<cp_size<<std::endl;
    for(int i =0;i<10;++i){
      std::cout<<i<<": "<<cp_x[i]<<" "<<cp_y[i]<<" "<<cp_z[i]<<std::endl;
    }
    for(int i =cp_size-10;i<cp_size;++i){
      std::cout<<i<<": "<<cp_x[i]<<" "<<cp_y[i]<<" "<<cp_z[i]<<std::endl;
    }

  }
  void privatize_cp_usage_color(){
    cp_usage = new bool[cp_size];
    cp_color = new int[cp_size];
  }
  void privatize_cp_usage_color_done(){
    delete [] cp_usage;
    delete [] cp_color;
  }
  void move_cores_together();
};

struct FunctionParam{
  ZMR* zmr_cores;
  const ZMR* zmr_sdss;
  std::vector<Cluster>* clusters;
};

struct DescreteVariables{
  std::vector<float> mass_infall;
  std::vector<float> radius_disrupt;
  void sort(){
    std::sort(mass_infall.begin(),mass_infall.end());
    std::sort(radius_disrupt.begin(),radius_disrupt.end());
  }
  float get_m_infall(int val){
    val += mass_infall.size()/2.0;
    if(val>=mass_infall.size())
      return mass_infall.back();
    else if(val<0)
      return mass_infall.at(0);
    else
      return mass_infall.at(val);
  }
 float get_r_disrupt(int val){
    val += radius_disrupt.size()/2.0;
    if(val>=radius_disrupt.size())
      return radius_disrupt.back();
    else if(val<0)
      return radius_disrupt.at(0);
    else
      return radius_disrupt.at(val);
  }
};
ZMR zmr_cores;
ZMR zmr_sdss;

DescreteVariables dv;
  

std::map<int,Cores> all_cores;
std::map<int,dtk::ChainingMeshIndex> all_core_cms;
std::map<int,HaloCat> all_halocats;
std::map<int,Cores>  all_halo_cores;
std::map<int,CorePrtcls>  all_core_prtcls;
std::map<int,AccumPrtcls> all_accum_prtcls;
std::vector<Cluster> all_clusters;

struct LockedR_Disrupt{
  std::vector<float> m_infall;
  std::vector<float> r_disrupt;
};

//Params
dtk::Param param;
int rank;
int nproc;
std::vector<int> steps;
float step;
//std::string fof_loc;
std::string sod_loc;
bool sod_hdf5; //If false, sod_loc is a gio file. If true, it's a hdf5 file.
std::string core_loc;
std::string core_central_loc;
std::string zmr_loc;
std::string core_prtcl_loc;
std::string accum_prtcl_loc;
std::string cluster_loc;

bool use_central_infall;
bool force_central_as_galaxy;
bool force_center_on_central;
bool no_core_particles = true;

bool read_clusters_from_file;
bool write_clusters_to_file;
bool write_core_clusters_to_file;
float rL;
size_t chaining_mesh_grid_size;
float cluster_radial_volume;
bool scaled_cluster_radial_volume;
float scaled_cluster_radial_factor;

int min_disrupt_fof;
float m_infall;
float r_disrupt;
float r_fof; 
float r_merger;

bool fit_r_fof;
bool fit_r_merger;
int  fit_var_num;
bool lock_r_disrupt;

int diff_projection_num;

bool cost_gal_density;
bool cost_dgal_dr;
bool cost_gal_accum;
bool cost_abundance;
float expected_comov_abundance;
int  cost_min_cluster_num;
float cost_missed_centrals;

bool run_mcmc;
int omp_num_threads;
int chain_length;
std::vector<float> m_infall_constraints;
std::vector<float> r_disrupt_constraints;
std::vector<float> r_merger_constraints;
std::vector<float> r_fof_constraints;

float sigma_init_boost;
float sigma_m_infall;
float sigma_r_disrupt;
float sigma_r_fof;
float sigma_r_merger;

double simplex_step_size;
double simplex_cnvrg_size;
int max_iterations;
//the start/end/bin count for the bins
bool calc_likelihood;
std::vector<float> mi_bins_info;
std::vector<float> rd_bins_info;
std::vector<float> rm_bins_info;

std::vector<float> mi_bins;
std::vector<float> rd_bins;
std::vector<float> rm_bins;

int cm_num;
float ll;
//std::default_random_engine rnd_gen;
//std::uniform_real_distribution<float> rnd_theta(0,M_PI/4);
dtk::StepRedshift stepz(200,0,500);

void load_param(char* file_name);
void load_cores();
// void load_cores(std::string file_name,Cores& cores);
void load_halo_cats();
void load_zmr_sdss(std::string,ZMR& zmr_sdss);
void broadcast_zmr(ZMR& zmr);
void make_clusters();
void defrag_halos(int64_t* htags, int64_t size);
void move_together(float x0,float* x, int64_t size);
void make_zmr(const std::vector<Cluster>& clstrs,float m_infall,float r_disrupt,
	      float r_fof,float r_merger,ZMR& zmr,bool verbose=true);
void make_zmr(const std::vector<Cluster>& clstrs,CoreParam cp, ZMR& zmr,bool verbose=true);

double calc_diff(const ZMR& zmr1, const ZMR& zmr2, const std::map<int,Cores>& all_cores,const CoreParam& cp);
double calc_diff_gal_density(const ZMR& zmr1, const ZMR& zmr2);
double calc_diff_dgal_dr(const ZMR& zmr1, const ZMR& zmr2);
double calc_diff_gal_accum(const ZMR& zmr1, const ZMR& zmr2);
double calc_diff_abundance(const std::map<int,Cores>& all_cores,const CoreParam& cp);
float get_locked_r_disrupt(float m_infall,Cores cores);
double param_cost(double m_i,double r_d,double r_m,double r_fof);
double param_cost(CoreParam cp);
void calculate_likelihood_grid(std::vector<float> mi_bins,
			       std::vector<float> rd_bins,
			       std::vector<float> rm_bins,
			       std::string       out_loc,
			       CoreParam& max_lkhd_param);
CoreParam find_min(std::vector<Cluster>& clstrs, ZMR& zmr_cores, const ZMR& zmr_sdss,
		   float mi0,float rd0,float rfof0,float rm0);
template<typename T>
void write_array(std::ofstream& file,std::string name,T* data, int64_t size){
  file<<name;
  for(int64_t i =0;i<size;++i){
    file<<" "<<data[i]; 
  }
  file<<std::endl;
}
template<typename T>
void write_array(std::ofstream& file,std::string name,std::vector<T>& data){
  write_array(file,name,&data[0],data.size());
}
template<typename T>
void write_value(std::ofstream& file,std::string name,T data){
  file<<name<<" "<<data<<std::endl;
}


float fake_log_likelihood(CoreParam cp){
  float sigma_mi =1;
  float sigma_rd =1;
  float mi = (11.5-std::log10(cp.m_infall))/sigma_mi;
  float rd = (std::log10(0.06)-std::log10(cp.r_disrupt))/sigma_rd;;
  // std::cout<<"M_infall: "<<cp.m_infall<<"/"<<sigma_mi<<" -> "<<mi<<std::endl;
  // std::cout<<"R_disrupt: "<<cp.r_disrupt<<"/"<<sigma_rd<<" -> "<<rd<<std::endl;
  return mi*mi+rd*rd;
}

float rand_float();
void write_out_clusters(std::vector<Cluster>& clstrs, CoreParam cp);
void test_chaining_mesh(){
  ChainingMesh<float> cm;
  int num = cm_num;
  float rl = 100;
  int grid_i = 100;
  int grid_j = 100;
  int grid_k = 100;
  float idx = rl/grid_i;
  float jdy = rl/grid_j;
  float kdz = rl/grid_k;
  cm.set_limits(0,rl,0,rl,0,rl,grid_i,grid_j,grid_k,false);
  for(int i =0;i<10;++i){

    float x0=1000;//5+i*10;
    
    std::cout<<"\n"<<i<<": "<<x0<<std::endl;
    std::cout<<cm.index_i(x0)<<" "
	     <<cm.index_j(x0)<<" "
	     <<cm.index_k(x0)<<std::endl;
    std::vector<int> ijk = cm.index3d_from_pos(x0,x0,x0);
    std::cout<<ijk[0]<<" "
	     <<ijk[1]<<" "
	     <<ijk[2]<<std::endl;
    std::cout<<cm.index1d_from_pos(x0,x0,x0)<<std::endl;
  }
  float*x,*y,*z;
  bool* use;
  int* c,*w;
  x = new float[num];
  y = new float[num];
  z = new float[num];
  w = new int[num];
  use = new bool[num];
  c   = new int[num];
  std::srand(time(NULL)); 
  //std::srand(1000);
  for(int i =0;i<num;++i){
    x[i] = dtk::rand<float>()*5.0+.1;
    y[i] = dtk::rand<float>()*5.0+.1;
    z[i] = dtk::rand<float>()*.0+.1;
    use[i]=true;
    c[i]=i;
  }
  std::vector<std::vector<float> > fof_x;
  std::vector<std::vector<float> > fof_y;
  std::vector<std::vector<float> > fof_z;
  std::cout<<"a problem here?"<<std::endl;
  cm.set_data(x,y,z,num);
  std::cout<<"not here"<<std::endl;
  const size_t* srt = cm.get_srt();
  const int* cell_assignment = cm.get_cell_assignment();
  dtk::reorder(x,num,srt);
  dtk::reorder(y,num,srt);
  dtk::reorder(z,num,srt);
  std::cout<<"starting new fof finding..."<<std::endl;
  dtk::AutoTimer t1;
  //find_fof(ll,x,y,z,use,c,num);
  std::cout<<"done. time: "<<t1<<std::endl;
  dtk::AutoTimer t12;
  find_fof_cm(ll,x,y,z,use,c,num,cm);
  std::cout<<"done2. time: "<<t12<<std::endl;
  //for(int i =0;i<num;++i){
  //    std::cout<<i<<" "<<c[i]<<"= "<<x[i]<<" "<<y[i]<<" "<<z[i]<<" "<<std::endl;
  //  }
  int result_size=cm_num;
  
  std::cout<<"starting old fof_finding..."<<std::endl;
  dtk::AutoTimer t2;
  n2_merger3d(x,y,z,w,
	      &result_size,
	      ll,
	      NULL);
  std::cout<<"done. time: "<<t2<<std::endl;
  std::cout<<"starting newer fof_finding..."<<std::endl;
  dtk::AutoTimer t3;
  find_fof_srt(ll,x,y,z,c,cm_num);
  std::cout<<"done: "<<t3<<std::endl;
  //exit(0);  
  std::ofstream file("tmp/fof_test.param");
  write_array(file,"x",x,num);
  write_array(file,"y",y,num);
  write_array(file,"z",z,num);
  write_array(file,"c",c,num);
  file.close();
  //exit(0);
  
}
float log_vary(float data, float var){
  return std::pow(10,std::log10(data)+var);
}
CoreParam vary_params(CoreParam& cp,gsl_rng* r,double sigma,int step){
  CoreParam cp2 = cp;
  int variable_to_vary = step%fit_var_num; //vary only one variable
  switch(variable_to_vary){
  case -1: //vary all parameters
    cp2.m_infall = log_vary(cp2.m_infall,gsl_ran_gaussian(r,sigma));
    cp2.r_disrupt = log_vary(cp2.r_disrupt,gsl_ran_gaussian(r,sigma));
    cp2.r_fof += gsl_ran_gaussian(r,sigma);
    cp2.r_merger += gsl_ran_gaussian(r,sigma);
    break;
  case 0 :   
    cp2.m_infall = log_vary(cp2.m_infall,gsl_ran_gaussian(r,sigma));
    break;
  case 1:   
    cp2.r_disrupt = log_vary(cp2.r_disrupt,gsl_ran_gaussian(r,sigma));
    break;
  case 2:
    if(fit_r_fof)
      cp2.r_fof += gsl_ran_gaussian(r,sigma);
    else
      cp2.r_merger += gsl_ran_gaussian(r,sigma);
  case 3: 
    cp2.r_merger += gsl_ran_gaussian(r,sigma);
  }
  return cp2;
}
CoreParam vary_params(CoreParam& cp, gsl_rng* r, int step){
 CoreParam cp2 = cp;
 int variable_to_vary = step%fit_var_num; //vary only one variable
 switch(variable_to_vary){
 case -1: //vary all parameters for initinal mcmc scattering of the chain start positions
   cp2.m_infall = log_vary(cp2.m_infall,gsl_ran_gaussian(r,sigma_m_infall*sigma_init_boost));
   cp2.r_disrupt = log_vary(cp2.r_disrupt,gsl_ran_gaussian(r,sigma_r_disrupt*sigma_init_boost));
   cp2.r_fof = log_vary(cp2.r_fof,gsl_ran_gaussian(r,sigma_r_fof*sigma_init_boost));
   cp2.r_merger = log_vary(cp.r_merger,gsl_ran_gaussian(r,sigma_r_disrupt*sigma_init_boost));
   break;
 case 0 :   
   cp2.m_infall = log_vary(cp2.m_infall,gsl_ran_gaussian(r,sigma_m_infall));
   if(lock_r_disrupt){
     cp2.r_disrupt = get_locked_r_disrupt(cp2.m_infall,all_cores[steps[0]]);
   }
   break;
 case 1:   
   if(!lock_r_disrupt)
     cp2.r_disrupt = log_vary(cp2.r_disrupt,gsl_ran_gaussian(r,sigma_r_disrupt));
   else if(fit_r_fof)
     cp2.r_fof = log_vary(cp2.r_fof,gsl_ran_gaussian(r,sigma_r_fof));
   else
     cp2.r_merger = log_vary(cp2.r_merger,gsl_ran_gaussian(r,sigma_r_merger));
   break;
 case 2:
   if(!lock_r_disrupt){
     if(fit_r_fof){
       cp2.r_fof = log_vary(cp2.r_fof,gsl_ran_gaussian(r,sigma_r_fof));
     }
     else{
       cp2.r_merger = log_vary(cp2.r_merger,gsl_ran_gaussian(r,sigma_r_merger));
     }
   }
   else{
     cp2.r_merger = log_vary(cp2.r_merger,gsl_ran_gaussian(r,sigma_r_merger));
   }
   break;
 case 3: 
   cp2.r_merger = log_vary(cp2.r_merger,gsl_ran_gaussian(r,sigma_r_merger));
   break;
 }
 return cp2;
}
void run_MCMC(CoreParam cp,int steps,int seed,bool verbose,
	      float* mcmc_m_i,
	      float* mcmc_r_d,
	      float* mcmc_r_m,
	      float* mcmc_r_f,
	      float* mcmc_val,
	      int*   mcmc_step){
  dtk::AutoTimer t1;
  ZMR zmr_mcmc_core;
  zmr_mcmc_core.copy_bins(zmr_sdss);
  make_zmr(all_clusters,cp,zmr_mcmc_core,verbose);
  double current_state_log_likelihood;
  //current_state_log_likelihood = fake_log_likelihood(cp);//calc_diff(zmr_mcmc_core,zmr_sdss)+param_cost(cp);
  current_state_log_likelihood = calc_diff(zmr_mcmc_core,zmr_sdss,all_cores,cp)/2.0+param_cost(cp);
  float current_val;
  const gsl_rng_type * T;
  gsl_rng * r;
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set(r,seed);
  dtk::Timer t;
  dtk::AutoTimer t3;
  dtk::AutoTimer t2;
  for(int i=0;i<steps;++i){
    if(rank==0 and omp_get_thread_num() ==0){
      if(steps/10000 > 0  and i%(steps/10000)==0){
	t.stop();
	std::cout<<"\t\t"<<float(i)/float(steps)<<"  ("<<t2<<" / "<<t3<<")\r";
	t2.start();
      }
      else if(steps/100 > 0 and i%(steps/100)==0){
	t.stop();
	std::cout<<"\t\t"<<float(i)/float(steps)<<"  ("<<t2<<" / "<<t3<<")\r";
	t2.start();
      }
    }
    t.start();
    if(verbose)
      std::cout<<"\n\nStep: "<<i<<"  "<<current_state_log_likelihood<<std::endl;
    //double mcmc_sigma = 0.005;
    CoreParam cp2 = vary_params(cp,r,i);
    make_zmr(all_clusters,cp2,zmr_mcmc_core,verbose);
    double new_state_log_likelihood;
    //new_state_log_likelihood = fake_log_likelihood(cp2);
    new_state_log_likelihood = calc_diff(zmr_mcmc_core,zmr_sdss,all_cores,cp)/2.0+param_cost(cp2);
    double diff_state_log_likelihood = new_state_log_likelihood-current_state_log_likelihood;
    if(verbose)
      std::cout<<"diff loglh: "<<diff_state_log_likelihood<<std::endl;
    if(diff_state_log_likelihood <= 0.0){
      cp = cp2; //accept the new state
      if(verbose)
	std::cout<<"auto accept"<<std::endl;
      current_state_log_likelihood = new_state_log_likelihood;
    }
    else{
      double p = std::exp(-diff_state_log_likelihood);
      if(p>gsl_rng_uniform(r)){
	cp=cp2;
	if(verbose)
	  std::cout<<"accept: "<<p<<std::endl;
	current_state_log_likelihood = new_state_log_likelihood;
      }
      else{
	if(verbose)
	  std::cout<<"reject: "<<p<<std::endl;
      }
    }
    mcmc_m_i[i]=cp.m_infall;
    mcmc_r_d[i]=cp.r_disrupt;
    if(fit_r_merger)
      mcmc_r_m[i]=cp.r_merger;
    if(fit_r_fof)
      mcmc_r_f[i]=cp.r_fof;
    mcmc_val[i]=current_state_log_likelihood;
    mcmc_step[i]=i;
    t.stop();
    if(verbose)
      std::cout<<"****\n time: "<<t<<"  cost: "<<current_state_log_likelihood<<std::endl;
  }
  /*std::ofstream file("mcmc/test.param");
  write_array(file,"m_i",mcmc_m_i);
  write_array(file,"r_d",mcmc_r_d);
  write_array(file,"r_m",mcmc_r_d);
  write_array(file,"r_f",mcmc_r_d);
  file.close();*/
  gsl_rng_free(r);
  if(verbose)
    std::cout<<"done with mcmc: "<<t1<<std::endl;
}
struct MCMC_Args{
  CoreParam cp;
  int steps;
  int seed;
  bool verbose;
  float* mcmc_m_i;
  float* mcmc_r_d;
  float* mcmc_r_m;
  float* mcmc_r_f;
  int*   mcmc_step;
};
void* run_MCMC_pthread(void* ptr){
  MCMC_Args* args = (MCMC_Args*) ptr;
  //std::cout<<"new thread: "<<args->cp.m_infall<<" "<<args->cp.r_disrupt<<" "<<std::endl;
  //std::cout<<"steps: "<<args->steps<<std::endl;
  throw; //alloc mcmc_mi etc in run MCMC all pthreads
  run_MCMC(args->cp,args->steps,args->seed,args->verbose,
	   args->mcmc_m_i,
	   args->mcmc_r_d,
	   args->mcmc_r_m,
	   args->mcmc_r_f,
	   NULL, //Note: this will break the code!! but it will compile
	   args->mcmc_step);

  return NULL;
}
void run_MCMC_all_pthreads(CoreParam cp,int thread_num,int steps,bool verbose){
  if(verbose)
    std::cout<<"Starting all pthreads mcmc"<<std::endl;
  dtk::AutoTimer t;
  const gsl_rng_type * T;
  gsl_rng * r;
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set(r,time(NULL));
  int return_vals[thread_num];
  pthread_t threads[thread_num];
  for(int i =0;i<thread_num;++i){
    MCMC_Args args;
    args.cp = vary_params(cp,r,0.01,-1); //the -1 varies all parameters instead of by step iteration
    args.steps = steps;
    args.seed = 0;
    args.verbose = verbose;
    void* ptr = &args;
    std::cout<<"Starting new thread"<<std::endl;
    return_vals[i]=pthread_create(&threads[i],NULL,run_MCMC_pthread,ptr);
    //usleep(5000);
  }
  for(int i =0;i<thread_num;++i){
    pthread_join(threads[i],NULL);
  }
  gsl_rng_free(r);
  if(verbose)
    std::cout<<"\tdone all pthread mcmc. time: "<<t<<std::endl;
}
void run_MCMC_all_openmp(CoreParam cp,int steps, bool verbose){
 if(verbose and rank ==0)
   std::cout<<"Starting all openmc mcmc len: "<<steps<<std::endl;
  dtk::AutoTimer t;
  const gsl_rng_type * T;
  gsl_rng * r;
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set(r,time(NULL)*rank);
  int max_threads = omp_get_max_threads();
  CoreParam thread_cp[max_threads];
  if(verbose and rank==0)
    std::cout<<"Max num threads: "<<max_threads<<std::endl;
  for(int i =0;i<max_threads;++i){
    thread_cp[i] = vary_params(cp,r,-1);
  }
  float* mcmc_m_i = new float[steps*max_threads];
  float* mcmc_r_d = new float[steps*max_threads];
  float* mcmc_r_m;
  float* mcmc_r_f;
  int*   mcmc_thread_id = new int[steps*max_threads];
  int*   mcmc_step = new int[steps*max_threads];
  float* mcmc_val = new float[steps*max_threads];
  if(fit_r_merger)
    mcmc_r_m = new float[steps*max_threads];
  else 
    mcmc_r_m = NULL;
  if(fit_r_fof)
    mcmc_r_f = new float[steps*max_threads];
  else
    mcmc_r_f = NULL;
#pragma omp parallel firstprivate(all_clusters)
  {

    int tnum = omp_get_thread_num();
    int numthreads = omp_get_num_threads();
    int thread_id = numthreads*rank+tnum;
    //std::cout<<"this is thread "<<tnum<<"/"<<numthreads<<std::endl;
    for(int i =0;i<all_clusters.size();++i){
      all_clusters[i].privatize_cp_usage_color(); //these arrays are 
      //used in the zmr calculation so they must be private for each 
      //thread.
    }
    
    std::fill(mcmc_thread_id+(tnum*steps),
	      mcmc_thread_id+(tnum*steps+steps),
	      thread_id);
	      

    run_MCMC(thread_cp[tnum],steps,tnum,false,
	     mcmc_m_i+(tnum*steps),
	     mcmc_r_d+(tnum*steps),
	     mcmc_r_m+(tnum*steps),
	     mcmc_r_f+(tnum*steps),
	     mcmc_val+(tnum*steps),
	     mcmc_step + (tnum*steps));
  }
  std::string file_loc = dtk::make_str("output/",param.get_file_name(),"/mcmc.gio");
  std::cout<<"saving mcmc result to "<<file_loc<<std::endl;
  dtk::ensure_file_path(file_loc);
  std::cout<<"we wrote it out"<<std::endl;
  gio::GenericIO giofile(MPI_COMM_WORLD,file_loc);
  giofile.setNumElems(max_threads*steps);
  giofile.addVariable("mcmc_mass_infall", mcmc_m_i);
  giofile.addVariable("mcmc_r_disrupt", mcmc_r_d);
  if(fit_r_merger)
    giofile.addVariable("mcmc_r_merger", mcmc_r_m);
  if(fit_r_fof)
    giofile.addVariable("mcmc_r_fof", mcmc_r_f);
  giofile.addVariable("mcmc_walker_id",mcmc_thread_id);
  giofile.addVariable("mcmc_walker_step",mcmc_step);
  giofile.addVariable("mcmc_value",mcmc_val);
  giofile.write();
  giofile.close();
  gsl_rng_free(r);
  if(verbose and rank==0)
    std::cout<<"\tdone all openmp mcmc. time: "<<t<<std::endl;

}
void test_mcmc(CoreParam cp){
  /*CoreParam cp;
  cp.m_infall  = m_infall;
  cp.r_disrupt = r_disrupt;
  cp.r_fof     = r_fof;
  cp.r_merger  = r_merger;*/
  //run_MCMC(cp,20,0,true);
}
void broadcast_clusters(Cluster& cluster){
  //broadcast halo properties 
  MPI_Bcast(&cluster.htag,1,MPI_INT64_T,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.sod_mass,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.sod_radius,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.x,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.y,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.z,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.redshift,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.step,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.z_i,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.m_i,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.core_size,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.cp_size,1,MPI_INT,0,MPI_COMM_WORLD);
  if(rank != 0){ //alloc space for core properties & core particles
    cluster.alloc_cores(cluster.core_size);
    cluster.alloc_cp(cluster.cp_size);
  }
  //broadcast cores
  MPI_Bcast(&cluster.core_id,cluster.core_size,MPI_INT64_T,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.core_htag,cluster.core_size,MPI_INT64_T,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.core_x,cluster.core_size,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.core_y,cluster.core_size,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.core_z,cluster.core_size,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.core_r,cluster.core_size,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.core_m,cluster.core_size,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.core_step,cluster.core_size,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.core_cp_offset,cluster.core_size,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.core_cp_size,cluster.core_size,MPI_INT,0,MPI_COMM_WORLD);
  //broadcast core particles
  MPI_Bcast(&cluster.cp_x,cluster.core_size,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.cp_y,cluster.core_size,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&cluster.cp_z,cluster.core_size,MPI_FLOAT,0,MPI_COMM_WORLD);
}
void broadcast_all_clusters(){
  for(int i=0;i<all_clusters.size();++i){
    broadcast_clusters(all_clusters[i]);
  }
}
void broadcast_zmr_old(ZMR& zmr){
  MPI_Bcast(&zmr.z_size,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&zmr.m_size,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&zmr.r_size,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&zmr.zm_size,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&zmr.mr_size,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&zmr.zmr_size,1,MPI_INT,0,MPI_COMM_WORLD);

  dtk::broadcast_vector(zmr.zm_Ngal,0,MPI_COMM_WORLD);
  dtk::broadcast_vector(zmr.zm_Ngal_err,0,MPI_COMM_WORLD);
  dtk::broadcast_vector(zmr.zm_Ngal_var,0,MPI_COMM_WORLD);

  dtk::broadcast_vector(zmr.zmr_gal_counts,0,MPI_COMM_WORLD);

  dtk::broadcast_vector(zmr.zmr_gal_density,0,MPI_COMM_WORLD);
  dtk::broadcast_vector(zmr.zmr_gal_density_err,0,MPI_COMM_WORLD);
  dtk::broadcast_vector(zmr.zmr_gal_density_var,0,MPI_COMM_WORLD);

  dtk::broadcast_vector(zmr.zmr_gal_accum,0,MPI_COMM_WORLD);
  dtk::broadcast_vector(zmr.zmr_gal_accum_err,0,MPI_COMM_WORLD);
  dtk::broadcast_vector(zmr.zmr_gal_accum_var,0,MPI_COMM_WORLD);

  dtk::broadcast_vector(zmr.zmr_dgal_dr,0,MPI_COMM_WORLD);
  dtk::broadcast_vector(zmr.zmr_dgal_dr_err,0,MPI_COMM_WORLD);
  dtk::broadcast_vector(zmr.zmr_dgal_dr_var,0,MPI_COMM_WORLD);

  dtk::broadcast_vector(zmr.zm_counts,0,MPI_COMM_WORLD);
 
  MPI_Bcast(&zmr.centrals_found,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&zmr.centrals_expected,1,MPI_INT,0,MPI_COMM_WORLD);
}

void write_cluster(H5::H5File& file, Cluster& cluster, int num){
  std::string ss = "/cluster${num}/";
  ss = dtk::rep_str(ss,"${num}",num);
  dtk::write_hdf5(file,dtk::make_str(ss,"htag"),cluster.htag);
  dtk::write_hdf5(file,dtk::make_str(ss,"sod_mass"),cluster.sod_mass);
  dtk::write_hdf5(file,dtk::make_str(ss,"sod_radius"),cluster.sod_radius);
  dtk::write_hdf5(file,dtk::make_str(ss,"x"),cluster.x);
  dtk::write_hdf5(file,dtk::make_str(ss,"y"),cluster.y);
  dtk::write_hdf5(file,dtk::make_str(ss,"z"),cluster.z);
  dtk::write_hdf5(file,dtk::make_str(ss,"redshift"),cluster.redshift);
  dtk::write_hdf5(file,dtk::make_str(ss,"step"),cluster.step);
  dtk::write_hdf5(file,dtk::make_str(ss,"z_i"),cluster.z_i);
  dtk::write_hdf5(file,dtk::make_str(ss,"m_i"),cluster.m_i);
  dtk::write_hdf5(file,dtk::make_str(ss,"core_size"),cluster.core_size);
  dtk::write_hdf5(file,dtk::make_str(ss,"cp_size"),cluster.cp_size);
  dtk::write_hdf5(file,dtk::make_str(ss,"core_id"),        cluster.core_id,        cluster.core_size);
  dtk::write_hdf5(file,dtk::make_str(ss,"core_htag"),      cluster.core_htag,      cluster.core_size);
  dtk::write_hdf5(file,dtk::make_str(ss,"core_x"),         cluster.core_x,         cluster.core_size);
  dtk::write_hdf5(file,dtk::make_str(ss,"core_y"),         cluster.core_y,         cluster.core_size);
  dtk::write_hdf5(file,dtk::make_str(ss,"core_z"),         cluster.core_z,         cluster.core_size);
  dtk::write_hdf5(file,dtk::make_str(ss,"core_r"),         cluster.core_r,         cluster.core_size);
  dtk::write_hdf5(file,dtk::make_str(ss,"core_m"),         cluster.core_m,         cluster.core_size);
  dtk::write_hdf5(file,dtk::make_str(ss,"core_is_central"),cluster.core_is_central,cluster.core_size);
  dtk::write_hdf5(file,dtk::make_str(ss,"core_step"),      cluster.core_step,      cluster.core_size);
  if(!no_core_particles){
    dtk::write_hdf5(file,dtk::make_str(ss,"core_cp_offset"), cluster.core_cp_offset, cluster.core_size);
    dtk::write_hdf5(file,dtk::make_str(ss,"core_cp_size"),   cluster.core_cp_size,   cluster.core_size);
    dtk::write_hdf5(file,dtk::make_str(ss,"cp_x"),  cluster.cp_x,   cluster.cp_size);
    dtk::write_hdf5(file,dtk::make_str(ss,"cp_y"),  cluster.cp_y,   cluster.cp_size);
    dtk::write_hdf5(file,dtk::make_str(ss,"cp_z"),  cluster.cp_z,   cluster.cp_size);
  }

  //cluster.print_prop();
  //cluster.print_cores();
  //cluster.print_cp();

}
void read_cluster(H5::H5File& file, Cluster& cluster, int num){
 std::string ss = "cluster${num}/";
 ss = dtk::rep_str(ss,"${num}",num);
 // std::cout<<"reading in the non-array values"<<std::endl;
 dtk::read_hdf5(file,dtk::make_str(ss,"htag"),cluster.htag);
 dtk::read_hdf5(file,dtk::make_str(ss,"sod_mass"),cluster.sod_mass);
 dtk::read_hdf5(file,dtk::make_str(ss,"sod_radius"),cluster.sod_radius);
 dtk::read_hdf5(file,dtk::make_str(ss,"x"),cluster.x);
 dtk::read_hdf5(file,dtk::make_str(ss,"y"),cluster.y);
 dtk::read_hdf5(file,dtk::make_str(ss,"z"),cluster.z);
 dtk::read_hdf5(file,dtk::make_str(ss,"redshift"),cluster.redshift);
 dtk::read_hdf5(file,dtk::make_str(ss,"step"),cluster.step);
 dtk::read_hdf5(file,dtk::make_str(ss,"z_i"),cluster.z_i);
 dtk::read_hdf5(file,dtk::make_str(ss,"m_i"),cluster.m_i);
 dtk::read_hdf5(file,dtk::make_str(ss,"core_size"),cluster.core_size);
 dtk::read_hdf5(file,dtk::make_str(ss,"cp_size"),cluster.cp_size);
 //std::cout<<"Now starting on the core arrays"<<std::endl;
 //true means allocate the array
 dtk::read_hdf5(file,dtk::make_str(ss,"core_id"),        cluster.core_id,        true);
 dtk::read_hdf5(file,dtk::make_str(ss,"core_htag"),      cluster.core_htag,      true);
 dtk::read_hdf5(file,dtk::make_str(ss,"core_x"),         cluster.core_x,         true);
 dtk::read_hdf5(file,dtk::make_str(ss,"core_y"),         cluster.core_y,         true);
 dtk::read_hdf5(file,dtk::make_str(ss,"core_z"),         cluster.core_z,         true);
 dtk::read_hdf5(file,dtk::make_str(ss,"core_r"),         cluster.core_r,         true);
 dtk::read_hdf5(file,dtk::make_str(ss,"core_m"),         cluster.core_m,         true);
 dtk::read_hdf5(file,dtk::make_str(ss,"core_is_central"),cluster.core_is_central,true);
 dtk::read_hdf5(file,dtk::make_str(ss,"core_step"),     cluster.core_step,      true);
 if(!no_core_particles){
   dtk::read_hdf5(file,dtk::make_str(ss,"core_cp_offset"),cluster.core_cp_offset, true);
   dtk::read_hdf5(file,dtk::make_str(ss,"core_cp_size"),  cluster.core_cp_size,   true);
   
   // std::cout<<"Now starting on the cp arrays"<<std::endl;
   dtk::read_hdf5(file,dtk::make_str(ss,"cp_x"),  cluster.cp_x,   true);
   dtk::read_hdf5(file,dtk::make_str(ss,"cp_y"),  cluster.cp_y,   true);
   dtk::read_hdf5(file,dtk::make_str(ss,"cp_z"),  cluster.cp_z,   true);
 }
 cluster.cp_usage = new bool[cluster.cp_size];
 cluster.cp_color = new int[cluster.cp_size];
 
 //cluster.print_prop();
 //cluster.print_cores();
 //cluster.print_cp();
}
void write_clusters(std::string loc,std::vector<Cluster>& clusters){
  std::cout<<"starting to write clusters"<<std::endl;
  //dtk::ensure_file_path(loc);
  std::cout<<"ensure path done"<<std::endl;
  H5::H5File hfile(loc,H5F_ACC_TRUNC);
  dtk::write_hdf5(hfile,"/cluster_num",clusters.size());
  std::cout<<"Done writing cluster num"<<std::endl;
  for(int i= 0;i<clusters.size();++i){
    dtk::AutoTimer t;
    write_cluster(hfile,clusters[i],i);
    if(i%10000){
      std::cout<<"\t\t"<<i<<"/"<<clusters.size()<<": "<<dtk::fraction(i,clusters.size())<<std::endl;
      std::cout<<"\t\t"<<t<<std::endl;
    }
  }
  std::cout<<"done write clusters"<<std::endl;
}
void read_clusters(std::string loc, std::vector<Cluster>& clusters){
  if(rank==0)
    std::cout<<"reading in teh clusters"<<std::endl;
  H5::H5File hfile(loc,H5F_ACC_RDONLY);
  size_t cluster_num;
  dtk::read_hdf5(hfile,"/cluster_num",cluster_num);
  clusters.resize(cluster_num);
  for(int i=0;i<cluster_num;++i){
    read_cluster(hfile,clusters[i],i);
  }
  hfile.close();
}
void write_clusters_fast(std::string loc, std::vector<Cluster>& clusters){
  std::cout<<"Writing clusters to file quickly..."<<std::endl;
  dtk::AutoTimer t;
  size_t cluster_num = clusters.size();
  std::vector<int64_t> htag(cluster_num);
  std::vector<float>   sod_mass(cluster_num), sod_radius(cluster_num),
    x(cluster_num), y(cluster_num), z(cluster_num), redshift(cluster_num);
  std::vector<int> step(cluster_num), z_i(cluster_num), m_i(cluster_num);
  std::vector<size_t> core_offset(cluster_num), core_size(cluster_num);
  
  size_t current_core_offset = 0;
  std::cout<<"going over clusters"<<std::endl;
  for(size_t i=0;i<cluster_num;++i){
    htag.at(i) = clusters.at(i).htag;
    sod_mass.at(i) = clusters.at(i).sod_mass;
    sod_radius.at(i) = clusters.at(i).sod_radius;
    x.at(i)=clusters.at(i).x;
    y.at(i)=clusters.at(i).y;
    z.at(i)=clusters.at(i).z;
    redshift.at(i) = clusters.at(i).redshift;
    z_i.at(i) = clusters.at(i).z_i;
    m_i.at(i) = clusters.at(i).m_i;
    core_offset.at(i)=current_core_offset;
    core_size.at(i)=clusters.at(i).core_size;
    current_core_offset+=clusters.at(i).core_size;
  }

  size_t total_core_size = current_core_offset;
  std::cout<<"allocing space for cores: "<<total_core_size<<std::endl;
  std::vector<int64_t> core_id(total_core_size), core_htag(total_core_size);
  std::vector<float> core_x(total_core_size), core_y(total_core_size), core_z(total_core_size),
    core_r(total_core_size), core_m(total_core_size);
  std::vector<int> core_is_central(total_core_size), core_step(total_core_size);
  for(size_t i=0;i<cluster_num;++i){
    std::cout<<"cores for "<<i<<"/"<<cluster_num<<std::endl;    
    size_t offset = core_offset.at(i);
    size_t core_num = core_size.at(i);
    std::cout<<"\toffset: "<<offset<<"/"<<total_core_size<<" number of cores: "<<core_num<<std::endl;
    dtk::copy_n(clusters.at(i).core_id, core_num, &core_id[offset]); 
    dtk::copy_n(clusters.at(i).core_htag, core_num, &core_htag[offset]);
    dtk::copy_n(clusters.at(i).core_x, core_num, &core_x[offset]); 
    dtk::copy_n(clusters.at(i).core_y, core_num, &core_y[offset]); 
    dtk::copy_n(clusters.at(i).core_z, core_num, &core_z[offset]); 
    dtk::copy_n(clusters.at(i).core_r, core_num, &core_r[offset]); 
    dtk::copy_n(clusters.at(i).core_m, core_num, &core_m[offset]); 
    dtk::copy_n(clusters.at(i).core_is_central, core_num, &core_is_central[offset]); 
    dtk::copy_n(clusters.at(i).core_step, core_num, &core_step[offset]); 
  } 
  H5::H5File hfile(loc,H5F_ACC_TRUNC);
  std::cout<<"writing out cluster number"<<std::endl;
  dtk::write_hdf5(hfile,"/cluster_num", cluster_num);
  // clusters
  std::cout<<"writing out clusters"<<std::endl;
  dtk::write_hdf5(hfile,"/cluster/htag",      htag);
  dtk::write_hdf5(hfile,"/cluster/sod_mass",  sod_mass);
  dtk::write_hdf5(hfile,"/cluster/sod_radius",sod_radius);
  dtk::write_hdf5(hfile,"/cluster/x",         x);
  dtk::write_hdf5(hfile,"/cluster/y",         y);
  dtk::write_hdf5(hfile,"/cluster/z",         z);
  dtk::write_hdf5(hfile,"/cluster/z_i",       z_i);
  dtk::write_hdf5(hfile,"/cluster/m_i",       m_i);

  dtk::write_hdf5(hfile,"/cluster/redshift",  redshift);
  dtk::write_hdf5(hfile,"/cluster/step",      step);
  dtk::write_hdf5(hfile,"/cluster/core_offset",core_offset);
  dtk::write_hdf5(hfile,"/cluster/core_size",  core_size);
  // Cores
  std::cout<<"writing out cores"<<std::endl;
  dtk::write_hdf5(hfile,"/cores/core_id",  core_id);
  dtk::write_hdf5(hfile,"/cores/core_htag",  core_htag);
  dtk::write_hdf5(hfile,"/cores/core_x",  core_x);
  dtk::write_hdf5(hfile,"/cores/core_y",  core_y);
  dtk::write_hdf5(hfile,"/cores/core_z",  core_z);
  dtk::write_hdf5(hfile,"/cores/core_r",  core_r);
  dtk::write_hdf5(hfile,"/cores/core_m",  core_m);
  dtk::write_hdf5(hfile,"/cores/core_is_central",  core_is_central);
  dtk::write_hdf5(hfile,"/cores/core_step",  core_step);
  std::cout<<"\tdone: "<<t<<std::endl;
}
void read_clusters_fast(std::string loc, std::vector<Cluster>& clusters){
  std::cout<<"Reading clusters quickly via new method: ..."<<std::endl;
  dtk::AutoTimer t;
  std::vector<int64_t> htag, core_id, core_htag;
  std::vector<size_t> core_offset, core_size;
  std::vector<int>  step, core_is_central, core_step, z_i, m_i;
  std::vector<float> sod_mass, sod_radius, x, y, z, redshift,
    core_x, core_y, core_z, core_r, core_m;
  size_t cluster_num;
  H5::H5File hfile(loc,H5F_ACC_RDONLY);
  dtk::read_hdf5(hfile, "/cluster_num", cluster_num);
  // cluster info
  dtk::read_hdf5(hfile, "/cluster/htag", htag);
  dtk::read_hdf5(hfile, "/cluster/sod_mass", sod_mass);
  dtk::read_hdf5(hfile, "/cluster/sod_radius", sod_radius);
  dtk::read_hdf5(hfile, "/cluster/x", x);
  dtk::read_hdf5(hfile, "/cluster/y", y);
  dtk::read_hdf5(hfile, "/cluster/z", z);
  dtk::read_hdf5(hfile, "/cluster/z_i", z_i);
  dtk::read_hdf5(hfile, "/cluster/m_i", m_i);
  dtk::read_hdf5(hfile, "/cluster/redshift", redshift);
  dtk::read_hdf5(hfile, "/cluster/core_offset", core_offset);
  dtk::read_hdf5(hfile, "/cluster/core_size", core_size);
  dtk::read_hdf5(hfile, "/cluster/step", step);
  // core info
  dtk::read_hdf5(hfile, "/cores/core_id", core_id);
  dtk::read_hdf5(hfile, "/cores/core_htag", core_htag);
  dtk::read_hdf5(hfile, "/cores/core_x", core_x);
  dtk::read_hdf5(hfile, "/cores/core_y", core_y);
  dtk::read_hdf5(hfile, "/cores/core_z", core_z);
  dtk::read_hdf5(hfile, "/cores/core_r", core_r);
  dtk::read_hdf5(hfile, "/cores/core_m", core_m);
  dtk::read_hdf5(hfile, "/cores/core_is_central", core_is_central);
  dtk::read_hdf5(hfile, "/cores/core_step", core_step);
  clusters.resize(cluster_num);
  std::cout<<"Cluster number: "<<cluster_num<<std::endl;
  for(int i =0;i<cluster_num;++i){
    Cluster& clstr = clusters[i];
    clstr.htag = htag[i];
    clstr.sod_mass = sod_mass[i];
    clstr.sod_radius = sod_radius[i];
    clstr.x = x[i];
    clstr.y = y[i];
    clstr.z = z[i];
    clstr.redshift = redshift[i];
    clstr.core_size = core_size[i];
    clstr.step = step[i];
    size_t clstr_core_offset = core_offset[i];
    // std::cout<<"cluster["<<i<<"] pos: "<<clstr.x<<" "<<clstr.y<<" "<<clstr.z<<std::endl;
    // std::cout<<"\t num cores: "<<clstr.core_size<<std::endl;
    // core data
    clstr.core_id = new int64_t[clstr.core_size];
    clstr.core_htag = new int64_t[clstr.core_size];
    clstr.core_x = new float[clstr.core_size];
    clstr.core_y = new float[clstr.core_size];
    clstr.core_z = new float[clstr.core_size];
    clstr.core_r = new float[clstr.core_size];
    clstr.core_m = new float[clstr.core_size];
    clstr.core_is_central = new int[clstr.core_size];
    clstr.core_step = new int[clstr.core_size];
    // load in core data
    dtk::copy_n(&core_id[clstr_core_offset], clstr.core_size, clstr.core_id);
    dtk::copy_n(&core_htag[clstr_core_offset], clstr.core_size, clstr.core_htag);
    dtk::copy_n(&core_x[clstr_core_offset], clstr.core_size, clstr.core_x);
    dtk::copy_n(&core_y[clstr_core_offset], clstr.core_size, clstr.core_y);
    dtk::copy_n(&core_z[clstr_core_offset], clstr.core_size, clstr.core_z);
    dtk::copy_n(&core_m[clstr_core_offset], clstr.core_size, clstr.core_m);
    dtk::copy_n(&core_r[clstr_core_offset], clstr.core_size, clstr.core_r);
    dtk::copy_n(&core_is_central[clstr_core_offset], clstr.core_size, clstr.core_is_central);
    dtk::copy_n(&core_step[clstr_core_offset], clstr.core_size, clstr.core_step);
    clstr.move_cores_together();
    // for(int j = 0;j<clstr.core_size;++j){
    //   // std::cout<<"\t\t["<<j<<"]: "<<clstr.core_x[j]<<" "<<clstr.core_y[j]<<" "<<clstr.core_z[j]<<std::endl;
    //   // std::cout<<"\t\t\t"<<clstr.core_m[j]<<" "<<clstr.core_r[j]<<std::endl;
    // }
  }
  std::cout<<"\tdone: "<<t<<std::endl;
}
void bcast_cluster(Cluster& cluster, MPI_Datatype mpi_cluster_type,int root, MPI_Comm comm){
  int myrank;
  MPI_Comm_rank(comm,&myrank);
  MPI_Bcast(&cluster,1,mpi_cluster_type,root,comm);

  //bcasting cores
  if(myrank != root){
    cluster.core_id =   new int64_t[cluster.core_size];
    cluster.core_htag = new int64_t[cluster.core_size];
    cluster.core_x = new float[cluster.core_size];
    cluster.core_y = new float[cluster.core_size];
    cluster.core_z = new float[cluster.core_size];
    cluster.core_r = new float[cluster.core_size];
    cluster.core_m = new float[cluster.core_size];
    cluster.core_step =      new int[cluster.core_size];
    cluster.core_cp_offset = new int[cluster.core_size];
    cluster.core_cp_size =   new int[cluster.core_size];
  }
  MPI_Bcast(cluster.core_id,        cluster.core_size, MPI_INT64_T, root,comm);
  MPI_Bcast(cluster.core_htag,      cluster.core_size, MPI_INT64_T, root,comm);
  MPI_Bcast(cluster.core_x,         cluster.core_size, MPI_FLOAT,   root,comm);
  MPI_Bcast(cluster.core_y,         cluster.core_size, MPI_FLOAT,   root,comm);
  MPI_Bcast(cluster.core_z,         cluster.core_size, MPI_FLOAT,   root,comm);
  MPI_Bcast(cluster.core_r,         cluster.core_size, MPI_FLOAT,   root,comm);
  MPI_Bcast(cluster.core_m,         cluster.core_size, MPI_FLOAT,   root,comm);
  MPI_Bcast(cluster.core_step,      cluster.core_size, MPI_INT,     root,comm);
  MPI_Bcast(cluster.core_cp_offset, cluster.core_size, MPI_INT,     root,comm);
  MPI_Bcast(cluster.core_cp_size,   cluster.core_size, MPI_INT,     root,comm);
  //bcasting core particles
  if(myrank != root){
    cluster.cp_x     = new float[cluster.cp_size];
    cluster.cp_y     = new float[cluster.cp_size];
    cluster.cp_z     = new float[cluster.cp_size];
    cluster.cp_usage = new  bool[cluster.cp_size];
    cluster.cp_color = new   int[cluster.cp_size];
  }
  MPI_Bcast(cluster.cp_x,         cluster.cp_size, MPI_FLOAT,   root,comm);
  MPI_Bcast(cluster.cp_y,         cluster.cp_size, MPI_FLOAT,   root,comm);
  MPI_Bcast(cluster.cp_z,         cluster.cp_size, MPI_FLOAT,   root,comm);
}
void bcast_clusters(std::vector<Cluster>& clusters, int root, MPI_Comm comm){
  if(nproc == 1){
    std::cout<<"No need to bcast clusters"<<std::endl;
    return;
  }

  if(rank==0)
    std::cout<<"Starting to Bcast clusters"<<std::endl;
  dtk::AutoTimer t;
  const int num_elems = 12;
  MPI_Datatype types[num_elems] = {MPI_INT64_T,//htag
			    MPI_FLOAT,//mass
			    MPI_FLOAT,//rad
			    MPI_FLOAT,//x
			    MPI_FLOAT,//y
			    MPI_FLOAT,//z
			    MPI_FLOAT,//redshift
			    MPI_INT,//step
			    MPI_INT,//z_i
			    MPI_INT,//m_i
			    MPI_INT64_T,//core_size
			    MPI_INT64_T};//cp_size
  int lengths[num_elems] = {1,1,1,1, 1,1,1,1, 1,1,1,1};
  MPI_Aint offsets[num_elems];
  offsets[0] = offsetof(Cluster,htag);
  offsets[1] = offsetof(Cluster,sod_mass);
  offsets[2] = offsetof(Cluster,sod_radius);
  offsets[3] = offsetof(Cluster,x);
  offsets[4] = offsetof(Cluster,y);
  offsets[5] = offsetof(Cluster,z);
  offsets[6] = offsetof(Cluster,redshift);
  offsets[7] = offsetof(Cluster,step);
  offsets[8] = offsetof(Cluster,z_i);
  offsets[9] = offsetof(Cluster,m_i);
  offsets[10] = offsetof(Cluster,core_size);
  offsets[11] = offsetof(Cluster,cp_size);

  MPI_Datatype mpi_cluster_type;
  MPI_Type_create_struct(num_elems,lengths,offsets,types,&mpi_cluster_type);
  MPI_Type_commit(&mpi_cluster_type);

  int cluster_num=clusters.size();
  MPI_Bcast(&cluster_num,1,MPI_INT,root,comm);
  clusters.resize(cluster_num);
  for(int i=0;i<cluster_num;++i){
    bcast_cluster(clusters[i],mpi_cluster_type,root,comm);
  }
  if(rank==0)
    std::cout<<"\tdone bcasting clusters. time: "<<t<<std::endl;
}
void write_core_params(CoreParam cp, std::string file_loc);
float get_locked_r_disrupt(float m_infall,Cores cores);
int core_fit(char* param_fname);

//Main...you know the thing that runs
//everything...
int main(int argc, char** argv){
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  dtk::AutoTimer t1;
  for(int i =1;i<argc;++i){
    std::cout<<"running: "<<argv[i]<<std::endl;
    core_fit(argv[i]);
  }
  std::cout<<"\n\nAll param files done: "<<t1<<std::endl;
  MPI_Finalize();
}
int core_fit(char* param_fname){


  dtk::Timer t;t.start();
  if(rank==0){
    std::cout<<"Fitting cores"<<std::endl;
    std::cout<<"num ranks: "<<nproc<<std::endl;
  }
  load_param(param_fname);
  std::cout<<"tring this out.."<<std::endl;
  std::string file_loc = dtk::make_str("output/",param.get_file_name(),"/mcmc.gio");
  dtk::ensure_file_path(file_loc);
  std::cout<<"done"<<std::endl;
  load_zmr_sdss(zmr_loc,zmr_sdss);
  if(read_clusters_from_file){
    if(rank == 0){
      read_clusters_fast(cluster_loc,all_clusters);
      std::cout<<"all cluster size: "<<all_clusters.size()<<std::endl;
    }
    if(cost_abundance | lock_r_disrupt)
      load_cores();
  }
  else{
    //we need to get the halo props & core info 
    if(rank == 0){
      load_halo_cats();
      load_cores();
      make_clusters();
      if(write_clusters_to_file && !read_clusters_from_file)
	//write_clusters(cluster_loc,all_clusters);
	write_clusters_fast(cluster_loc,all_clusters);
    }
  }
  bcast_clusters(all_clusters,0,MPI_COMM_WORLD);
  CoreParam fit_cp;
  fit_cp.m_infall = m_infall;
  fit_cp.r_disrupt = r_disrupt;
  fit_cp.r_merger = r_merger;
  fit_cp.r_fof    = r_fof;
  std::cout<<"----:"<<std::endl;
  make_zmr(all_clusters,m_infall,r_disrupt,r_fof,r_merger,zmr_cores,false);
  std::cout<<"not make zmr"<<std::endl;
  double cost1 = calc_diff(zmr_sdss, zmr_cores,all_cores,fit_cp);
  double cost2 = calc_diff(zmr_cores,zmr_sdss,all_cores,fit_cp);
  std::cout<<"it's not the cost"<<std::endl;
  if(max_iterations > 0){
    fit_cp = find_min(all_clusters, zmr_cores, zmr_sdss, std::log10(m_infall),
		      std::log10(r_disrupt), std::log10(r_fof), std::log10(r_merger));
  }
  make_zmr(all_clusters,fit_cp,zmr_cores,false);
  std::cout<<"this is the cost: "<<calc_diff(zmr_sdss, zmr_cores,all_cores,fit_cp)<<std::endl;

  std::cout<<"fit parameters: "<<std::endl;
  std::cout<<"\tM_infall: "<<fit_cp.m_infall<<std::endl;
  std::cout<<"\tR_disrupt: "<<fit_cp.r_disrupt<<std::endl;
  if(fit_r_fof)
    std::cout<<"\tR_fof: "<<fit_cp.r_fof<<std::endl;
  if(fit_r_merger)
    std::cout<<"\tR_merger: "<<fit_cp.r_merger<<std::endl;
  std::cout<<"writing out zmrs"<<std::endl;
  std::string sdss_file_loc = dtk::make_str("output/",param.get_file_name(),"/zmr_sdss.param");
  std::string cores_file_loc = dtk::make_str("output/",param.get_file_name(),"/zmr_cores.param");
  // dtk::ensure_file_path(sdss_file_loc);
  // dtk::ensure_file_path(cores_file_loc);
  zmr_sdss.write_txt(sdss_file_loc);
  zmr_cores.write_txt(cores_file_loc);
  write_core_params(fit_cp,dtk::make_str("output/",param.get_file_name(),"/fit_core_params.hdf5"));
  std::cout<<"wrote out zmrs..."<<std::endl;
  if(write_core_clusters_to_file)
    write_out_clusters(all_clusters,fit_cp);
  if(run_mcmc)
    run_MCMC_all_openmp(fit_cp,chain_length,true);
  if(calc_likelihood){
    CoreParam max_lkhd_param;
    calculate_likelihood_grid(mi_bins,rd_bins,rm_bins,dtk::make_str("output/",param.get_file_name(),"/lgrid.param"), max_lkhd_param);
    make_zmr(all_clusters, max_lkhd_param, zmr_cores, false);
    std::string lkhd_cores_file_loc = dtk::make_str("output/", param.get_file_name(), "/zmr_lkhd_cores.param");
    zmr_cores.write_txt(lkhd_cores_file_loc);
    //write_core_param(max_lkhd_param
  }
  

  //test_mcmc(fit_cp);  
  
  while(false){

    CoreParam cparam = find_min(all_clusters, zmr_cores, zmr_sdss, std::log10(m_infall),
				std::log10(r_disrupt), std::log10(r_fof), std::log10(r_merger));
    std::cout<<"result param: "<<cparam.m_infall<<" "<<cparam.r_disrupt<<" "<<cparam.r_fof
	     <<" "<<cparam.r_merger<<std::endl;
    make_zmr(all_clusters,cparam,zmr_cores);
     //write_out_clusters(all_clusters,cparam);
    if(rank==0)
      std::cout<<"feel free change the fitting params in the file"<<std::endl;
    //if(!dtk::ask_continue())
      break;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  t.stop();
  if(rank==0)
    std::cout<<" total time "<<t<<std::endl;
}




void load_param(char* file_name){
  dtk::Timer t;t.start();
  if(rank==0){
    std::cout<<"Loading params"<<std::endl;
    param.load(file_name);
    std::cout<<"done loading param"<<std::endl;
  }
  dtk::broadcast_param(param,0,MPI_COMM_WORLD);

  //fof_loc = param.get<std::string>("fof_loc");
  sod_loc = param.get<std::string>("sod_loc");
  if(param.has("sod_hdf5"))
    sod_hdf5 = param.get<bool>("sod_hdf5");
  else
    sod_hdf5 = false;
  core_loc = param.get<std::string>("core_loc");
  core_central_loc = param.get<std::string>("core_central_loc");
  zmr_loc = param.get<std::string>("zmrh5_loc");
  steps = param.get_vector<int>("steps");
  step  = param.get<int>("step");
  core_prtcl_loc = param.get<std::string>("core_prtcl_loc");
  accum_prtcl_loc = param.get<std::string>("accum_prtcl_loc");
  cluster_loc = param.get<std::string>("cluster_loc");
  
  use_central_infall = param.get<bool>("use_central_infall");
  force_central_as_galaxy = param.get<bool>("force_central_as_galaxy");
  force_center_on_central = param.get<bool>("force_center_on_central");
  if(use_central_infall || force_central_as_galaxy || force_center_on_central){
    std::cout<<"The use of use_central_infall, force_central_as_galaxy, force_center_on_central is not allowed"<<std::endl;
    std::cout<<"Maybe force_center_on_central may work, but needs to be looked into"<<std::endl;
    throw;
  }
  read_clusters_from_file = param.get<bool>("read_clusters_from_file");
  write_clusters_to_file = param.get<bool>("write_clusters_to_file");
  write_core_clusters_to_file = param.get<bool>("write_core_clusters_to_file");

  rL = param.get<float>("rL");
  chaining_mesh_grid_size = param.get<int>("chaining_mesh_grid_size");
  cluster_radial_volume = param.get<float>("cluster_radial_volume");
  scaled_cluster_radial_volume = param.get<bool>("scaled_cluster_radial_volume");
  scaled_cluster_radial_factor = param.get<float>("scaled_cluster_radial_factor");
  
  min_disrupt_fof = param.get<float>("min_disrupt_fof");

  m_infall = pow(10,param.get<float>("m_infall"));
  
  r_disrupt = param.get<float>("r_disrupt");
  r_fof     = param.get<float>("r_fof");
  r_merger  = param.get<float>("r_merger");
  
  fit_r_fof    = param.get<bool>("fit_r_fof");
  fit_r_merger = param.get<bool>("fit_r_merger");
  lock_r_disrupt = param.get<bool>("lock_r_disrupt");
  diff_projection_num  = param.get<int>("diff_projection_num");
  
  cost_gal_density     = param.get<bool>("cost_gal_density");
  cost_dgal_dr         = param.get<bool>("cost_dgal_dr");
  cost_gal_accum       = param.get<bool>("cost_accum_gal");
  cost_abundance       = param.get<bool>("cost_abundance");
  cost_min_cluster_num = param.get<int>("cost_min_cluster_num");
  cost_missed_centrals = param.get<float>("cost_missed_centrals");
  expected_comov_abundance = param.get<float>("expected_comov_abundance");
  run_mcmc             = param.get<bool>("run_mcmc");
  chain_length         = param.get<int>("chain_length");
  omp_num_threads       = param.get<int>("omp_num_threads");
  m_infall_constraints = param.get_vector<float>("m_infall_constraints");
  r_disrupt_constraints= param.get_vector<float>("r_disrupt_constraints");
  r_merger_constraints = param.get_vector<float>("r_merger_constraints");
  r_fof_constraints    = param.get_vector<float>("r_fof_constraints");
  for(int i =0;i<2;++i){
    m_infall_constraints[i] = std::pow(10,m_infall_constraints[i]);
    r_merger_constraints[i] = std::log10(r_merger_constraints[i]);
  }
  sigma_init_boost = param.get<float>("sigma_init_boost");
  sigma_m_infall = param.get<float>("sigma_m_infall");
  sigma_r_disrupt = param.get<float>("sigma_r_disrupt");
  sigma_r_fof = param.get<float>("sigma_r_fof");
  sigma_r_merger = param.get<float>("sigma_r_merger");
  
  calc_likelihood = param.get<bool>("calc_likelihood");
  mi_bins_info = param.get_vector<float>("mi_bins_info");
  rd_bins_info = param.get_vector<float>("rd_bins_info");
  rm_bins_info = param.get_vector<float>("rm_bins_info");
  mi_bins = dtk::logspace(mi_bins_info[0],mi_bins_info[1],mi_bins_info[2]);
  rd_bins = dtk::linspace(rd_bins_info[0],rd_bins_info[1],rd_bins_info[2]);
  rm_bins = dtk::linspace(rm_bins_info[0],rm_bins_info[1],rm_bins_info[2]);

  simplex_step_size = param.get<double>("simplex_step_size");
  simplex_cnvrg_size= param.get<double>("simplex_cnvrg_size");
  max_iterations = param.get<int>("max_iterations");
  std::cout<<"we got here"<<std::endl;
  fit_var_num = 2;
  if(lock_r_disrupt)
    --fit_var_num;
  if(fit_r_fof)
    ++fit_var_num;
  if(fit_r_merger)
    ++fit_var_num;
  if(omp_num_threads != -1)
    omp_set_num_threads(omp_num_threads);
  t.stop();
  if(rank==0)
    std::cout<<"\tdone. time: "<<t<<std::endl;
}
void load_cores(){
  dtk::Timer t;t.start();
  std::cout<<"loading core catalog"<<std::endl;
  std::cout<<steps.size()<<std::endl;
  for(int i=0;i<steps.size();++i){
    std::string file = dtk::rep_str(core_loc,"${step}",steps[i]);
    Cores cores;
    dtk::read_gio_quick(file,"fof_halo_tag",cores.host_htag, cores.size);
    dtk::read_gio_quick(file,"core_tag",cores.ctag,     cores.size);
    dtk::read_gio_quick(file,"radius",cores.radius,     cores.size);
    dtk::read_gio_quick(file,"x",cores.x,               cores.size);
    dtk::read_gio_quick(file,"y",cores.y,               cores.size);
    dtk::read_gio_quick(file,"z",cores.z,               cores.size);
    dtk::read_gio_quick(file,"infall_mass",cores.infall_mass,  cores.size);
    dtk::read_gio_quick(file,"infall_step",cores.infall_step,  cores.size);
    dtk::read_gio_quick(file,"infall_fof_halo_tag",cores.infall_htag,cores.size);
    cores.is_central = new int[cores.size];
    cores.compact = new int[cores.size];
    if(use_central_infall || force_central_as_galaxy || force_center_on_central){
      int64_t central_size;
      int64_t* core_central_id;
      int64_t* srt;
      file = dtk::rep_str(core_central_loc,"${step}",steps[i]);
      dtk::read_gio_quick(file,"core_tag",core_central_id,central_size);
      dtk::read_gio_quick(file,"core_central_htag",cores.central_htag,central_size);
      dtk::read_gio_quick(file,"core_central_step",cores.central_step,central_size);
      dtk::read_gio_quick(file,"core_central_mass",cores.central_mass,central_size);
      dtk::read_gio_quick(file,"core_is_central",cores.is_central,central_size);
      srt = dtk::arg_sort(core_central_id,central_size);
      for(int j =0;j<cores.size;++j){
	int64_t indx = dtk::search_sorted(core_central_id,cores.ctag[j],srt,static_cast<size_t>(central_size));
	if(indx != -1 && core_central_id[indx] == cores.ctag[j]){
	  if(core_central_id[indx] == 287418457){
	    std::cout<<"\nreplacing "<<cores.infall_mass[j]<<" with "<<cores.central_mass[indx]<<std::endl;
	    std::cout<<cores.infall_htag[j]<<" -> "<<cores.central_htag[indx]<<std::endl;
	    std::cout<<cores.infall_step[j]<<" -> "<<cores.central_step[indx]<<std::endl;
	    std::cout<<cores.ctag[j]<<" =?= "<<core_central_id[indx]<<std::endl;
	    std::cout<<cores.is_central[indx]<<std::endl;
	  }
	  if(use_central_infall)
	    cores.infall_mass[j]=cores.central_mass[indx];
	  //cores.infall_htag[j]=core_central_htag[indx];
	  //cores.infall_step[j]=core_central_step[indx];
	}
	else{
	  std::cout<<"Hmmm? core not found in central core catalog?!"<<std::endl;
	  dtk::pause();
	}
      }

      delete [] core_central_id;
      delete [] srt;
    }
    // New cores now have the correct infall mass for centrals
    // else{
    //   cores.central_mass = new float[cores.size];
    //   std::cout<<"Do not use core central infall. "<<std::endl;
    // }
    //TODO: Deal with the fragments correctly/test the effect. 
    defrag_halos(cores.host_htag,cores.size);
    defrag_halos(cores.infall_htag,cores.size);
    std::cout<<"defrag is fine"<<std::endl;    
    //Note: Do we need to reorder by htag?!
    // int64_t* srt = dtk::new_sort_array<int64_t>(cores.size);
    // std::sort(srt,srt+cores.size,dtk::Sorter<int64_t>(cores.host_htag));
    // dtk::reorder(cores.host_htag,cores.size,srt);
    // dtk::reorder(cores.ctag,cores.size,srt);
    // dtk::reorder(cores.x,cores.size,srt);
    // dtk::reorder(cores.y,cores.size,srt);
    // dtk::reorder(cores.z,cores.size,srt);
    // dtk::reorder(cores.radius,cores.size,srt);
    // dtk::reorder(cores.infall_mass,cores.size,srt);
    // dtk::reorder(cores.infall_htag,cores.size,srt);
    // dtk::reorder(cores.infall_step,cores.size,srt);
    // if(use_central_infall || force_central_as_galaxy || force_center_on_central){
    //   dtk::reorder(cores.central_mass,cores.size,srt);
    //   dtk::reorder(cores.central_htag,cores.size,srt);
    //   dtk::reorder(cores.central_step,cores.size,srt);
    // }
    // dtk::reorder(cores.is_central,cores.size,srt);

    float  len_xyz[3] = {rL, rL, rL};
    size_t len_ijk[3] = {chaining_mesh_grid_size, chaining_mesh_grid_size, chaining_mesh_grid_size};
    dtk::ChainingMeshIndex cm(len_xyz, len_ijk, 3);
    float* data_xyz[3] = {cores.x, cores.y, cores.z};
    cm.place_onto_mesh(data_xyz, cores.size);
    // std::cout<<"min/max radius: "<<dtk::min(cores.radius,cores.size)
    // 	     <<"/"<<dtk::max(cores.radius,cores.size)<<std::endl;
    // std::cout<<"min/max infall mass: "<<dtk::min(cores.infall_mass,cores.size)
    // 	     <<"/"<<dtk::max(cores.infall_mass,cores.size)<<std::endl;


    all_cores[steps[i]]=cores;
    all_core_cms[steps[i]]=cm;
    std::cout<<"here is okay"<<std::endl;
  }
  t.stop();
  std::cout<<"\tdone. time: "<<t<<std::endl;

}
void load_core_prtcls(int step){
  if(all_core_prtcls.count(step)==0){
    dtk::Timer t;t.start();
    std::cout<<"Loading cores step "<<step<<"...";
    //not loaded yet -> lets load that step
    std::string file = dtk::rep_str(core_prtcl_loc,"${step}",step);
    CorePrtcls cp;
    dtk::read_gio_quick(file,"fof_halo_tag",cp.htag,cp.size);
    dtk::read_gio_quick(file,"id",          cp.pid,cp.size);
    int64_t* srt = dtk::new_sort_array<int64_t>(cp.size);
    std::sort(srt,srt+cp.size,dtk::Sorter<int64_t>(cp.htag));
    dtk::reorder(cp.htag,cp.size,srt);
    dtk::reorder(cp.pid,cp.size,srt);
    delete [] srt;
    all_core_prtcls[step]=cp;
    t.stop();
    std::cout<<t<<std::endl;
  }
}
void load_accum_prtcls(int step){
  if(all_accum_prtcls.count(step)==0){
    dtk::Timer t;t.start();
    std::cout<<"Loading accum step "<<step<<"...";
    //not loaded yet -> lets load that step
    std::string file = dtk::rep_str(accum_prtcl_loc,"${step}",step);
    AccumPrtcls ap;
    dtk::read_gio_quick(file,"id", ap.pid,ap.size);
    dtk::read_gio_quick(file,"x",  ap.x,  ap.size);
    dtk::read_gio_quick(file,"y",  ap.y,  ap.size);
    dtk::read_gio_quick(file,"z",  ap.z,  ap.size);
    int64_t* srt = dtk::new_sort_array<int64_t>(ap.size);
    std::sort(srt,srt+ap.size,dtk::Sorter<int64_t>(ap.pid));
    dtk::reorder(ap.pid,ap.size,srt);
    dtk::reorder(ap.x,ap.size,srt);
    dtk::reorder(ap.y,ap.size,srt);
    dtk::reorder(ap.z,ap.size,srt);
    delete [] srt;
    all_accum_prtcls[step]=ap;
    t.stop();
    std::cout<<t<<std::endl;
  }
}
void load_halo_cats(){
  dtk::Timer t;t.start();
  std::cout<<"Loading halos"<<std::endl;
  for(int i=0;i<steps.size();++i){
    //int64_t fof_size,sod_size;
    ///int64_t* fof_htag,sod_htag;
    //float* x,y,z,fof_mass,sod_mass,sod_radius;
    //std::string fof_file = dtk::rep_str(fof_loc,"${step}",step[i]);
    HaloCat hc;
    std::string sod_file = dtk::rep_str(sod_loc,"${step}",steps[i]);
    // If 
    if(sod_hdf5){
      std::cout<<sod_file<<std::endl;
      std::cout<<hc.htag<<" "<<hc.size<<std::endl;
      dtk::read_hdf5(sod_file, "fof_halo_tag",      hc.htag,       hc.size, true); 
      dtk::read_hdf5(sod_file, "fof_halo_center_x", hc.x,          hc.size, true);
      dtk::read_hdf5(sod_file, "fof_halo_center_y", hc.y,          hc.size, true);
      dtk::read_hdf5(sod_file, "fof_halo_center_z", hc.z,          hc.size, true);
      dtk::read_hdf5(sod_file, "sod_halo_mass_m200m",     hc.sod_mass,   hc.size, true);
      dtk::read_hdf5(sod_file, "sod_halo_radius_r200m",   hc.sod_radius, hc.size, true);
    }
    else{
      dtk::read_gio_quick(sod_file,"fof_halo_tag",      hc.htag,       hc.size);
      dtk::read_gio_quick(sod_file,"fof_halo_center_x", hc.x,          hc.size);
      dtk::read_gio_quick(sod_file,"fof_halo_center_y", hc.y,          hc.size);
      dtk::read_gio_quick(sod_file,"fof_halo_center_z", hc.z,          hc.size);
      dtk::read_gio_quick(sod_file,"sod_halo_mass",     hc.sod_mass,   hc.size);
      dtk::read_gio_quick(sod_file,"sod_halo_radius",   hc.sod_radius, hc.size);
    }
    all_halocats[steps[i]]=hc;
  }
  t.stop();
  std::cout<<"\tdone. time: "<<t<<std::endl;
}
void load_zmr_sdss(std::string zmr_loc,ZMR& zmr_sdss){
  if(rank==0)
    std::cout<<"Loading redmapper zmr..."<<std::endl;
  dtk::AutoTimer t;
  if(rank == 0){
    zmr_sdss.load(zmr_loc);
  }
  broadcast_zmr(zmr_sdss);
  /*  for(int i =0;i<nproc;++i){
    if(i==rank){
      std::cout<<"Rank: "<<i<<std::endl;
      zmr_sdss.print_non_zero();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    }*/
  zmr_cores.copy_bins(zmr_sdss);
  if(rank==0)
    std::cout<<"\tdone. time: "<<t<<std::endl;
}
void broadcast_zmr(ZMR& zmr){
  MPI_Bcast(&zmr.z_size,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&zmr.m_size,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&zmr.r_size,1,MPI_INT,0,MPI_COMM_WORLD);
  if(rank != 0){
    zmr.alloc(zmr.z_size,zmr.m_size,zmr.r_size);
    zmr.z_bins.resize(zmr.z_size+1);
    zmr.m_bins.resize(zmr.m_size+1);
    zmr.r_bins.resize(zmr.r_size+1);
  }
  MPI_Bcast(&zmr.z_bins[0], zmr.z_bins.size(), MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&zmr.m_bins[0], zmr.m_bins.size(), MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&zmr.r_bins[0], zmr.r_bins.size(), MPI_FLOAT,0,MPI_COMM_WORLD);
  
  MPI_Bcast(&zmr.zm_Ngal[0],     zmr.zm_Ngal.size(),     MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&zmr.zm_Ngal_err[0], zmr.zm_Ngal_err.size(), MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&zmr.zm_Ngal_var[0], zmr.zm_Ngal_var.size(), MPI_FLOAT,0,MPI_COMM_WORLD);

  MPI_Bcast(&zmr.zmr_gal_density[0],     zmr.zmr_gal_density.size(),     MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&zmr.zmr_gal_density_err[0], zmr.zmr_gal_density_err.size(), MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&zmr.zmr_gal_density_var[0], zmr.zmr_gal_density_var.size(), MPI_FLOAT,0,MPI_COMM_WORLD);

  MPI_Bcast(&zmr.zmr_gal_accum[0],     zmr.zmr_gal_accum.size(),     MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&zmr.zmr_gal_accum_err[0], zmr.zmr_gal_accum_err.size(), MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&zmr.zmr_gal_accum_var[0], zmr.zmr_gal_accum_var.size(), MPI_FLOAT,0,MPI_COMM_WORLD);

  MPI_Bcast(&zmr.zmr_dgal_dr[0],     zmr.zmr_dgal_dr.size(),         MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&zmr.zmr_dgal_dr_err[0], zmr.zmr_dgal_dr_err.size(),     MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&zmr.zmr_dgal_dr_var[0], zmr.zmr_dgal_dr_var.size(),     MPI_FLOAT,0,MPI_COMM_WORLD);

  MPI_Bcast(&zmr.zm_counts[0],  zmr.zm_counts.size(),  MPI_INT64_T,0,MPI_COMM_WORLD);
  MPI_Bcast(&zmr.zmr_counts[0], zmr.zmr_counts.size(), MPI_INT64_T,0,MPI_COMM_WORLD);
}
void make_clusters(){
  std::cout<<"Making clusters"<<std::endl;
  dtk::Timer t;t.start();
  for(int i =0;i<steps.size();++i){
    std::cout<<"\tstep="<<steps[i]<<std::endl;
    for(int j =0;j< all_halocats.at(steps[i]).size;++j){
      //std::cout<<all_halocats.at(i).htag[j]<<" "<<all_halocats.at(i).sod_mass[j]<<" "<<all_halocats.at(i).sod_radius[j]<<" " <<all_halocats.at(i).x[j]<<" " <<all_halocats.at(i).y[j]<<" "<<all_halocats.at(i).z[j]<<" "<<steps.at(i)<<std::endl;
      if(all_halocats.at(steps[i]).sod_mass[j] >1e14){
	dtk::AutoTimer t;
	all_clusters.push_back(Cluster(all_halocats.at(steps[i]).htag[j],
				       all_halocats.at(steps[i]).sod_mass[j],
				       all_halocats.at(steps[i]).sod_radius[j],
				       all_halocats.at(steps[i]).x[j],
				       all_halocats.at(steps[i]).y[j],
				       all_halocats.at(steps[i]).z[j],
				       steps.at(i)));
			       //all_clusters.push_back(clstr);
	//break; //only one cluster
	std::cout<<"\tcluster["<<j<<"]/"<<all_halocats.at(steps[i]).size<<std::endl;
	std::cout<<"\t"<<dtk::fraction(j,all_halocats.at(steps[i]).size)<<std::endl;
	std::cout<<"\ttime: "<<t<<std::endl;
      }
      if(all_clusters.size() > 5000){
	break;
      }
    }
  }
  dv.sort();
  t.stop();
  std::cout<<"done making clusters. Time: "<<t<<std::endl;
}
void defrag_halos(int64_t *htags, int64_t size){
  for(int64_t i =0;i<size;++i){
    //    bool print = htags[i]<0;
    //if(print)
    //  std::cout<<htags[i]<<" -> ";
    
    int64_t htag = std::abs(htags[i]); //get rid of negative sign
    htags[i] = htag & 0x0000FFFFffffFFFF;
    //if(htags[i] == 535735191)
    //std::cout<<"it's fucking here!!!"<<std::endl;
    //if(print)
    //  std::cout<<htag<<std::endl;
  }
}
void move_together(float x0,float rL, float* x, int64_t size){
  for(int64_t i=0;i<size;++i){
    float diff = x0-x[i];
    if(diff > rL/2.0)
      x[i] += rL;
    else if(diff < -rL/2.0)
      x[i] -= rL;
  }
}
void make_zmr(const std::vector<Cluster>& clstrs, float m_infall,float r_disrupt,float r_fof,float r_merger,
	      ZMR& zmr,bool verbose){
  zmr.check_same_bins(zmr_sdss);
  if(verbose)
    std::cout<<"Making a zmr "<<m_infall<<" "<<r_disrupt<<" "<<r_fof<<" "<<r_merger<<std::endl;
  dtk::Timer t;t.start();
  std::vector<float> r2_bins(zmr.r_bins.size()); //calculate the squared bin size to avoid the sqrt in the 
  for(size_t i =0;i<r2_bins.size();++i){  //radial bin finding. 
    r2_bins[i]=zmr.r_bins[i]*zmr.r_bins[i];
  }
  zmr.zero_out();
  float Ngal;
  std::vector<float> r_cnt(zmr.r_size);
  if(verbose)
    std::cout<<"Makign zmr"<<std::endl;
  //#pragma omp parallel for
  for(int i =0;i<clstrs.size();++i){
    //zero out the values;
    Ngal = 0;
    for(int j =0;j<r_cnt.size();++j)
      r_cnt[j]=0.0;
    int z_i = dtk::find_bin(zmr.z_bins,clstrs[i].redshift);
    int m_i = dtk::find_bin(zmr.m_bins,clstrs[i].sod_mass);
    // std::cout<<"z_i: "<<z_i<<"  m_i: "<<m_i<<std::endl;
    float Ngal_r2_lim =1.0; //galaxies within sqrt(1.0)*r200 are counted for Ngal;
    Galaxies gal;
    // std::cout<<"Printing cluster info"<<std::endl;
    // clstrs[i].print_prop();
    // clstrs[i].print_cores();
    gal.reserve(clstrs[i].core_size);
    bool massive = clstrs[i].sod_mass > 1.1e+14;
    //std::cout<<clstrs[i].sod_mass<<std::endl;
    if(fit_r_merger){
      if(verbose && massive && false){
	std::cout<<"mass: "<<clstrs[i].sod_mass<<std::endl;
	std::cout<<"we are merging shit!"<<std::endl;
	std::cout<<r_merger<<std::endl;
      }
      clstrs[i].get_compact_merged_galaxies(m_infall,r_disrupt,r_merger,gal, verbose=false);
    }
    else{
      // std::cout<<"Standard compact cores"<<std::endl;
      clstrs[i].get_compact_galaxies(m_infall,r_disrupt,gal);
    }
    if(fit_r_fof){
      //dtk::Timer t1,t2;t1.start();
      clstrs[i].get_fof_galaxies(m_infall,r_disrupt,r_fof,gal);
      //t1.stop();t2.start();
      //clstrs[i].get_fof_galaxies2(m_infall,r_disrupt,r_fof,gal);
      //t2.stop();
      //std::cout<<"old: "<<t1<<" new: "<<t2<<std::endl;
    }

    //clstrs[i].get_radial_bins(m_infall,r_disrupt,r_fof,r2_bins,Ngal_r2_lim,r_cnt,Ngal);
    //add the galaxies to the zmr by binning them radially. x
    clstrs[i].get_radial_bins(gal,zmr.r_bins,r2_bins,Ngal_r2_lim,r_cnt,Ngal);
    //#pragma omp critcal
    {
      zmr.add_cluster(z_i,m_i,r_cnt,Ngal);
    }
    //gal.print(10);
  }
  zmr.finalize_data();

  t.stop();
  if(verbose){
    //zmr.print_non_zero();
    std::cout<<"\tdone making zmr. time: "<<t<<std::endl;
  }
}
void make_zmr(const std::vector<Cluster>& clstrs,CoreParam cp, ZMR& zmr,bool verbose){
  make_zmr(clstrs,cp.m_infall,cp.r_disrupt,cp.r_fof,cp.r_merger,zmr,verbose);
}

double calc_diff(const ZMR& zmr1,const  ZMR& zmr2,const std::map<int,Cores>& all_cores,const CoreParam& cp){
  double cost=0.0;
  if(!zmr1.check_same_bins(zmr2)){
    std::cout<<"Different bins?"<<std::endl;
    throw;
  }
  //All the costs are calculated as X^2/goodness of fit. 
  if(cost_gal_density)
    cost +=calc_diff_gal_density(zmr1,zmr2);
  if(cost_dgal_dr)
    cost +=calc_diff_dgal_dr(zmr1,zmr2);
  if(cost_gal_accum)
    cost +=calc_diff_gal_accum(zmr1,zmr2);
  if(cost_abundance)
    cost +=calc_diff_abundance(all_cores,cp);
  //return need a factor of 1/2 to go from X^2->log likelihood
  return cost/2.0;
}
double calc_diff_gal_density(const ZMR& zmr1, const ZMR& zmr2){
  //  std::cout<<"gal density"<<std::endl;
  double cost=0.0;
  int z_size = zmr1.z_size;
  int m_size = zmr1.m_size;
  int r_size = zmr1.r_size;
  for(int z_i=0;z_i<z_size;++z_i)
    for(int m_i=0;m_i<m_size;++m_i){
      int indx = zmr1.index(z_i,m_i);
      if(zmr1.zm_counts[indx] > cost_min_cluster_num && zmr2.zm_counts[indx] > cost_min_cluster_num){ //only calculate the difference where
	for(int r_i=0;r_i<r_size;++r_i){//both zmrs have clusters. 
	  int indx = zmr1.index(z_i,m_i,r_i);
	  float diff = std::fabs(zmr1.zmr_gal_density[indx]-zmr2.zmr_gal_density[indx]);
	  diff = diff*diff;
	  float err1 = zmr1.zmr_gal_density_err[indx];
	  float err2 = zmr2.zmr_gal_density_err[indx];
	  float err_sq = err1*err1 + err2*err2;
	  cost+= diff/err_sq;
	  if(std::isnan(cost)){
	    std::cout<<"is nan"<<std::endl;
	    std::cout<<"diff_sq"<<diff<<std::endl;
	    std::cout<<"err1"<<err1<<std::endl;
	    std::cout<<"err2"<<err2<<std::endl;
	    std::cout<<"err_sq"<<std::endl;
	    dtk::pause();
	  }
	}
      }
    }
  return cost;
}
double calc_diff_dgal_dr(const ZMR& zmr1, const ZMR& zmr2){
  //std::cout<<"dgal dr"<<std::endl;
  double cost=0.0;
  int z_size = zmr1.z_size;
  int m_size = zmr1.m_size;
  int r_size = zmr1.r_size;
  for(int z_i=0;z_i<z_size;++z_i)
    for(int m_i=0;m_i<m_size;++m_i){
      int indx = zmr1.index(z_i,m_i);
      if(zmr1.zm_counts[indx] > cost_min_cluster_num && zmr2.zm_counts[indx] > cost_min_cluster_num){ //only calculate the difference where
	for(int r_i=0;r_i<r_size;++r_i){//both zmrs have clusters. 
	  int indx = zmr1.index(z_i,m_i,r_i);
	  float diff = std::fabs(zmr1.zmr_dgal_dr[indx]-zmr2.zmr_dgal_dr[indx]);
	  diff = diff*diff;
	  float err1 = zmr1.zmr_dgal_dr_err[indx];
	  float err2 = zmr2.zmr_dgal_dr_err[indx];
	  float err_sq = err1*err1 + err2*err2;
	  cost+= diff/err_sq;
	}
      }
    }
  return cost;
}
double calc_diff_gal_accum(const ZMR& zmr1, const ZMR& zmr2){
  std::cout<<"gal accum"<<std::endl;
  double cost=0.0;
  int z_size = zmr1.z_size;
  int m_size = zmr1.m_size;
  int r_size = zmr1.r_size;
  for(int z_i=0;z_i<z_size;++z_i)
    for(int m_i=0;m_i<m_size;++m_i){
      int indx = zmr1.index(z_i,m_i);
      if(zmr1.zm_counts[indx] > cost_min_cluster_num && zmr2.zm_counts[indx] > cost_min_cluster_num){ //only calculate the difference where
	for(int r_i=0;r_i<r_size;++r_i){//both zmrs have clusters. 
	  int indx = zmr1.index(z_i,m_i,r_i);
	  float diff = std::fabs(zmr1.zmr_gal_accum[indx]-zmr2.zmr_gal_accum[indx]);
	  diff = diff*diff;
	  float err1 = zmr1.zmr_gal_accum_err[indx];
	  float err2 = zmr2.zmr_gal_accum_err[indx];
	  float err_sq = err1*err1 + err2*err2;
	  cost+= diff/err_sq;
	}
      }
    }
  return cost;
}
double calc_diff_abundance(const std::map<int,Cores>& all_cores,const CoreParam& cp){
  double cost = 0.0;
  float box_vol = rL*rL*rL;
  float expected_sum = expected_comov_abundance*box_vol; 
  // std::cout<<"parameters: "<<cp.m_infall<<" "<<cp.r_disrupt<<std::endl;  
  // std::cout<<"expected_sum: "<<expected_sum<<std::endl;
  //  dtk::pause();
  for(int i =0;i<steps.size();++i){
    const Cores cores = all_cores.at(steps[i]);
    int* compact = new int[cores.size];
    int ssum;
    for(int j =0;j<cores.size;++j){
      compact[j]=(cores.radius[j]<cp.r_disrupt) && (cores.infall_mass[j]>cp.m_infall);
    }		
    float offset =0;
    float sum = std::accumulate(compact,compact+cores.size,offset);
    // std::cout<<"sum:"<<sum<<std::endl;
    float sum_err = std::sqrt(sum+1.0);
    // std::cout<<"sum_err: "<<sum_err<<std::endl;
    // std::cout<<"diff: "<<expected_sum-sum<<std::endl;
    float err = (expected_sum-sum)/sum_err;
    // std::cout<<"err: "<<err<<std::endl;
    // std::cout<<"err^2: "<<err*err<<std::endl;
    cost += err*err;
    delete [] compact;
  }
  //std::cout<<"cost: "<<cost<<std::endl;
  // std::cout<<"\n"<<std::endl;
  // dtk::pause();
  return cost;
}

void calculate_likelihood_grid(std::vector<float> mi_bins, std::vector<float> rd_bins,
			       std::vector<float> rm_bins, std::string       out_loc, 
			       CoreParam& max_lkhd_param){
  std::cout<<"Calculating likelihoods..."<<std::endl;
  std::cout<<"bins size: "<<mi_bins.size()<<" "<<rd_bins.size()<<" "<<rm_bins.size()<<std::endl;
  dtk::AutoTimer t;
  std::vector<double> result(mi_bins.size()*rd_bins.size()*rm_bins.size());
  uint indx =0;
  dtk::AutoTimer t2;
#pragma omp parallel for firstprivate(zmr_cores, t2), private(indx), schedule(dynamic)
  for(uint mi=0;mi<mi_bins.size();++mi){
    zmr_cores.copy_bins(zmr_sdss);
    // if(omp_get_thread_num()==0)
    //   std::cout<<"\r\t"<<mi<<"/"<<mi_bins.size()<<"\t*/"<<rd_bins.size()<<" time: "<<t2;
    t2.start();
    for(uint rd=0;rd<rd_bins.size();++rd){
      // if(omp_get_thread_num()==0)
      // 	std::cout<<"\r\t"<<mi<<"/"<<mi_bins.size()<<"\t"<<rd<<"/"<<rd_bins.size();
      float r_dspt_lock_val;
      if(lock_r_disrupt)
	r_dspt_lock_val= get_locked_r_disrupt(mi_bins.at(mi),all_cores[401]);
      for(uint rm=0;rm<rm_bins.size();++rm){
	if(omp_get_thread_num()==0)
	  std::cout<<"\r\t"<<mi<<"/"<<mi_bins.size()<<"\t"<<rd<<"/"<<rd_bins.size()<<"\t"<<rm<<"/"<<rm_bins.size();

	indx = rm + rd*rm_bins.size() + mi*rm_bins.size()*rd_bins.size();
	if(!lock_r_disrupt){
	  make_zmr(all_clusters,mi_bins.at(mi),rd_bins.at(rd),0.0,rm_bins.at(rm),zmr_cores,false);
	}
	else{
	  make_zmr(all_clusters,mi_bins.at(mi),r_dspt_lock_val,0.0,rm_bins.at(rm),zmr_cores,false);
	}
	CoreParam cp;
	cp.m_infall = mi_bins[mi];
	cp.r_disrupt=rd_bins[rd];
	cp.r_merger= rm_bins[rm];
	cp.r_fof = 0.0;
	float log_lkhd;
	//log_lkhd= fake_log_likelihood(cp);
	log_lkhd = calc_diff(zmr_sdss,zmr_cores,all_cores,cp);
	result.at(indx)=log_lkhd;
	//++indx;
      }
    }
  }
  int max_lkhd_indx = 0;
  float max_lkhd_val = result[max_lkhd_indx];
  for(int i =0; i< result.size();++i){
    if(result[i]<max_lkhd_val){
      max_lkhd_indx=i;
      max_lkhd_val = result[i];
    }
  }
  int mi_max_indx, rd_max_indx, rm_max_indx;
  mi_max_indx = max_lkhd_indx/(rm_bins.size()*rd_bins.size());
  rd_max_indx = (max_lkhd_indx/rm_bins.size())%rd_bins.size();
  rm_max_indx = max_lkhd_indx%rm_bins.size();
  max_lkhd_param.m_infall  = mi_bins[mi_max_indx];
  max_lkhd_param.r_disrupt = rd_bins[rd_max_indx];
  max_lkhd_param.r_merger  = rm_bins[rm_max_indx];
  max_lkhd_param.r_fof     = -1.0;
  std::cout<<std::endl;
  std::cout<<"log likelihood: "<<max_lkhd_val<<std::endl;
  std::cout<<"m_infall best: "<<mi_bins[mi_max_indx]<<std::endl;
  std::cout<<"r_disrup best: "<<rd_bins[rd_max_indx]<<std::endl;
  std::cout<<"r_merger best: "<<rm_bins[rm_max_indx]<<std::endl;
  std::ofstream file(out_loc.c_str());
  write_array(file,"result",result);
  write_array(file,"mi_bins",mi_bins);
  write_array(file,"rd_bins",rd_bins);  
  write_array(file,"rm_bins",rm_bins);
  file.close();
  std::cout<<"done time: "<<t<<std::endl;
}
double param_cost(double m_i,double r_d,double r_m,double r_fof){
  double massive_cost = 100000;
  double result = 0.0;
  if(m_i < m_infall_constraints[0] || m_i > m_infall_constraints[1]){
    result += massive_cost;
    float m_under = m_i - m_infall_constraints[0];
    float m_over =  m_i - m_infall_constraints[1];
    if(m_under < 0)
      result += -std::log10(m_i)*100;
    if(m_over > 0)
      result += std::log10(m_i)*100;
  }
  if(r_d < r_disrupt_constraints[0] || r_d > r_disrupt_constraints[1])
    result += massive_cost;
  if(fit_r_merger && (r_m < r_merger_constraints[0] || r_m > r_merger_constraints[1])){
    result += massive_cost; //only use constraints if this is an active variable

  }
  if(fit_r_fof && (r_fof < r_fof_constraints[0] || r_fof> r_fof_constraints[1]))
    result += massive_cost;//only use constraints if this is an active variable
  return result;
}
double param_cost(CoreParam cp){
  param_cost(cp.m_infall,cp.r_disrupt,cp.r_merger,cp.r_fof);
}
double central_cost(ZMR& zmr){
  return (zmr.centrals_expected - zmr.centrals_found)*cost_missed_centrals;

}
double min_func_gsl(const gsl_vector *x, void* param){
  FunctionParam* fparam = (FunctionParam*) param;
  double m_i,r_d,r_m,r_fof;
  int indx=0;
  m_i = std::pow(10,gsl_vector_get(x,indx++));
  if(!lock_r_disrupt)
    r_d = std::pow(10,gsl_vector_get(x,indx++));
  else
    r_d = get_locked_r_disrupt(m_i,all_cores[step]);
  //m_i = dv.get_m_infall(gsl_vector_get(x,indx++));
  //  r_d = dv.get_r_disrupt(gsl_vector_get(x,indx++));
  //std::cout<<"m_i: "<<m_i<<" r_d: "<<r_d<<std::endl;
  if(fit_r_fof){
    r_fof = std::pow(10,gsl_vector_get(x,indx++));
  }
  if(fit_r_merger){
    r_fof = std::pow(10,gsl_vector_get(x,indx++));
  }
  double cost = param_cost(m_i,r_d,r_m,r_fof);
  make_zmr(*(fparam->clusters),m_i,r_d,r_fof,r_m,*(fparam->zmr_cores),false);
  cost += central_cost(*(fparam->zmr_cores));
  CoreParam cp;
  cp.m_infall  = m_i;
  cp.r_disrupt = r_d;
  return cost+calc_diff(*(fparam->zmr_cores),*(fparam->zmr_sdss),all_cores,cp);//fake_log_likelihood(cp);//

}

CoreParam find_min(std::vector<Cluster>& clstrs, ZMR& zmr_cores, const ZMR& zmr_sdss,
		   float mi0,float rd0,float rfof0, float rm0){
  std::cout<<"Finding min..."<<std::endl;
  dtk::AutoTimer t;
  CoreParam result;
  FunctionParam fparam; //the constant parameters of the function
  fparam.zmr_cores = &zmr_cores;
  fparam.zmr_sdss = &zmr_sdss;
  fparam.clusters = &clstrs;
  std::cout<<"the vals are :"<<mi0<<" "<<rd0<<" "<<rfof0<<" "<<rm0<<std::endl;
  int indx=0;
  gsl_vector* x, *step_size; //the variable values in teh function
  x = gsl_vector_alloc(fit_var_num);
  gsl_vector_set(x,indx++,mi0);
  if(!lock_r_disrupt)
    gsl_vector_set(x,indx++,rd0);
  if(fit_r_fof)
    gsl_vector_set(x,indx++,rfof0);
  if(fit_r_merger)
    gsl_vector_set(x,indx++,rm0);
  std::cout<<fit_var_num<<std::endl;
  for(int i =0;i<fit_var_num;++i){
    std::cout<<"["<<i<<"]: "<<gsl_vector_get(x,i)<<std::endl;
    
  }
  gsl_multimin_function gmf; //the struct describing the functoin to minimize.
  gmf.f = min_func_gsl; //the function to minimize
  gmf.n = fit_var_num; //the number of variables to search over in x.
  gmf.params = (void*) &fparam; //the constant params of the function
  const gsl_multimin_fminimizer_type* T = gsl_multimin_fminimizer_nmsimplex2rand;
  gsl_multimin_fminimizer* state= NULL;
  size_t iter=0;
  int status;
  double size;
  step_size = gsl_vector_alloc(fit_var_num);  
  gsl_vector_set_all(step_size,simplex_step_size);
  gsl_vector_set(step_size,0,simplex_step_size*4.0);
  state = gsl_multimin_fminimizer_alloc(T,fit_var_num);
  gsl_multimin_fminimizer_set(state,&gmf,x,step_size);
  dtk::Timer t2;
  do{
    t2.start();
    ++iter;
    if(iter>max_iterations)
      break;

    status = gsl_multimin_fminimizer_iterate(state);
    if(status)
      break; //some failure...
    size = gsl_multimin_fminimizer_size(state);
    status = gsl_multimin_test_size(size,simplex_cnvrg_size);
    if(status == GSL_SUCCESS){
      std::cout<<"Converged!!\n"<<std::endl;
    }
    indx=0;
    t2.stop();
    std::cout<<"\n["<<iter<<"] cost= "<<state->fval
	     <<" time= "<<t2<<"\n\t\t "
	     <<gsl_vector_get(state->x,indx++)<<" ";
    if(!lock_r_disrupt)
      std::cout<<gsl_vector_get(state->x,indx++)<<" ";
    if(fit_r_fof)
      std::cout<<gsl_vector_get(state->x,indx++)<<" ";
    if(fit_r_merger)
      std::cout<<gsl_vector_get(state->x,indx++)<<" ";
      std::cout<<std::endl;
  }while(status==GSL_CONTINUE && iter<max_iterations);
  indx = 0;
  result.m_infall  = std::pow(10,gsl_vector_get(state->x,indx++));
  if(!lock_r_disrupt)
    result.r_disrupt = std::pow(10,gsl_vector_get(state->x,indx++));
  if(fit_r_fof)
    result.r_fof  = std::pow(10,gsl_vector_get(state->x,indx++));
  else
    result.r_fof = -1;
  if(fit_r_merger)
    result.r_merger  = std::pow(10,gsl_vector_get(state->x,indx++));
  else
    result.r_merger = -1;
  gsl_vector_free(x);
  gsl_vector_free(step_size);
  gsl_multimin_fminimizer_free(state);
  std::cout<<"done minimizing time: "<<t<<std::endl;
  return result;
}
float rand_float(){
  return (float) rand()/RAND_MAX;
}
void write_out_clusters(std::vector<Cluster>& clstrs, CoreParam cp){
  float m_infall = cp.m_infall;
  float r_disrupt = cp.r_disrupt;
  float r_merger  = cp.r_merger;
  float r_fof     = cp.r_fof;
  for(int i=0;i<clstrs.size();++i){
    Galaxies gal;
    if(fit_r_merger)
      clstrs[i].get_compact_merged_galaxies(m_infall,r_disrupt,r_merger,gal);
    else
      clstrs[i].get_compact_galaxies(m_infall,r_disrupt,gal);
    if(fit_r_fof)
      clstrs[i].get_fof_galaxies(m_infall,r_disrupt,r_fof,gal);

    std::stringstream ss;
    ss<<"tmp/clstr"<<i<<".param";
    clstrs[i].write_out(ss.str(),cp,gal);
  }
}
int glob_cnt=0;
float move_together(float dx, float rL){
  if(dx > rL/2.0)
    return dx - rL;
  if(dx < -rL/2.0)
    return dx + rL;
  return dx;
}
Cluster::Cluster(int64_t htag,float sod_mass,float sod_radius,float x,float y, float z,int step):
  htag(htag),sod_mass(sod_mass),sod_radius(sod_radius),redshift(stepz.z_from_step(step)),
  x(x),y(y),z(z),step(step),core_size(0),cp_size(0){
  //Q: why is sod_radius scaled this way??
  //A: it looks like this sod is only used for scaled_cluster_radial selection of cores. 
  sod_radius = sod_radius/(1.0/(redshift+1.0));
  //find the cores that belong to this halo/clusters
  std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
  std::cout<<"halo "<<htag<<" mass: "<<sod_mass<<std::endl;
  Cores corecat = all_cores[step]; std::cout<<__LINE__<<std::endl;
  dtk::ChainingMeshIndex& cmi = all_core_cms.at(step);std::cout<<__LINE__<<std::endl;
  std::vector<size_t> my_cores;std::cout<<__LINE__<<std::endl;
  float halo_pos[3] = {x,y,z};std::cout<<__LINE__<<std::endl;
  // Don't interate over all the cores, only the ones within the box
  if(false){
    for(int i =0;i<corecat.size;++i){
      float dx = fabs(move_together(corecat.x[i] - x,rL));
      float dy = fabs(move_together(corecat.y[i] - y,rL));
      float dz = fabs(move_together(corecat.z[i] - z,rL));
      float dr = sqrt(dx*dx + dy*dy + dz*dz);
      if(scaled_cluster_radial_volume && dr<scaled_cluster_radial_factor*sod_radius){
	my_cores.push_back(i);
      }
      else if(dr<cluster_radial_volume){
	//std::cout<<(corecat.host_htag[i]==htag)<<" "<<corecat.host_htag[i]<<" "<<dr<<" "<<std::endl;
	my_cores.push_back(i);
      }
    }
  }
  else if(false){
    //Old method of getting galaxy cores. Only checks local cell
    size_t central_cell = cmi.get_cell_id_from_position(halo_pos);
    size_t* cell_element;
    size_t cell_element_size;
    cmi.get_cell_element_indexes(central_cell, cell_element, cell_element_size);
    for(int i=0;i<cell_element_size;++i){
      size_t index = cell_element[i];
      float dx = fabs(move_together(corecat.x[index] - x, rL));
      float dy = fabs(move_together(corecat.y[index] - y, rL));
      float dz = fabs(move_together(corecat.z[index] - z, rL));
      float dr = sqrt(dx*dx + dy*dy + dz*dz);
      if(scaled_cluster_radial_volume && dr<scaled_cluster_radial_factor*sod_radius){
	my_cores.push_back(index);
      }
      else if(dr<cluster_radial_volume){
	//std::cout<<(corecat.host_htag[i]==htag)<<" "<<corecat.host_htag[i]<<" "<<dr<<" "<<std::endl;
	my_cores.push_back(index);
      }
    }
  }
  else{
    // Fastest & correct method of getting cores
    std::cout<<__LINE__<<std::endl;
    std::cout<<cluster_radial_volume<<std::endl;
    std::cout<<halo_pos<<std::endl;
    std::cout<<halo_pos[0]<<std::endl;
    my_cores  = cmi.query_elements_within(halo_pos, cluster_radial_volume);
  }
  std::cout<<"\tfound cores:"<<my_cores.size()<<std::endl;
  alloc_cores(my_cores.size());
  std::cout<<"\t alloc was fine.."<<std::endl;
  //copying the cores to clusters local info
  bool center_replaced =false;
  float central_mass = 0;
  for(int i =0;i<my_cores.size();++i){
    core_id[i]=corecat.ctag[my_cores[i]];
    core_htag[i]=corecat.infall_htag[my_cores[i]];
    core_step[i]=corecat.infall_step[my_cores[i]];
    core_x[i] =corecat.x[my_cores[i]];
    core_y[i] =corecat.y[my_cores[i]];
    core_z[i] =corecat.z[my_cores[i]];
    core_m[i] =corecat.infall_mass[my_cores[i]];
    core_r[i] =corecat.radius[my_cores[i]];
    core_is_central[i] =  corecat.is_central[my_cores[i]];
    dv.mass_infall.push_back(core_m[i]);
    dv.radius_disrupt.push_back(core_r[i]);
    // std::cout<<"\t\tcore_infall_mass: "<<core_m[i]<<" radius: "<<core_r[i]<<std::endl<<" radial_distance: ";
    // float dx = fabs(move_together(core_x[i] - x, rL));
    // float dy = fabs(move_together(core_y[i] - y, rL));
    // float dz = fabs(move_together(core_z[i] - z, rL));
    // float dr = sqrt(dx*dx + dy*dy + dz*dz);
    // std::cout<<dr<<std::endl;
    // std::cout<<core_x[i]<<" "<<core_y[i]<<" "<<core_z[i]<<std::endl;
    // std::cout<<x<<" "<<y<<" "<<z<<std::endl;
    if(force_center_on_central){
      if(corecat.is_central[my_cores[i]]==1 && corecat.host_htag[my_cores[i]] == htag &&
	 corecat.central_mass[my_cores[i]] > central_mass){
	central_mass = corecat.central_mass[my_cores[i]];
	center_replaced = true;
	std::cout<<"Central found!"<<std::endl;
	std::cout<<"\ncore id: "<<core_id[i]<<std::endl;
	std::cout<<"host htag: "<<corecat.host_htag[my_cores[i]]<<std::endl;
	std::cout<<"step: "<<core_step[i]<<std::endl;
	std::cout<<"pos: "<<core_x[i]<<" "<<core_y[i]<<" "<<core_z[i]<<std::endl;
	std::cout<<"mass: "<<std::log10(core_m[i])<<std::endl;
	std::cout<<"central: "<<core_is_central[i]<<std::endl;
	std::cout<<corecat.host_htag[my_cores[i]]<<"=="<<htag<<std::endl;
	x = core_x[i];
	y = core_y[i];
	z = core_z[i];
      }
    }
  }
  if(force_center_on_central && !center_replaced){
    std::cout<<"A central core wasn't found for this halo?!"<<std::endl;
    dtk::pause();
  }
  //Moving all the core particles across period boundary conditions to
  //be in one area
  move_together(x,rL,core_x,core_size);
  move_together(y,rL,core_y,core_size);
  move_together(z,rL,core_z,core_size);

  //If we don't need core particles, we don't need to load them :)
  if(no_core_particles)
    return;
  //getting the pids for each core
  std::vector<std::vector<int64_t> > cores_pids;
  for(int i =0;i<my_cores.size();++i){
    std::vector<int64_t> pids;
    pids.reserve(20);
    load_core_prtcls(core_step[i]);
    CorePrtcls& cprtcls = all_core_prtcls.at(core_step[i]);
    int indx = std::lower_bound(cprtcls.htag,cprtcls.htag+cprtcls.size,core_htag[i])-cprtcls.htag;
    while(indx<cprtcls.size){
      if(cprtcls.htag[indx]==core_htag[i]){
	pids.push_back(cprtcls.pid[indx]);
	++indx;
      }
      else
	break;
    }
    cores_pids.push_back(pids);
  }
  //getting the locations of the cores we want to update the position;
  std::vector<int> accum_indx;
  load_accum_prtcls(step);
  AccumPrtcls& ap =all_accum_prtcls[step];
  for(int i =0;i<my_cores.size();++i){

    int cnt = 0; //how many we found
    core_cp_offset[i]=accum_indx.size();

    for(int j =0;j<cores_pids.at(i).size();++j){
      int indx = std::lower_bound(ap.pid,ap.pid+ap.size,cores_pids[i][j])-ap.pid;
      if(ap.pid[indx]==cores_pids[i][j]){
	++cnt;
	accum_indx.push_back(indx);
      }
    }
    core_cp_size[i]=cnt;
  }
  std::cout<<"accum prtcls found: "<<accum_indx.size()<<std::endl;
  //actually recording the postions;
  cp_size = accum_indx.size();
  alloc_cp(accum_indx.size());
  for(int i =0;i<accum_indx.size();++i){
    cp_x[i] = ap.x[accum_indx[i]];
    cp_y[i] = ap.y[accum_indx[i]];
    cp_z[i] = ap.z[accum_indx[i]];
  }
  move_together(x,rL,cp_x,cp_size);
  move_together(y,rL,cp_y,cp_size);
  move_together(z,rL,cp_z,cp_size);
  //cp_cm.set_limits(cp_x,cp_y,cp_z,cp_size,cm_num,cm_num,cm_num);
  //cp_cm.set_data(cp_x,cp_y,cp_z,cp_size);
}
void Cluster::alloc_cores(int size){
  core_size = size;
  core_id = new int64_t[size];
  core_htag = new int64_t[size];
  core_step = new int[size];
  core_x = new float[size];
  core_y = new float[size];
  core_z = new float[size];
  core_r = new float[size];
  core_m = new float[size];
  core_is_central = new int[size];
  core_cp_offset = new int[size];
  core_cp_size = new int[size];
}
void Cluster::alloc_cp(int size){
  cp_x = new float[size];
  cp_y = new float[size];
  cp_z = new float[size];
  cp_usage = new bool[size];
  cp_color = new int[size];
}
void Cluster::write_out(std::string file_loc,CoreParam cp,Galaxies& gal){
  std::ofstream file(file_loc.c_str());
  file<<"halo_x "<<x<<std::endl;
  file<<"halo_y "<<y<<std::endl;
  file<<"halo_z "<<z<<std::endl;
  file<<"halo_mass "<<sod_mass<<std::endl;
  file<<"halo_radius "<<sod_radius<<std::endl;

  write_array(file,"core_x",core_x,core_size);
  write_array(file,"core_y",core_y,core_size);
  write_array(file,"core_z",core_z,core_size);
  write_array(file,"core_m",core_m,core_size);
  write_array(file,"core_r",core_r,core_size);

  write_array(file,"cprtcl_x",cp_x,cp_size);
  write_array(file,"cprtcl_y",cp_y,cp_size);
  write_array(file,"cprtcl_z",cp_z,cp_size);

  write_array(file,"gal_x",gal.x);
  write_array(file,"gal_y",gal.y);
  write_array(file,"gal_z",gal.z);
  write_array(file,"gal_w",gal.w);
  write_array(file,"gal_type",gal.type);

  write_array(file,"dis_cp_x",dis_cp_x);
  write_array(file,"dis_cp_y",dis_cp_y);
  write_array(file,"dis_cp_z",dis_cp_z);

  std::vector<float> r2_bins = zmr_sdss.r_bins;
  std::vector<float>   r_cnt(r2_bins.size()-1);
  float Ngal;
  for(int i =0;i<r2_bins.size();++i){
    r2_bins[i]=r2_bins[i]*r2_bins[i];
    r_cnt[i]=0;
  }
  get_radial_bins(gal,zmr_sdss.r_bins,r2_bins,1.0,r_cnt,Ngal);
  write_array(file,"r_bins",zmr_sdss.r_bins);
  write_array(file,"r_cnt",r_cnt);
  write_value(file,"Ngal",Ngal);
  write_value(file,"m_infall",cp.m_infall);
  write_value(file,"r_disrupt",cp.r_disrupt);
  write_value(file,"r_merger",cp.r_merger);
  write_value(file,"r_fof",cp.r_fof);
  
}
void Cluster::get_compact_galaxies(float m_infall, float r_disrupt, Galaxies& gal) const{
  // std::cout<<"====get_compact_galaxies===="<<m_infall<<"  "<<r_disrupt<<std::endl;
  for(int i =0;i<core_size;++i){
    // std::cout<<i<<" "<<core_m[i]<<" "<<core_r[i];
    if(core_m[i]>m_infall && (core_r[i]<r_disrupt || (force_central_as_galaxy && core_is_central[i]))){
      gal.x.push_back(core_x[i]);
      gal.y.push_back(core_y[i]);
      gal.z.push_back(core_z[i]);
      gal.w.push_back(1);
      gal.type.push_back(0);
      // std::cout<<" accepted"<<std::endl;
    }
    else{
      // std::cout<<" rejected"<<std::endl;
    }
  }
}
void Cluster::get_compact_merged_galaxies(float m_infall,float r_disrupt,float r_merger,
					  Galaxies& gal, bool verbose ) const{
  //ostd::cout<<"core merger core size: "<<core_size;
  std::vector<float> new_x,new_y,new_z;
  std::vector<int> new_w;
  for(int i =0;i<core_size;++i){
    if(core_m[i]>m_infall && (core_r[i]<r_disrupt || (force_central_as_galaxy && core_is_central[i]))){
      new_x.push_back(core_x[i]);
      new_y.push_back(core_y[i]);
      new_z.push_back(core_z[i]);
      new_w.push_back(0.0);
    }
  }
  int result_size=new_x.size();
  if(verbose){
    std::cout<<"m_infall: "<<m_infall<<std::endl;
    std::cout<<"r_disrupt: "<<r_disrupt<<std::endl;
    std::cout<<"r_merger: "<<r_merger<<std::endl<<std::endl;
    std::cout<<" n2 start size: "<<result_size;
  }
  n2_merger3d(&new_x[0],
	      &new_y[0],
	      &new_z[0],
	      &new_w[0],
	      &result_size,
	      r_merger,
	      NULL);
  if(verbose)
    std::cout<<" -> end size "<<result_size<<std::endl;
  for(int i=0;i<result_size;++i){
    gal.x.push_back(new_x[i]);
    gal.y.push_back(new_y[i]);
    gal.z.push_back(new_z[i]);
    gal.w.push_back(new_w[i]);
    //if(new_w[i]
    gal.type.push_back(0);
  }

}
void Cluster::get_fof_galaxies(float m_infall, float r_disrupt, float r_fof,Galaxies& gal) const{
  //Not implemented
  throw;
  for(int i =0;i<core_size;++i){
    bool use;
    if(core_m[i]>m_infall && core_r[i]>r_disrupt)
      use = true;
    else
      use = false;
    for(int j=0;j<core_cp_size[i];++j)
      cp_usage[core_cp_offset[i]+j]=use;
  }
  std::vector<float> new_x(cp_x,cp_x+cp_size);
  std::vector<float> new_y(cp_y,cp_y+cp_size);
  std::vector<float> new_z(cp_z,cp_z+cp_size);
  std::vector<int>   new_w(cp_size);

  int64_t* srt = dtk::new_sort_array<int64_t>(cp_size);
  std::sort(srt,srt+cp_size,dtk::Sorter<bool>(cp_usage));
  dtk::reorder(new_x,srt);
  dtk::reorder(new_y,srt);
  dtk::reorder(new_z,srt);
  dtk::reorder(cp_usage,cp_size,srt);
  int64_t indx = std::lower_bound(cp_usage,cp_usage+cp_size,true)-cp_usage;
  int result_size = cp_size-indx;
  //  std::cout<<"dis_cp size: "<<result_size<<std::endl;
  for(int i =0;i<result_size;++i){
    // dis_cp_x.push_back(new_x[indx+i]);
    // dis_cp_y.push_back(new_y[indx+i]);
    // dis_cp_z.push_back(new_z[indx+i]);
  }
  n2_merger3d(&new_x[indx],
	      &new_y[indx],
	      &new_z[indx],
	      &new_w[indx],
	      &result_size,
	      r_fof,
	      NULL);
  for(int i=0;i<result_size;++i){
    //if there are enough partilces in this fof group, 
    //call it a bcg like galaxy. 
    //std::cout<<i<<"/"<<result_size<<": "<<new_w[indx+i]<<std::endl;
    
    if(new_w[indx+i]>=min_disrupt_fof){
      //std::cout<<gal.x[indx+i]<<" "<<gal.y[indx+i]<<" "<<gal.z[indx+i]<<std::endl;
      gal.x.push_back(new_x[indx+i]);
      gal.y.push_back(new_y[indx+i]);
      gal.z.push_back(new_z[indx+i]);
      gal.w.push_back(0.0);
      gal.type.push_back(1);
    }
  }
  delete [] srt;
}
void Cluster::get_fof_galaxies2(float m_infall,float r_disrupt, float r_fof, Galaxies& gal) const{
  //Not implemented
  throw;
  for(int i =0;i<core_size;++i){
    bool use;
    if(core_m[i]>m_infall && core_r[i]>r_disrupt)
      use = true;
    else
      use = false;
    for(int j=0;j<core_cp_size[i];++j)
      cp_usage[core_cp_offset[i]+j]=use;
  }
  // dtk::reorder(cp_usage,cp_size,cp_cm.get_srt());
  // find_fof_cm(r_fof,cp_x,cp_y,cp_z,cp_usage,cp_color,cp_size,cp_cm);
  std::vector<std::vector<float> > fof_x,fof_y,fof_z;
  /*std::vector<int> reordered_colors(cp_color,cp_color+cp_size);
  int* srt = dtk::new_sort_array<int>(cp_size);
  std::sort(srt,srt+cp_size,dtk::Sorter<int>(cp_color));
  dtk::reorder(reordered_colors,srt);
  std::vector<size_t> group_start,group_size;
  dtk::group_by_id(reordered_colors,group_start,group_size);
  std::vector<float> fof_xs,fof_ys,fof_zs;
  for(int i =0;i<group_start.size();++i){
    int group_id = cp_color[srt[group_start[i]]];
    if(group_id == -1) //unused cores
      continue;
    for(int j =0;j<group_size[i];++j){
      fof_xs.push_back(cp_x[srt[group_start[i]+j]]);
      fof_ys.push_back(cp_y[srt[group_start[i]+j]]);
      fof_zs.push_back(cp_z[srt[group_start[i]+j]]);
    }
    if(fof_xs.size()>=min_disrupt_fof){
      gal.x.push_back(dtk::average(fof_xs));
      gal.y.push_back(dtk::average(fof_ys));
      gal.z.push_back(dtk::average(fof_zs));
      gal.w.push_back(fof_xs.size());
      gal.type.push_back(1);
    }
    fof_xs.clear();
    fof_ys.clear();
    fof_zs.clear();
    }
    delete [] srt;*/
}
void Cluster::get_radial_bins(float m_infall,float r_disrupt,float r_fof, 
			      std::vector<float>& r2_bins, float Ngal_r2_lim,
			      std::vector<float>& r_cnt, float& Ngal) const{
  //Note: r2_bins are the radial bin edges squared. <- so we can avoid an unneed square root.
  //same for Ngal_r2_lim.
  //Galaxies gal;
  //get_compact_galaxies(m_infall,r_disrupt,gal);
  //get_bcg_galaxies(m_infall,r_disrupt,r_fof,gal);
  //get_radial_bins(gal,r2_bins,Ngal_r2_lim,r_cnt,Ngal);
  throw; //we should not use this function
}
void Cluster::get_radial_bins(Galaxies& gal, std::vector<float>& r_bins,
			      std::vector<float>& r2_bins, float Ngal_r2_lim,
			      std::vector<float>& r_cnt, float& Ngal) const{
  float sqr_r200 = sod_radius*sod_radius;
  /*  for(int i =0;i<r2_bins.size();++i){
    std::cout<<i<<": "<<sqrt(r2_bins[i])<<std::endl;
    }*/
  for(int i =0;i<gal.x.size();++i){
    float dx = gal.x[i]-x;
    float dy = gal.y[i]-y;
    float dz = gal.z[i]-z;
    //TODO: Remove temp factor of 1.2
    float dr = sqrt(dx*dx + dy*dy + dz*dz)/(sod_radius); 
    
    if(false){//only one projection along z-axis
      float dr2 = dx*dx+dy*dy;
      int indx = dtk::find_bin(r2_bins,dr2);
      if(indx!=-1){
	r_cnt.at(indx) += 1.0;
      }
      if(dr2<Ngal_r2_lim){
	Ngal+=1.0;
      }
    }
    else if(false){ //various projection count
      double weight = 1.0/diff_projection_num;
      for(int j =0;j<diff_projection_num;++j){
	float theta = rand_float()*M_PI/2.0;
	float proj_dr = std::sin(theta)*dr;
	float proj_dr2 = proj_dr*proj_dr;
	int indx = dtk::find_bin(r2_bins,proj_dr2);
	if(indx !=-1){
	  r_cnt.at(indx) += weight;
	  //std::cout<<"i="<<indx<<"  "<<sqrt(r2_bins[indx])<<" < "<<sqrt(r2)<<" < "<<sqrt(r2_bins[indx+1])<<std::endl;
	}
	if(proj_dr2 < Ngal_r2_lim)
	  Ngal += weight;
      }
    }
    else if(false){//analytical projections //WRONG!!!!
      throw; //THis should never be used!!! I'm not deleting just to keep a 
      //record of it. 
      if(dr <= r_bins[0]) //the galaxy doesn't make it to any bin. 
	continue;
      float start_angle = std::asin(r_bins[0]/dr);
      //start_angle = 0;
      std::cout<<"radius: "<<dr<<std::endl;
      float sum = 0;
      for(int i =0;i<r_bins.size()-1;++i){
	std::cout<<"["<<i<<"]: "<<r_bins[i]<<"->"<<r_bins[i+1]<<std::endl;
	if(dr<r_bins[i+1]){
	  r_cnt[i] += (M_PI/2.0-start_angle)/(M_PI/2.0);
	  std::cout<<"["<<i<<"]b:"<<(M_PI/2.0-start_angle)/(M_PI/2.0)<<std::endl;;
	  sum +=(M_PI/2.0-start_angle)/(M_PI/2.0);
	  break; //the remainder of the projected angle goes int this bin
	}
	else{
	  float end_angle = std::asin(r_bins[i+1]/dr);
	  //std::cout<<"\t\tstart_angle: "<<start_angle<<std::endl;
	  //std::cout<<"\t\tend_angle: "<<end_angle<<std::endl;
	  r_cnt[i] += (end_angle-start_angle)/(M_PI/2.0);
	  std::cout<<"["<<i<<"]a:"<<(end_angle-start_angle)/(M_PI/2.0)<<std::endl;
	  sum +=(end_angle-start_angle)/(M_PI/2.0);
	  start_angle = end_angle;
	}
      }
      std::cout<<"sum: "<<sum<<std::endl;
    }
    //Correct code in this section
    else{
      if(dr < r_bins[0]){
	//doesn't even make into any radial bin
      }
      else{
	float start_angle;
	if(dr==0){
	  start_angle = 0;// M_PI*0.5; Why was this M_PI?
	}
	else{
	  start_angle= std::asin(r_bins[0]/dr);
	}
	//solid angle = 2pi(1-cos(theta))
	//theta 0->pi, in our case 0->pi/2
	//A half sphere is 2pi, but we will consider it to be 1
	//So the area for each radial bin is (area(outer bin edge) - area(inner bin edge) 
	// => (1-cos(out)) -(1-cos(in))
	// => -cos(out) + cos(in)
	// std::cout<<"radius: "<<dr<<std::endl; 
	// std::cout<<"start_angle: "<<start_angle<<std::endl;
	float sum =0;
	for(int i =0;i<r_bins.size()-1;++i){
	  // std::cout<<"["<<i<<"]: "<<r_bins[i]<<"->"<<r_bins[i+1]<<std::endl;
	  if(dr<r_bins[i+1]){
	    float weight = 0.0+std::cos(start_angle); //(area(pi/2) - area(start_angle))/tot_area
	    //float weight = -1.0+std::cos(start_angle); //(area(pi/2) - area(start_angle))/tot_area
	    r_cnt[i] = weight;
	    // std::cout<<"["<<i<<"]a:"<<weight<<std::endl;
	    // std::cout<<"start_angle: "<<start_angle<<std::endl;
	    // if(std::isnan(weight)){
	    //   std::cout<<"weight is nan a: "<<weight<<std::endl;
	    //   dtk::pause();
	    // }
	    // if(dr==0){
	    //   dtk::pause();
	    // }
	    sum+=weight;
	    break;
	  }
	  else{
	    float end_angle = std::asin(r_bins[i+1]/dr);
	    float weight = -std::cos(end_angle) +std::cos(start_angle); //(area(end)-area(start))/tot_area
	    // std::cout<<"["<<i<<"]b:"<<weight<<std::endl;
	    // if(std::isnan(weight)){
	    //   std::cout<<"weight is nan b: "<<weight<<std::endl;
	    //   dtk::pause();
	    // }
	    start_angle = end_angle;
	    r_cnt[i] += weight;
	    
	    sum+=weight;
	  }
	}
	// std::cout<<"sum: "<<sum<<std::endl;
	// dtk::pause();
      }
    }
  }
}
void Cluster::get_fof_only_galaxies(float r_disrupt,Galaxies& gal){
  std::vector<float> new_x(cp_x,cp_x+cp_size);
  std::vector<float> new_y(cp_y,cp_y+cp_size);
  std::vector<float> new_z(cp_z,cp_z+cp_size);
  std::vector<int>   new_w(cp_size);
  int result_size = cp_size;
  for(int i =0;i<cp_size;++i){
    dis_cp_x.push_back(new_x[i]);
    dis_cp_y.push_back(new_y[i]);
    dis_cp_z.push_back(new_z[i]);
  }
  n2_merger3d(&new_x[0],
	      &new_y[0],
	      &new_z[0],
	      &new_w[0],
	      &result_size,
	      r_fof,
	      NULL);
  for(int i=0;i<result_size;++i){
    //if there are enough partilces in this fof group, 
    //call it a bcg like galaxy. 
    //std::cout<<i<<"/"<<result_size<<": "<<new_w[indx+i]<<std::endl;
    
    if(new_w[i]>=min_disrupt_fof){
      //std::cout<<gal.x[indx+i]<<" "<<gal.y[indx+i]<<" "<<gal.z[indx+i]<<std::endl;
      gal.x.push_back(new_x[i]);
      gal.y.push_back(new_y[i]);
      gal.z.push_back(new_z[i]);
      gal.type.push_back(1);
    }
  }
}
void Cluster::write_out_lb_cores(std::string file_loc){
  std::ofstream file(file_loc.c_str());
  file<<"core_x\tcore_y\tcore_z\tcore_radius\tcore_infall_mass\tcore_infall_time"<<std::endl;
  for(int i =0;i<core_size;++i){
    file<<core_x[i]<<"\t"<<core_y[i]<<"\t"<<core_z[i]<<"\t"<<core_r[i]<<"\t"<<core_m[i]<<"\t"<<stepz.z_from_step(core_step[i])<<std::endl;
  }
  file.close();
}
void Cluster::write_out_lb_prop(std::string file_loc){
  std::ofstream file(file_loc.c_str());
  file<<"halo_x\thalo_y\thalo_z\thalo_m200\thalo_r200\thalo_redshift"<<std::endl;
  file<<x<<"\t"<<y<<"\t"<<z<<"\t"<<sod_mass<<"\t"<<sod_radius<<"\t"<<redshift<<std::endl;
  file.close();
}
void Cluster::move_cores_together(){
  move_together(x,rL,core_x,core_size);
  move_together(y,rL,core_y,core_size);
  move_together(z,rL,core_z,core_size);
}
void ZMR::write_txt(std::string file_loc){
  std::ofstream file(file_loc.c_str());
  write_array(file,"z_bins",z_bins);
  write_array(file,"m_bins",m_bins);
  write_array(file,"r_bins",r_bins);
  write_array(file,"zm_Ngal",zm_Ngal);
  write_array(file,"zm_Ngal_err",zm_Ngal_err);
  write_array(file,"zm_Ngal_var",zm_Ngal_var);
  write_array(file,"zmr_gal_counts",zmr_gal_counts);
  write_array(file,"zmr_gal_density",zmr_gal_density);
  write_array(file,"zmr_gal_density_err",zmr_gal_density_err);
  write_array(file,"zmr_gal_density_var",zmr_gal_density_var);
  write_array(file,"zmr_dgal_dr",zmr_dgal_dr);
  write_array(file,"zmr_dgal_dr_err",zmr_dgal_dr_err);
  write_array(file,"zmr_dgal_dr_var",zmr_dgal_dr_var);
  write_array(file,"zmr_gal_accum",zmr_gal_accum);
  write_array(file,"zmr_gal_accum_err",zmr_gal_accum_err);
  write_array(file,"zmr_gal_accum_var",zmr_gal_accum_var);
  write_array(file,"zm_counts",zm_counts);
  write_array(file,"zmr_counts",zmr_counts);
  file.close();
}
void write_core_params(CoreParam cp, std::string file_loc){
  H5::H5File file(file_loc,H5F_ACC_TRUNC);
  dtk::write_hdf5(file,"/m_infall",cp.m_infall);
  dtk::write_hdf5(file,"/r_disrupt",cp.r_disrupt);
  if(fit_r_fof)
    dtk::write_hdf5(file,"/r_fof",cp.r_fof);
  if(fit_r_merger)
    dtk::write_hdf5(file,"/r_merger",cp.r_merger);
}
int get_num_compact_cores(float m_infall, float r_disrupt, Cores cores){
  int num=0;
  for(int i =0;i<cores.size;++i){
    num += (cores.infall_mass[i]>m_infall) && (cores.radius[i]<r_disrupt);
  }
  return num;
}
float get_locked_r_disrupt(float m_infall,Cores cores){
  float r_min =0;
  float r_max =1.0;
  float box_vol = rL*rL*rL;
  float expected_sum = expected_comov_abundance*box_vol; 
  int r_min_num = get_num_compact_cores(m_infall,r_min,cores);
  int r_max_num = get_num_compact_cores(m_infall,r_max,cores);
  if(expected_sum > r_max_num){
    // std::cout<<"locked_r_disrupt Error: expected_sum ["<<expected_sum<<"] is larger than maximum possible ["<<r_max_num<<"]"<<std::endl;
    // throw;
    return 1.0; //We expect more cores than we actually have. 
  }
  // std::cout<<"we are starting..."<<std::endl;
  do{
    // std::cout<<expected_sum<<std::endl;
    // std::cout<<r_min_num<<" -> "<<r_max_num<<std::endl;
    // std::cout<<r_min<<"->"<<r_max<<std::endl;
    float r_half = (r_min+r_max)/2.0;
    int r_half_num = get_num_compact_cores(m_infall,r_half,cores);
    // std::cout<<"\t"<<r_half<<": "<<r_half_num<<std::endl;
    if(r_half_num < expected_sum){
      r_min = r_half;
      r_min_num = r_half_num;
    }
    else{ // r_half_num < r_max_num
      r_max = r_half;
      r_max_num = r_half_num;
    }
    // std::cout<<"\t"<<r_min_num<<" -> "<<r_max_num<<std::endl;
    // std::cout<<"\t"<<r_min<<"->"<<r_max<<std::endl;
    // dtk::pause();
  }while((r_max-r_min)>0.001);
  return (r_min+r_max)/2.0;
} 
