#ifndef CATALOG_HPP
#define CATALOG_HPP
#include <iostream>

#include "dtk/all.hpp"


class Catalog{
  
  class Data{
    size_t size;
    void* get_data();
  };
  template<class T>
  class DataSpec{
    
  };
  std::vector<std::pair<std::string,std::string> > var_int32;
  std::vector<std::pair<std::string,std::string> > var_int64;
  std::vector<std::pair<std::string,std::string> > var_float32;
  std::vector<std::pair<std::string,std::string> > var_float64;

  std::vector<std::vector<int> >     data_int32;
  std::vector<std::vector<int64_t> > data_int64;
  std::vector<std::vector<float> >   data_float32;
  std::vector<std::vector<double> >  data_float64;
  
  std::string file_name;

public:
  Catalog();
  Catalog(std::string file_name);
  void read_gio();
  void merge(Catalog c1,Catalog c2);
  template <typename T>
  void add_variable_name(std::string var_name){
    //add the variable to catalog as it's own name;
    add_variable_name<T>(var_name,var_name); 
  }
  template <typename T>
  void add_variable_name(std::string var_name,std::string as_name);
  
};

template <> void Catalog::add_variable_name<int>    (std::string var_name,std::string as_name);
template <> void Catalog::add_variable_name<int64_t>(std::string var_name,std::string as_name);
template <> void Catalog::add_variable_name<float>  (std::string var_name,std::string as_name);
template <> void Catalog::add_variable_name<double> (std::string var_name,std::string as_name);


#endif// CATALOG_HPP
