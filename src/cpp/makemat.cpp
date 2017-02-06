#include "newick.hpp"
#include "TaxonSet.hpp"
#include <algorithm>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cctype>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <map>
#include <cassert>

namespace po = boost::program_options;
using namespace std;


int main(int argc, char** argv) {

  bool bygene=false;

  po::options_description desc("ASTRID v2.0");

  //desc.add_options()
  //  ("input", po::value<string>(), "gene trees in newick format");

  desc.add_options()
    ("matrix", po::value<string>(), "output for matrix");

  desc.add_options()
    ("taxlist", po::value<string>(), "output for taxon mapping file");  

  desc.add_options()
    ("n_missing", po::value<string>(), "number of missing elements");  
  
  desc.add_options()
    ("taxcutoff", po::value<int>()->default_value(0), "minimum taxa in a tree");  
  
  desc.add_options()
    ("nanplaceholder", po::value<string>()->default_value("--"), "representation for nan in matrix");  

  desc.add_options()
    ("inputmap", po::value<string>()->default_value(""), "input of species/gene_copies mapping file");

  desc.add_options()
    ("bygene", po::bool_switch(&bygene),"AGIDs are calculated by gene and then averaged by the number of genes");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    
  
  bool multi=vm["inputmap"].as<string>().compare("")!=0;
  vector<string> trees;
  
  //Tree input 
  //ifstream infile(vm["input"].as<string>());
  string line;
  while (getline(cin, line)) {
    trees.push_back(line);
  }

  //Taxa parsing
  unordered_set<string> taxa;
  for (string& tree : trees) {
    int ntaxa = newick_to_ts(tree, taxa);
    if (ntaxa < vm["taxcutoff"].as<int>()) {
      tree.clear();
    }
  }
  TaxonSet ts(taxa.size()); 
  for (string t : taxa) {
    ts.add(t);
  }
  ts.freeze();

  //Variables that will be used inside the conditional statement 
  dm_type* pdist_mat;
  dm_type* pmask_mat;
  int mat_size;
  TaxonSet* finalset;
  TaxonSet* ss;

  if(multi) {
    //Species parsing
    unordered_set<string> species;
    map<string,string> gcopy_sstring_map;
    int ts_to_ss[taxa.size()];

    typedef boost::tokenizer<boost::char_separator<char> > tokenizer; 
    boost::char_separator<char> sep(" ");
    
    ifstream infile(vm["inputmap"].as<string>());
    while (getline(infile,line)) {
      tokenizer tokens(line, sep);
      int n_tok=0;
      vector<string> temp;
      temp.assign(tokens.begin(),tokens.end());
      assert(temp.size()==2);
      gcopy_sstring_map[temp[0]]=temp[1];
      species.insert(temp[1]);
    }

    ss= new TaxonSet(species.size());
    for (string s : species) {
      ss->add(s);
    }
    ss->freeze();
    finalset=ss;

    //Taxon id to species id array for quick conversion
    for (const auto& t : taxa) {
      ts_to_ss[ts[t]]=(*ss)[gcopy_sstring_map[t]];
      assert(gcopy_sstring_map.count(t)!=0);//All gene copies have to be contained in the map
    }
   
    //Matrix initialization 
    dm_type dist_mat_temp(boost::extents[taxa.size()][taxa.size()]);
    dm_type mask_mat_temp(boost::extents[taxa.size()][taxa.size()]); 
    pdist_mat= new dm_type(boost::extents[species.size()][species.size()]);
    pmask_mat= new dm_type(boost::extents[species.size()][species.size()]);
    mat_size=species.size();

    int pi;
    int pj;

    if (bygene)
    {
      dm_type dist_mat_bylocus(boost::extents[species.size()][species.size()]);
      dm_type mask_mat_bylocus(boost::extents[species.size()][species.size()]);
      for (string& tree : trees)  {
        if (tree.size()) {
           //Reset matrices
           for (int i=0; i<taxa.size()-1; ++i) for (int j=i+1; j<taxa.size();++j) { //Only upper triangle
             dist_mat_temp[i][j]=0;
             mask_mat_temp[i][j]=0;
           } 
           for (int i=0; i<mat_size; ++i) for (int j=0; j<mat_size;++j) {
             dist_mat_bylocus[i][j]=0;
             mask_mat_bylocus[i][j]=0;
           }
  
           //Calculate GIDs
           newick_to_dm(tree, ts, dist_mat_temp, mask_mat_temp);
 
           //Summarize per species
           for (int i=0; i<taxa.size()-1; ++i) for (int j=i+1; j<taxa.size();++j) { //Only upper triangle
             pi=ts_to_ss[i];
             pj=ts_to_ss[j];
             dist_mat_bylocus[pi][pj]+=dist_mat_temp[i][j];
             mask_mat_bylocus[pi][pj]+=mask_mat_temp[i][j];
             dist_mat_bylocus[pj][pi]=dist_mat_bylocus[pi][pj];
             mask_mat_bylocus[pj][pi]=mask_mat_bylocus[pi][pj];
           }

          //Add mean to the matrix of means, with 1 or 0 to the mask indicating the number of means added per comparison
          for (int i=0; i<mat_size-1;++i) for (int j=i+1; j<mat_size;++j) { //Upper triangle
            if(mask_mat_bylocus[i][j]!=0){
              (*pdist_mat)[i][j]+=dist_mat_bylocus[i][j]/mask_mat_bylocus[i][j];
              (*pmask_mat)[i][j]+=1;
            }
          }
        }
      }
    
    //Filling up lower triangle
    for (int i=0; i<mat_size-1;++i) for (int j=i+1; j<mat_size;++j) { //Upper triangle
      (*pdist_mat)[j][i]=(*pdist_mat)[i][j];
      (*pmask_mat)[j][i]=(*pmask_mat)[i][j];
    }

  } else {
      for (string& tree : trees)  {
        if (tree.size()) {
           //Reset matrices
           for (int i=0; i<taxa.size()-1; ++i) for (int j=i+1; j<taxa.size();++j) { //Only upper triangle
             dist_mat_temp[i][j]=0;
             mask_mat_temp[i][j]=0;
           }
 
           //Calculate GIDs
           newick_to_dm(tree, ts, dist_mat_temp, mask_mat_temp);
           
           //Summarize per species
           for (int i=0; i<taxa.size()-1; ++i) for (int j=i+1; j<taxa.size();++j) { //Only upper triangle
             pi=ts_to_ss[i];
             pj=ts_to_ss[j];
             (*pdist_mat)[pi][pj]+=dist_mat_temp[i][j];
             (*pmask_mat)[pi][pj]+=mask_mat_temp[i][j];
             (*pdist_mat)[pj][pi]=(*pdist_mat)[pi][pj];//Transposed elements
             (*pmask_mat)[pj][pi]=(*pmask_mat)[pi][pj];
           }
        }
      }
    }

  } else {
    //Matrix initialization
    pdist_mat= new dm_type(boost::extents[taxa.size()][taxa.size()]);
    pmask_mat= new dm_type(boost::extents[taxa.size()][taxa.size()]);
    mat_size=taxa.size();
    finalset=&ts;

    for (string& tree : trees)  {
      if (tree.size())
         newick_to_dm(tree, ts, *pdist_mat, *pmask_mat);
    }

  }

  //Using references to make the code more readable
  dm_type& dist_mat=*pdist_mat;  
  dm_type& mask_mat=*pmask_mat;
  

  ofstream matfile(vm["matrix"].as<string>());
  ofstream taxlist(vm["taxlist"].as<string>());
  ofstream missingfile(vm["n_missing"].as<string>());    

  int n_missing = 0;
  
  matfile << mat_size << endl;
  for (int i = 0; i < mat_size; i++) {
    dist_mat[i][i] = 1;
    mask_mat[i][i] = 1;    
    matfile << i << " ";
    for (int j = 0; j < mat_size; j++)  {
      dist_mat[i][j] /= mask_mat[i][j];
      if (mask_mat[i][j])
	matfile << dist_mat[i][j]-1  << " ";
      else {
	matfile << vm["nanplaceholder"].as<string>()  << " ";
	n_missing ++;
      }
    }
    matfile << endl;
    taxlist << (*finalset)[i] << endl;
  }

  missingfile << n_missing << endl;
  
}
