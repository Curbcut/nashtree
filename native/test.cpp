#include "nashtree/inference.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
using nashtree::Forest;
void check(bool b) { if(!b) throw std::runtime_error("test assertion failed"); }
template<class F> void rejects(F f) { bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught); }
int main() {
 try {
  Forest f;f.features=1;f.original_samples=10;f.sample_ids={2,7};f.sample_trees={0,2,3};f.tree_nodes={0,3,4,5};
  f.nodes={{0,1,2,1,0.5,0},{-1,-1,-1,0,0,1},{-1,-1,-1,0,0,2},
    {-1,-1,-1,0,0,0.25},{-1,-1,-1,0,0,-3}};f.metadata="{\"fixture\":true}";
  auto dir=std::filesystem::temp_directory_path()/"nashtree-native-contract";
  check(std::filesystem::create_directory(dir));auto path=(dir/"test.ntf").string();
  f.save(path);auto g=Forest::load(path);check(g.nodes.size()==5 && g.sample_ids==f.sample_ids);
#ifndef _WIN32
  nashtree::sync_path(path);nashtree::sync_path(dir.string());
  rejects([&]{nashtree::sync_path((dir/"missing").string());});
#endif
  const double x[]={std::nextafter(0.5,-INFINITY),0.5,std::nextafter(0.5,INFINITY),NAN,-INFINITY,INFINITY};
  nashtree::MatrixView view{x,6,1,1,6,6};std::vector<double> out(18);
  g.predict_fixed(view,{7,2,7},out.data(),out.size());
  for(int i=0;i<6;++i){check(out[i]==-3);check(out[12+i]==-3);check(out[6+i]==((i==2||i==5)?2.25:1.25));}
  std::vector<double> mapped(6);g.predict_mapped(view,{2,7,2,7,2,7},mapped.data(),mapped.size());
  for(int i=0;i<6;++i) check(mapped[i]==(i%2?-3:out[6+i]));
  rejects([&]{g.predict_fixed(view,{9},mapped.data(),6);});
  rejects([&]{g.predict_fixed({x,6,1,SIZE_MAX,1,6},{2},mapped.data(),6);});
  rejects([&]{g.predict_mapped(view,{2},mapped.data(),6);});
  rejects([&]{g.save(path);});
  auto bad=g;bad.nodes[0].left=0;rejects([&]{bad.validate();});
  bad=g;bad.nodes[0].kind=2;rejects([&]{bad.validate();});
  bad=g;bad.nodes[0].right=1;rejects([&]{bad.validate();});
  bad=g;bad.sample_ids={2,2};rejects([&]{bad.validate();});
  rejects([&]{Forest::load(path,{128,100,100});});
  std::fstream corrupt(path,std::ios::binary|std::ios::in|std::ios::out);corrupt.seekp(-1,std::ios::end);corrupt.put('!');corrupt.close();
  rejects([&]{Forest::load(path);});
  std::filesystem::resize_file(path,32);rejects([&]{Forest::load(path);});
  // Repeated owned handle destruction is also exercised under sanitizers.
  std::filesystem::remove(path);f.save(path);for(int i=0;i<100;++i){auto h=Forest::load(path);h.predict_fixed(view,{2},mapped.data(),6);}
  std::filesystem::remove_all(dir);std::cout<<"native contract tests passed\n";
 } catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}
}
