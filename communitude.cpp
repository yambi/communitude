#include <algorithm>
#include <bitset>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sstream> 
#include <string>
#include <vector>
#include <queue>

using namespace std;
int n,m;
map<int,set<int>> edges;
set<int> vs;
set<int> restrict;
map<int,bool> fixset;
map<int,bool> removedset;
map<int,string> id2label;
int rep=1;
bool tc = false;//input true community or not
bool fix = false;//input fix sets or not
bool rinit = false;//use a random vertex as an initial state
string iterfile = "iter.txt";
string outfile="result.txt";
bool comflag=true;
bool modflag=true;
bool dsflag=true;
bool oqcflag=true;
bool conductanceflag=true;


void load(char* s){
  n=0;
  m=0;
  int maxn=0;
  ifstream ifs(s);
  int u,v;
  string line;
  //cout<<endl<<restrict.size()<<endl;
  while(getline(ifs,line)){
    istringstream ls(line);
    ls>>u>>v;
    if(restrict.size()>0){
      if(restrict.find(u)==restrict.end())continue;
      if(restrict.find(v)==restrict.end())continue;
    }
    if(ifs.fail())break;
    //cout<<u<<" "<<v<<endl;
    vs.insert(u);
    vs.insert(v);
    maxn=max(maxn,max(u,v));
    if(edges.find(u)==edges.end())edges[u]=set<int>();
    if(edges.find(v)==edges.end())edges[v]=set<int>();
    edges[u].insert(v);
    edges[v].insert(u);
    m++;
  }
  n=vs.size();
  cout<<"n="<<n<<", m="<<m<<endl;
  if(n<=maxn){
    cout<<"========== waraning! ============"<<endl;
    cout<<maxn<<":"<<n<<endl;
    cout<<"skip some label"<<endl<<endl;
  }
}

void load_label(char* s){
  ifstream ifs(s);
  int id;
  string label;
  while(!ifs.eof()){
    ifs>>id>>label;
    if(ifs.fail())break;
    id2label[id]=label;
  }
}

void load_restrict(char* s){
  ifstream ifs(s);
  int id;
  while(!ifs.eof()){
    ifs>>id;
    //cout<<id<<" ";
    if(ifs.fail())break;
    restrict.insert(id);
  }
}


vector<int> randomperm(){
  vector<int> p(n);
  int i=0;
  for(int v:vs){
    p[i]=v;
    i++;
  }
  mt19937 r{random_device{}()};
  shuffle(p.begin(),p.end(),r);
  return p;
}


void show(vector<int> v){
  for(int i:v)cout<<i<<" ";
  cout<<endl;
}
void show(set<int> v){
  for(int i:v)cout<<i<<" ";
  cout<<endl;
}

int intersection_size(set<int>& a, set<int>& b){
  vector<int> intersection;
  set_intersection(a.begin(), a.end(), b.begin(), b.end(), back_inserter(intersection));
  return intersection.size();
}
int union_size(set<int>& a, set<int>& b){
  vector<int> uni;
  set_union(a.begin(), a.end(), b.begin(), b.end(), back_inserter(uni));
  return uni.size();
}

double jaccard(set<int>& a, set<int>& b){
  return double(intersection_size(a,b))/union_size(a,b);
}

double sq(double d){return d*d;}

double com(int cn, int cm, int cdeg){
  double z;
  if(cn==n || cn<=1)z=0;
  else{
    double p = sq(cdeg/2.0/m);
    z = (double(cm)/m - p)/sqrt(p*(1.0-p));
  }
  return z;
}

double mod(int cn, int cm, int cdeg){
  double p = sq(cdeg/2.0/m);
  return  double(cm)/m - p;
}
double densest(int cn, int cm, int cdeg){
  return (cn==0)?0:2.0*cm/cn;
}
double oqc(int cn, int cm, int cdeg){
  return double(cm)-double(cn)*(cn-1.0)/2.0/3.0;
}
double conductance(int cn, int cm, int cdeg){
  if(cn==0 || cn==n)return -1.0E10;
  return -(cdeg-2.0*cm)/min(double(cdeg),2.0*m-cdeg);
}


set<int> peeling(double (*f)(int,int,int)){
  map<int,int> next;for(int v: vs)next[v]=-1;
  map<int,int> prev;for(int v: vs)prev[v]=-1;
  map<int,int> remain;for(int v: vs)remain[v]=1;
  map<int,int> degree;for(int v: vs)degree[v]=edges[v].size();
  vector<int> degnext(n,-1);
  set<int> res;
  int cm = m;
  int cdeg = 0;
  int k=n;//current minimum degree
  int r=0;//# of unfixed element
  for(int v:vs){
    if(!(fixset[v] or removedset[v])){
      if(degnext[degree[v]]>=0){
        int u = degnext[degree[v]];
        prev[u]=v;
        next[v]=u;
      }
      degnext[degree[v]]=v;
      k = min(degree[v],k);
      r++;
      cdeg+=degree[v];
    }
    else if(fixset[v]){
      res.insert(v);
      cdeg+=degree[v];
      remain[v]=0;
    }
    else{
      remain[v]=0;
    }
  }

  // //output degree table
  // for(int i=n-1;i>=0;--i){
  //   cout<<i<<": ";
  //   int v=degnext[i];
  //   while(v>=0){
  //     cout<<v<<" ";
  //     v=next[v];
  //   }
  //   cout<<endl;
  // }

    
  //peel
  vector<int> ordering;
  double maxf = 0;//maximum value
  int maxi=r, maxm, maxdeg;
  for(int i=r;i>0;--i){//# of remaining vertices
    double cf = f(i,cm,cdeg);
    //cout<<i<<": "<<zr<<", "<<cm<<", "<<cdeg<<endl;
    if(cf>maxf){
      maxf=cf;
      maxi=i;
      maxm=cm;
      maxdeg=cdeg;
    }

    while(degnext[k]<0)k++;
    int v = degnext[k];
    ordering.push_back(v);

    //remove v
    cm -= k;
    cdeg-=edges[v].size();
    remain[v]=0;
    if(prev[v]>=0)next[prev[v]]=next[v];
    else degnext[k]=next[v];
    if(next[v]>=0)prev[next[v]]=prev[v];
    //remove edges
    for(int u:edges[v]){
      if(remain[u]){
        //remove u from degree[u] list
        if(prev[u]>=0)next[prev[u]]=next[u];
        else degnext[degree[u]]=next[u];          
        if(next[u]>=0)prev[next[u]]=prev[u];
        //add u to degree[u]-1 list
        degree[u]--;
        int w = degnext[degree[u]];
        if(w>=0)prev[w]=u;
        next[u]=w;
        prev[u]=-1;
        degnext[degree[u]]=u;
      }
    }
    if(k>0)k--;
  }

  for(int i=0;i<maxi;++i){
    res.insert(ordering[r-i-1]);
  }
  //cout<<"f: "<<maxf<<", #nodes="<<maxi<<", #edges="<<maxm<<", deg.="<<maxdeg<<endl;
  return res;
}

set<int> local_search(double (*f)(int,int,int), string s){
  auto const begin = std::chrono::system_clock::now();
  set<int> initial;
  set<int> res;
  if(!fix && !rinit)initial=peeling(f);
  if(fix){
    initial=peeling(f);
  }

  double max_score=-DBL_MAX;
  int maxn,maxm,maxdeg;
  for(int k=0;k<rep;++k){
    set<int> community;
    int cn=0,cm=0,cdeg=0;
    if(rinit){
      random_device rd;
      mt19937 mt(rd());
      uniform_int_distribution<int> r(0,n-1);
      community.insert(r(mt));
    }
    else for(int i: initial)community.insert(i);
    
    for(int i=0;i<n;++i)if(community.find(i)!=community.end()){
        cn++;
        for(int j:edges[i]){
          if(community.find(j)!=community.end())cm++;
          cdeg++;
        }
      }
    cm/=2;
    int loop=-1;
    double score = f(cn,cm,cdeg);
    int repeatnum=0;
    ofstream ofs(iterfile,ios::app);
    //if(tc)ofs<<endl<<endl;
    ofs<<endl<<s<<endl;
    ofs.close();
    while(loop!=0){
      repeatnum++;
      ofstream ofs(iterfile,ios::app);
      ofs<<repeatnum<<": "<<score<<endl;
      ofs.close();
      cout<<repeatnum<<": score="<<score<<", cn="<<cn<<", cm="<<cm<<", exchagne="<<loop<<endl;
      loop=0;
      vector<int> rl = randomperm();
      for(int i:rl)if(!(fixset[i] or removedset[i])){
          int delta_m=0;
          for(int j:edges[i]){
            if(community.find(j)!=community.end())delta_m++;
          }
          //add
          if(community.find(i)==community.end()){
            double newscore = f(cn+1,cm+delta_m,cdeg+edges[i].size());
            if(newscore>score){
              score=newscore;
              cn++;
              cm+=delta_m;
              cdeg+=edges[i].size();
              community.insert(i);
              loop++;
            }
          }
          //remove
          else{
            double newscore = f(cn-1,cm-delta_m,cdeg-edges[i].size());
            if(newscore>score){
              score=newscore;
              cn--;
              cm-=delta_m;
              cdeg-=edges[i].size();
              community.erase(i);
              loop++;
            }
          }
        }
    }
    if(max_score<score){
      max_score=score;
      maxn=cn;
      maxm=cm;
      maxdeg=cdeg;
      res=community;
    }
    cout<<"iterated "<<repeatnum<<" times"<<endl<<endl;
  }
  auto const end = std::chrono::system_clock::now();
  auto t = std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
  cout<<t<<" ms"<<endl;
  char output[256];
  sprintf(output,"time.txt");
  ofstream ofs(output,ios::app);
  ofs <<s<<": "<< n <<" "<< m <<" "<< t << endl;
  ofs.close();

  cout<<"score:"<<max_score<<", cn:"<<res.size()<<endl;
  cout<<"cn:"<<maxn<<", cm:"<<maxm<<", cdeg:"<<maxdeg<<endl;
  if(res.size()<=100)show(res);
  return res;
}

double choose2(double i){
  return i*(i-1)/2;
}
double choose3(double i){
  return i*(i-1)*(i-2)/6;
}


void output(set<int> s, string output){
  ofstream ofs(output,ios::app);
  for(int i: s){
    if(id2label.find(i)==id2label.end())ofs<<i<<" ";
    else ofs<<id2label[i]<<" ";
  }
  ofs<<endl;
  ofs.close();
}



void usage(char* s){
  cout<<"usage: "<<s<<" network.dat [options]"<<endl;
  cout<<"  -t TRUE_COMMUNITY_FILE"<<endl;
  cout<<"  -f FIXED_ELEMENTS_FILE"<<endl;
  cout<<"  -o RESULT_FILE"<<endl;
  cout<<"  -r REPEAT_NUMBER"<<endl;
}


int main(int argc, char *argv[]){
  if(argc<2){
    usage(argv[0]);
    return 1;
  }

  set<int> trueset;
  //fixset=vector<bool>(n,false);
  for(int i=2;i<argc;++i){
    string s(argv[i]);
    if(s == "-f"){
      i++;
      fix=true;
      ifstream ifs(argv[i]);
      while(!ifs.eof()){
        int k;ifs>>k;
        fixset[k]=true;
      }
    }
    else if(s == "-rinit"){
      rinit=true;
    }
    else if(s == "-t"){
      i++;
      tc=true;
      ifstream ifs(argv[i]);
      while(!ifs.eof()){
        int k;ifs>>k;
        trueset.insert(k);
      }
    }
    else if(s == "-i"){
      i++;
      load_restrict(argv[i]);
    }
    else if(s == "-o"){
      i++;
      outfile=string(argv[i]);
    }
    else if(s == "-r"){
      i++;
      rep=atoi(argv[i]);
    }
    else if(s == "-l"){
      i++;
      load_label(argv[i]);
    }
    else if(s == "-all-off"){
      comflag=false;
      modflag=false;
      dsflag=false;
      oqcflag=false;
      conductanceflag=false;
    }
    else if(s == "-com-off"){
      comflag=false;
    }
    else if(s == "-mod-off"){
      modflag=false;
    }
    else if(s == "-ds-off"){
      dsflag=false;
    }
    else if(s == "-oqc-off"){
      oqcflag=false;
    }
    else if(s == "-conductance-off"){
      conductanceflag=false;
    }
    else if(s == "-com"){
      comflag=true;
    }
    else if(s == "-mod"){
      modflag=true;
    }
    else if(s == "-ds"){
      dsflag=true;
    }
    else if(s == "-oqc"){
      oqcflag=true;
    }
    else if(s == "-conductance"){
      conductanceflag=true;
    }
    else{
      usage(argv[0]);
      return 1;
    }
  }
  load(argv[1]);


  if(fix){
    cout<<"fixed: ";
    for(int i=0;i<n;++i){
      if(fixset[i])cout<<i<<" ";
    }
    cout<<endl;
  }
  set<int> res;    
  ofstream ofs(outfile,ios::app);

  if(comflag){
    cout<<endl<<"communitude"<<endl;
    remove("data/com.group");
    res=local_search(com,"communitude");
    output(res,"data/com.group");
    if(tc)ofs << jaccard(res,trueset) << " ";
    for(int v:res)removedset[v]=true;
    cout<<endl;
    removedset=map<int,bool>();
  }
  else{
    if(argc>=3)ofs << "* ";
  }

  if(modflag){
    cout<<endl<<"modularity"<<endl;
    remove("data/mod.group");
    res=local_search(mod,"modularity");
    output(res,"data/mod.group");
    if(tc)ofs << jaccard(res,trueset) << " ";
    for(int v:res)removedset[v]=true;
    cout<<endl;
    removedset=map<int,bool>();
  }
  else{
    if(tc)ofs << "* ";
  }

  if(dsflag){
    cout<<endl<<"densest subgraph"<<endl;
    remove("data/ds.group");
    res=local_search(densest,"densest subgraph");
    output(res,"data/ds.group");
    if(tc)ofs << jaccard(res,trueset) << " ";
    for(int v:res)removedset[v]=true;
    cout<<endl;
    removedset=map<int,bool>();
  }
  else{
    if(argc>=3)ofs << "* ";
  }

  if(oqcflag){
    cout<<endl<<"oqc"<<endl;
    remove("data/oqc.group");
    res=local_search(oqc,"OQC");
    output(res,"data/oqc.group");
    if(tc)ofs << jaccard(res,trueset) << " ";
    for(int v:res)removedset[v]=true;
    removedset=map<int,bool>();
  }
  else{
    if(tc)ofs << "* ";
  }

  if(conductanceflag){
    cout<<endl<<"conductance"<<endl;
    remove("data/conductance.group");
    res=local_search(conductance,"conductance");
    output(res,"data/conductance.group");
    if(tc)ofs << jaccard(res,trueset) << " ";
    for(int v:res)removedset[v]=true;
    removedset=map<int,bool>();
  }
  else{
    if(tc)ofs << "* ";
  }
  
  ofs << endl;
  ofs.close();
}
