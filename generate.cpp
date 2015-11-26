#include <iostream>
#include <fstream>
#include <ctime>
#include <random>

using namespace std;

void usage(string s){
  cout<<"usage: "<<s<<endl;
  cout<<"\t-n NUMBER_OF_VERTICES\t(default 100)"<<endl;
  cout<<"\t-c COMMUNITY_SIZE\t(default 10)}"<<endl;
  cout<<"\t-in CONNECTED_PROBABILITY_OF_INNER_PAIR\t(default 0.5)"<<endl;
  cout<<"\t-out CONNECTED_PROBABILITY_OF_OUTER_PAIR\t(default 0.1)"<<endl;
  cout<<"\t-single\t{single community mode (default)}"<<endl;
  cout<<"\t-multi\t{multi community mode}"<<endl;
  cout<<"\t-net NETWORK_FILE\t(default data/network.dat)"<<endl;
  cout<<"\t-com COMMUNITY_FILE\t(default data/community.dat)"<<endl;
  cout<<"\t-fix FIX_FILE\t(default data/fix.dat)"<<endl;  
}
int main(int argc, char* argv[]){
  int n=100;
  int m=0;
  int c=10;//community size
  int fixnum=5;
  double pin=0.5;
  double pout=0.1;
  bool single = true;
  string netfile="data/network.dat";
  string comfile="data/community.dat";
  string fixfile="data/fix.dat";
  
  for(int i=1;i<argc;++i){
    string s(argv[i]);
    if(s == "-n"){
      i++;
      n=atoi(argv[i]);
    }
    else if(s == "-c"){
      i++;
      c=atoi(argv[i]);
    }
    else if(s == "-in"){
      i++;
      pin=atof(argv[i]);
    }
    else if(s == "-out"){
      i++;
      pout=atof(argv[i]);
    }
    else if(s == "-fixnum"){
      i++;
      fixnum=atoi(argv[i]);
    }
    else if(s == "-single"){
      single=true;
    }
    else if(s == "-multi"){
      single=false;
    }
    else if(s == "-net"){
      i++;
      netfile=atoi(argv[i]);
    }
    else if(s == "-com"){
      i++;
      comfile=atoi(argv[i]);
    }
    else if(s == "-fix"){
      i++;
      fixfile=atoi(argv[i]);
    }
    else{
      usage(argv[0]);
      return 1;
    }
  }
  mt19937 r{random_device{}()};
  uniform_real_distribution<double> d( 0.0, 1.0 ) ;
  auto in = [&]{return d(r)<pin;};
  auto out = [&]{return d(r)<pout;};

  ofstream netfs(netfile);
  for(int i=0;i<n;++i){
    for(int j=i+1;j<n;++j){
      if((single and i<c and j<c and in()) or
         (single and (i>=c or j>=c) and out()) or
         (!single and (int(i/c)==int(j/c)) and in()) or
         (!single and (int(i/c)!=int(j/c)) and out())
         ){
        netfs<<i<<" "<<j<<endl;
        m++;
      }
    }
  }
  netfs.close();
  cout<<"network data: "<<comfile<<endl;

  ofstream comfs(comfile);
  for(int i=0;i<c;++i)comfs<<i<<endl;
  comfs.close();
  cout<<"community data: "<<comfile<<endl;

  if(!single){
    ofstream fixfs(fixfile);
    for(int i=0;i<fixnum;++i)fixfs<<i<<endl;
    fixfs.close();
    cout<<"fix data: "<<fixfile<<endl;
  }


  cout<<"n="<<n<<", m="<<m<<", c="<<c<<", pin="<<pin<<", pout="<<pout<<endl;
  return 0;
}
