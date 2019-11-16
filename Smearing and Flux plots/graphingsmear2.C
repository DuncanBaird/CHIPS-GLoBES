
using std::cout;
using std::endl;
using std::string;

void graphingsmear2(){

  // TFile* f = TFile::Open("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Graphing for 2019-10-24/20190824_000_confFile_20190824_104918_02POMs_NanobeaconB_6200mV_acceptanceTest_nano.root");
  // TTree *t1 =  (TTree*)f->Get("CLBOpt_tree");
  // TBranch *b1 = (TBranch*)t1->GetBranch("TimeStamp_ns");
  TCanvas *myCanvas = new TCanvas();
  myCanvas->SetGrid();

 
  
  int array_size = 86;
  int x_i = 0;




  string fOutFileName("smear.root"); 


  std::ifstream file1("smear.dat");

  double intensity;
  int i;
  int j ;
  double c[88][88];
  if (file1.is_open()) {
    file1.seekg(0);
    while (!file1.eof()) {
      file1 >> i >> j >> intensity;
      c[i][j]=intensity;
      cout<<c[i][j]<<'\n'; 
    }
    file1.close();
  } else cout << "Error, cannot open file 1\n";

}