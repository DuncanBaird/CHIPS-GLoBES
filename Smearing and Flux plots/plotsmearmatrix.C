#include <fstream>
void plotsmearmatrix(){

   TCanvas *myCanvas = new TCanvas();
  myCanvas->SetGrid();

ifstream file;
file.open("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/smear.dat");
if (!file.is_open()) {printf("Error\n");}

const int number_rows = 86;
const int number_columns = 88;
double my_array [number_rows][number_columns]{};
double value;
for (int i{}; i != number_rows; ++i) {
    for (int j{}; j != number_columns; ++j) {
      file >> value;
      if(value==85){
        my_array[i][j] = value;
      }else{
        my_array[i][j] = value;
      }
    }
}

double bins_array[number_rows+1] = {0,0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4};

for (int edge=0;edge<=number_rows;edge++){
  bins_array[edge]=bins_array[edge]+bins_array[edge-1];
  //printf("Bin edges: %f\n",bins_array[edge]);
}


TH2D *hist = new TH2D("nu mucc sk2","2D Histrogram of Smearing Matrix",number_rows,bins_array,number_columns-2,bins_array);//"distribution name","Name of plot",86,0,86,86,0,86);
  //hist->Fill(smear);
  hist->GetXaxis()->SetTitle("Real energy GeV (Matrix Row i)");
  hist->GetYaxis()->SetTitle("Reconstructed energy GeV (Matrix Column j)");


  for (int i=1;i<number_rows;i++){
     for(int j=3;j<number_columns;j++){
         
         hist->Fill(i,j-2,my_array[i][j]);
     }
  }
  hist->Draw("COLZ");

  gPad->Update();
  TPaveStats *st = (TPaveStats*)hist->FindObject("stats");

  //stat box position
  st->SetX2NDC(0.6); 
  st->SetY2NDC(0.15); 
  st->SetX1NDC(0.8); 
  st->SetY1NDC(0.35); 

  st->SetOptStat(111110110);

  myCanvas->Update();
  myCanvas->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/GraphSmearMatrixCorrectedNewAxis.svg");

  for (int bintest = 1;bintest<87;bintest++){
  //printf("bin %d low edge: %f\n",bintest,hist->GetXaxis()->GetBinLowEdge(bintest));
  //printf("value %f\n",my_array[bintest][5]);
  //printf("bin %d centre: %f\n",bintest,hist->GetXaxis()->GetBinCenter(bintest));
  //printf("bin %d up edge: %f\n",bintest,hist->GetXaxis()->GetBinUpEdge(bintest));
  }
}