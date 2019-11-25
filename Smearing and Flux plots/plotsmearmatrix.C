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

TH2D *hist = new TH2D("nu mucc sk2","2D Histrogram of Smearing Matrix",number_rows,0,number_rows-1,number_columns-2,0,number_columns-3);//"distribution name","Name of plot",86,0,86,86,0,86);
  //hist->Fill(smear);
  hist->GetXaxis()->SetTitle("Matrix index(row) Real");
  hist->GetYaxis()->SetTitle("Matrix index(column) Resonctructed");

  for (int i=0;i<number_rows;i++){
     for(int j=2;j<number_columns;j++){
         
         hist->Fill(i-2,j-2,my_array[i][j]);
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
  myCanvas->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/GraphSmearMatrixCorrected.pdf");

}