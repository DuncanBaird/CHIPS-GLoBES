#include <fstream>
void readsmeardata2(){

   TCanvas *myCanvas = new TCanvas();
  myCanvas->SetGrid();

ifstream file;
file.open("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/smear.dat");
if (!file.is_open()) {printf("Error\n");}

double my_array [88][88]{};
double value;
for (int i{}; i != 88; ++i) {
    for (int j{}; j != 88; ++j) {
      file >> value;
      if(value==85){
        my_array[i][j] = 0;
      }else{
        my_array[i][j] = value;
      }
    }
}

printf("%f\n",my_array[1][1]);

TH2D *hist = new TH2D("nu mucc sk2","2D Histrogram of Smearing Matrix",88,0,87,88,0,87);//"distribution name","Name of plot",86,0,86,86,0,86);
  //hist->Fill(smear);
  hist->GetXaxis()->SetTitle("Matrix index");
  hist->GetYaxis()->SetTitle("Matrix index");

  for (int i=0;i<88;i++){
     for(int j=0;j<88;j++){
         if(my_array[i][j]!=0.0){printf("%f\n",my_array[i][j]);}
         hist->Fill(i,j,my_array[i][j]);
     }
  }
  hist->Draw("COL");

  gPad->Update();
  myCanvas->Update();
  myCanvas->SaveAs("/home/duncan/Documents/CHIPS Repository/CHIPS-GLoBES/Smearing and Flux plots/GraphSmear.pdf");

}