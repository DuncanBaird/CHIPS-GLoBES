{
   // example of macro to read data from an ascii file and
   // create a root file with an histogram and a TTree
   gROOT->Reset();

   // the structure to hold the variables for the branch

   struct staff_t {
      Int_t i;
      Int_t j;
      Int_t value;
   };
   staff_t staff;
   // continued...
   // open the ASCII file
   FILE *fp = fopen("smear.dat","r");
   char line[81];
   // create a new ROOT file
   TFile *f = new TFile("staff.root","RECREATE");
   // create a TTree
   TTree *tree = new TTree("T","staff data from ascii file");
   // create one branch with all information from the stucture
   tree->Branch("staff",&staff);
   // fill the tree from the values in ASCII file
   while (fgets(&line,80,fp)) {
      sscanf(&line[0],"%d%d%d%d",&staff.i,&staff.j,
          &staff.value);
      tree->Fill();
   }
   // check what the tree looks like
   tree->Print();
   fclose(fp);
   f->Write();
}