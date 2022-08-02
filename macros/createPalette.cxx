   void createPalette()
   {
   	TF2 *f1 = new TF2("f2","0.1+(1-(x-2)*(x-2))*(1-(y-2)*(y-2))",1,3,1,3);
   	UInt_t Number = 3;
   	Double_t Red[Number]    = { 1.00, 0.00, 0.00};
   	Double_t Green[Number]  = { 0.00, 1.00, 0.00};
   	Double_t Blue[Number]   = { 1.00, 0.00, 1.00};
   	Double_t Length[Number] = { 0.00, 0.50, 1.00 };
   	Int_t nb=50;
   	TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
   }