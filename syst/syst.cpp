void syst(){

  int nsyst = 11;
  int nSFsyst = 17;
  nsyst += nSFsyst;
  nsyst--;
  double num1, num2, num3;
  num1 = 0.0; num2 = 1.0; num3 = -1.0;

  cout << fixed;
  cout.precision(1);
  for(int j=0;j<=nsyst;j++){   //Y-axis
    for(int i=0;i<=nsyst;i++){ //X-axis
      if(i==j) cout << std::setw(5) << num2;
      else     cout << std::setw(5) << num1; 
      if(i==nsyst) cout << " 3 0" << endl;
    }
    for(int i=0;i<=nsyst;i++){ //X-axis
      if(i==j) cout << std::setw(5) << num3;
      else     cout << std::setw(5) << num1; 
      if(i==nsyst) cout << " 3 0" << endl;
    }
  }
}
