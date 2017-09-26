void test(){
   int myjets = 3; //example

   int N_comb = 1;
   for(unsigned int i = 0; i < myjets; i++){
     N_comb *= 2;
   }
   double M_min = 9999999999.0;
   int j_count;
   for(int i = 1; i < N_comb-1; i++){
     cout << i << " 번째 루프" << endl;
     //TLorentzVector j_temp1, j_temp2;
     //TLorentzVector j_temp1, j_temp2;
     int itemp = i;
     j_count = N_comb/2;
     int count = 0;
     while(j_count > 0){
       if(itemp/j_count == 1){
         cout << "itemp/j_count is 1" << endl;
         //j_temp1 += myjets[count];
       } else {
         cout << "itemp/j_count is not 1" << endl;
         //j_temp2 += myjets[count];
       }
       itemp -= j_count*(itemp/j_count);
       j_count /= 2;
       count++;
     }
 //    double M_temp = j_temp1.M2()+j_temp2.M2();
 //    // smallest mass
 //    if(M_temp < M_min){
 //      M_min = M_temp;
 //      j1 = j_temp1;
 //      j2 = j_temp2;
 //    }
 //  }
   }
 //  if(j2.Pt() > j1.Pt()){
 //    TLorentzVector temp = j1;
 //    j1 = j2;
 //    j2 = temp;
 //  }
  //test
 } 
