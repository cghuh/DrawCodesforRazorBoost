//void Combination(vector<int> arr, vector<int> comb, vector<vector<int>> res, int r, int index, int depth)
//{
//  if (r == 0) {
//    vector<int>tmp;
//    cout << "print" << endl;
//    for(int i = 0; i < comb.size(); i++) {
//      tmp.push_back(comb[i]);
//      cout << comb[i] << " ";
//    }
//    cout << endl;
//    res.push_back(comb);
//    cout << "end" << endl;
//    for (int i=res.size()-1;i<res.size();i++) {
//      vector<int> tmp = res.at(i);
//      //cout << tmp.size() << endl;;
//      for (int j=0;j<tmp.size();j++) {
//        cout << tmp.at(j) << " ";
//      }
//      cout << endl;
//    }
//  }
//  else if (depth == arr.size()) return;
//  else {
//    comb[index] = arr[depth];
//    Combination(arr, comb, res, r - 1, index + 1, depth + 1);
//    Combination(arr, comb, res, r, index, depth + 1);
//  }
//}
//

void Combination(vector<int> arr, vector<int> comb, vector<vector<int>> res, int r, int index, int depth) {
  if (r == 0) {
    vector<int>tmp;
    //cout << "print" << endl;
    for(int i = 0; i < comb.size(); i++) {
      tmp.push_back(comb[i]);
      //cout << comb[i] << " ";
    }
    cout << endl;
    res.push_back(comb);
    //cout << "end" << endl;
    for (int i=res.size()-1;i<res.size();i++) {
      vector<int> tmp = res.at(i);
      //cout << tmp.size() << endl;;
      for (int j=0;j<tmp.size();j++) {
        //cout << tmp.at(j) << " ";
      }
      cout << endl;
    }
  }
  else if (depth == arr.size()) return;
  else {
    comb[index] = arr[depth];
    Combination(arr, comb, res, r - 1, index + 1, depth + 1);
    Combination(arr, comb, res, r, index, depth + 1);
  }
}

void test2() {
  int Num=4;
  vector<int> vec;
  for(int i=1;i<Num;++i) vec.push_back(i);

  for(int r=1;r<=Num;r++){
    vector<int> comb(r);
    vector<vector<int>> res;
    Combination(vec, comb, res, r, 0, 0);
/*
    for (int i=1;i<res.size();i++) {
      vector<int> tmp = res.at(i);
      for (int j=0;j<tmp.size();j++) {
       cout << tmp.at(j) << " ";
      }
      cout << endl;
    }
*/
  }
}
