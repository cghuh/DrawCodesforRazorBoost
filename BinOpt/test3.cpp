//https://rosettacode.org/wiki/Combinations#C.2B.2B
std::vector<int> comb(int N, int K) {
  std::string bitmask(K, 1);
  bitmask.resize(N, 0);
  std::vector<int> store;
  do {
    for (int i = 0; i < N; ++i) {
      if (bitmask[i]) store.push_back(i+1);
    }
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
  return store;
}

void test3() {
  const int N = 15;
  std::vector<int> store;
  for(int i=1;i<=N;i++){
    store = comb(N, i);
    for(int j=0;j<store.size();j++) {
      cout << " " << store.at(j);
      if((j+1)%i == 0) cout << endl;
    }
  }
}
