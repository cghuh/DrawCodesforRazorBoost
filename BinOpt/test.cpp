#include <iostream>
#include <vector>

void generateSubgroups(std::vector<int>& elements, std::vector<int>& subgroup,
                       int subgroup_size, int start) {
    if (subgroup.size() == subgroup_size) {
        // Print the subgroup
        std::cout << "{ ";
        for (int i = 0; i < subgroup.size(); ++i) {
            std::cout << subgroup[i] << " ";
        }
        std::cout << "}" << std::endl;
    } else {
        for (int i = start; i < elements.size(); ++i) {
            subgroup.push_back(elements[i]);
            generateSubgroups(elements, subgroup, subgroup_size, i + 1);
            subgroup.pop_back();
        }
    }
}

int test() {
    std::vector<int> elements = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    int subgroup_size = 10;
    std::vector<int> subgroup;

    for(int i=1;i<=10;i++){
      subgroup.clear();
      subgroup_size = i;
      std::cout << "All possible subgroups of size " << subgroup_size << ":" << std::endl;
      generateSubgroups(elements, subgroup, subgroup_size, 0);
    }

    return 0;
}
/*
const int MAX = 40;
int dp[MAX];

int countGroups(int position, int length, char *num) {
   if (position == length) return 1;
   if (dp[position] != -1) return dp[position];
   //cout << dp[position] << endl;

   dp[position] = 0;

   int res = 0;
   int temp = 0;

   //cout << position << ", " << length;
   for (int i = position; i < length; i++)	{
   cout << "i : " << i << endl;
   cout << (num[i] - '0') << endl;
   temp = countGroups(i + 1, length, num);
   res += temp;
   cout << "count : " << res << ", " << temp << ", " << i << ", " << length << endl;
}

dp[position] = res;

return res;
}

void test(){
  //char num[] = "1379";
  char num[] = "1234";
  //char num[] = "0123456789abcdefghijklmnopqrst";
  int len = strlen(num);
  //cout << len << endl;
  memset(dp, -1, sizeof(dp));
  //countGroups(0, 0, len, num);
  cout << countGroups(0, len, num) << endl;
}
*/
