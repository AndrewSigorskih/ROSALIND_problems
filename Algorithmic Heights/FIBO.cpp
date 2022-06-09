#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

const int NIL = -1;
vector<int> memo = {0, 1};

int fibo_rec_memo(int n)
{
    if (n < 0) return 0;
    if (memo[n] == NIL)
    {
        memo[n] = fibo_rec_memo(n-1) + fibo_rec_memo(n-2);
    }
    return memo[n];
}

int main()
{
    ifstream ist{"rosalind_fibo.txt"};
	  if (!ist) 
	  {
		    cout << "Cannot open input file!\n";
		    exit(1);
	  }
    int n;
    ist >> n;

    for (int i = 2; i <= n; i++)
        memo.push_back(NIL);
    
    cout << fibo_rec_memo(n) << "\n";

    return 0;
}
