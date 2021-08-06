#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

int
TheChangeProblem(int Money , vector <int> Coins)
{
    vector <int> MinNumCoins;
    MinNumCoins.push_back(0);
    for (int m = 1; m <= Money ; m++)
    {
        MinNumCoins.push_back(1000000);
        for (int i = 0; i < Coins.size(); i++)
        {
            if (m >= Coins[i])
            {
                if (MinNumCoins[m-Coins[i]] + 1 < MinNumCoins[m])
                {
                    MinNumCoins[m] = MinNumCoins[m - Coins[i]] + 1;
                }
                
            }
            
        }
        
    }
    return MinNumCoins[Money];
}

int main ()
{
    int Money = 25;
    vector <int> Coins = {2,3};
    cout << TheChangeProblem(Money , Coins);
}
