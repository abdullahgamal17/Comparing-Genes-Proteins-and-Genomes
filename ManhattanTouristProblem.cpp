#include <vector>
#include <string>
#include <iostream>
#include<algorithm>

using namespace std;

vector <vector <int>> 
ReadDown(unsigned int n , unsigned int m)
{
    vector <int> Column(m+1,0);
    vector <vector <int>> Down (n , Column );
    Column.clear();
    int a;
    for (int i = 0 ; i < n; i++)
    {
        for (int j = 0; j < m+1; j++)
        {
            cin >>a;
            Down[i][j] = a; 
        }
    }
    return Down;
}

vector <vector <int>> 
ReadRight(unsigned int n , unsigned int m)
{
    vector <int> Column(m,0);
    vector <vector <int>> Right (n+1 , Column );
    Column.clear();
    int a;
    for (int i = 0 ; i < n+1; i++)
    {
        for (int j = 0; j < m; j++)
        {
            cin >>a;
            Right[i][j] = a; 
        }
    }
    return Right;
}

int
Manhattan(unsigned int n , unsigned int m , vector <vector <int>> Down,vector <vector <int>> Right)
{
    vector <int> Column(m+1,0);
    vector <vector <int>> s(n+1 , Column );
    int row = n;
    int column = m;
    Column.clear();
    for (int i = 1; i < n+1; i++)
    {
        s[i][0] = s[i-1][0] + Down[i-1][0]; 
    }
    for (int j = 1; j < m+1; j++)
    {
        s[0][j] = s[0][j-1] + Right[0][j-1];
    }
    for (int i = 1; i < n+1; i++)
    {
        for (int j = 1; j < m+1; j++)
        {
            int a, b;
            a = 0;
            b = 0;
            a =s[i-1][j] + Down[i-1][j];
            b = s[i][j-1] + Right[i][j-1];
            int arr [2] = {a, b};
            s[i][j] = *max_element(arr , arr + 2);
        }
        
    }
    return s[n][m];
}


int main()
{
    int n ,m;
    cin >> n >> m ;
    string none;
    vector <vector <int>> Down = ReadDown(n,m);
    cin >> none;
    vector <vector <int>> Right = ReadRight(n,m);
    cout << Manhattan(n,m, Down,Right);
}