#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

string parse(vector <vector<int>> a)
{
    string aa = "{ {";
    for (int i = 0; i < a.size(); i++)
    {
        for(auto j:a[i])
        {

            aa += to_string(j)  ;
            aa +=  ","; 
        }
        aa += "} , {";
    }
    return aa;
}
string NumberToSymbol(int a)
{
    string value = "n";

    switch(a)
    {
        case 0 :
        value = "A";
        break;

         case 1 :
        value = "B";
        break;


        default :
        value = "r";
    }
    return value ;
}

string NumberToPattern(int index , int k)
{

    if (k == 1 )
    {
        return NumberToSymbol(index) ;
    }
    string PrefixPattern;
    unsigned int prefixIndex = index/2;
    int r = index%2;
    string symbol = NumberToSymbol(r);
    
    PrefixPattern = NumberToPattern(prefixIndex,k-1);
    return PrefixPattern + symbol;
}

unsigned short 
SymbolToNumber(char symbol)
{
    int value;
    switch (symbol)
    {
        case 'A':
        value = 0;
        break;

         case 'B':
        value = 1;
        break;
    }
    return value;
} 

unsigned int
PatternToNumber(string Pattern)
{
    if (Pattern.empty())
    {
        return 0 ;
    }
    char symbol = Pattern.back();
    string prefix = Pattern.substr(0,Pattern.length()-1);
    return 2*PatternToNumber(prefix) + SymbolToNumber(symbol);
}

int 
Test(int NumberOfAminoAcids , int MassOfAminoAcid)
{
    string arr[2] = {"A","B"};
    int arrr[2] = {2,3};
    int n = MassOfAminoAcid/2;
    vector <string> AllProtein(pow(2,n)*n,"");
    int c = 0;
    for(int j = 1 ; j <= n ;j++)
    {
        for (int i = 0; i < pow(2,j); i++)
        {
            AllProtein[c] = NumberToPattern(i,j);
            c++;
        }
    }
    vector <int> Masses(pow(2,n)*n,0);
    for (int i = 0; i < AllProtein.size(); i++)
    {
        int ProteinMass = 0;
        for (int j = 0; j < AllProtein[i].length(); j++)
        {
            if (AllProtein[i][j] == 'A')
            {
                ProteinMass += 2;
            }
            if (AllProtein[i][j] == 'B')
            {
                ProteinMass += 3;
            }
            
        }
        Masses[i] = ProteinMass;   
    }
    int count = 0;
    for(auto v: Masses)
    {
        if ( v == MassOfAminoAcid)
        {
            count++;
        }
    }
    return count;
}

int main()
{
    vector <int> a(20,0);
    vector <vector<int>> aa(20,a);
    for(int i = 0 ; i < aa.size() ; i++ )
    {
        vector <int> a(20,0);
        for (int j = 0; j < a.size(); j++)
        {
            int x;
            cin >> x;
            a[j] = x;
        }
        aa[i] = a;
    }
    cout << parse(aa);
}


