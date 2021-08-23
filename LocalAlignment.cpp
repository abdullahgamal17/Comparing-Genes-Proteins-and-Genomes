#include <iostream>
#include <vector>
#include <string> 
#include <cmath> 
#include <iomanip>
#include <algorithm>

using namespace std;


void BackTracking(vector <vector<char>> Matrix,string firstAminoAcid,string secondAminoAcid,pair<int,int> coordinates)
{
   int rows = Matrix.size();
   int columns = Matrix[0].size();
   vector<char> c,d;
   int a = coordinates.first-1;
   int b = coordinates.second-1;
    // 0 is diagonal , 1 Horizontal , 2 Vertical;
   for (int i = coordinates.first; i >= 0 ; )
   {
      for (int j = coordinates.second; j >=0 ;)
      {
         char S = Matrix[i][j];
         if(Matrix[i][j] == 'C')
         {
            c.push_back(firstAminoAcid[a]);
            d.push_back(secondAminoAcid[b]);
            a--;
            b--;
            i--;
            j--;
            continue;
         }
         else if(Matrix[i][j] == 'D')
         {
            i--;
            c.push_back(firstAminoAcid[a]);
            d.push_back('-');
            a--;
            continue;
         }
         else if(Matrix[i][j] == 'R')
         {
            j--;
            c.push_back('-');
            d.push_back(secondAminoAcid[b]);
            b--;
            continue;
         }
         if (i == 0 && j == 0)
         {
            break;
         }
         if (Matrix[i][j] == 'S')
         {
             break;
         }
         
         
      }
      break;
   }
   

for (int i = c.size()-1; i >=0 ; i--)
{
   cout << c[i];
}
cout << endl;
for (int i = d.size()-1; i >=0 ; i--)
{
   cout << d[i];
}
cout << endl;
   

   //return DiorGap; 
}

vector <vector<int>> inputMatrix(int rows,int columns)
{
    vector <vector<int>> PAM250(20);
    for (int i = 0; i < 20; i++)
    {
        vector <int> row(20);
        for (int j = 0; j < 20; j++)
        {
            int x;
            cin >> x;
            row[j] = x; 
        }
        PAM250[i] = row;
        
    }
    return PAM250;
}

int GetPos(char AminoAcid)//Finds The Location of Amino Acid in BLOSUM62 Matrix
{
   switch(AminoAcid)
   {
      case 'A' :
      return 0;
      
      case 'C' :
      return 1;

      case 'D' :
      return 2;

      case 'E' :
      return 3;

      case 'F' :
      return 4;

      case 'G' :
      return 5;
      
      case 'H' :
      return 6;

      case 'I' :
      return 7;

      case 'K' :
      return 8;

      case 'L' :
      return 9;

      case 'M' :
      return 10;
      
      case 'N' :
      return 11;

      case 'P' :
      return 12;

      case 'Q' :
      return 13;

      case 'R' :
      return 14;

      case 'S' :
      return 15;
      
      case 'T' :
      return 16;

      case 'V' :
      return 17;

      case 'W' :
      return 18;

      case 'Y' :
      return 19;
   }
}

template <class T>
void OutputMatrix(vector <vector<T>> Mat)
{
    int rows = Mat.size();
    int columns = Mat[0].size();
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            cout << setw(5) << Mat[i][j] << " ";
        }
        cout << endl;   
    }
    
}
/*
int maximum(vector<int> arr)
{
    int max = -50000;
    for (int i = 0; i < arr.size(); i++)
    {
        if (max < arr[i])
        {
            max = arr[i];
        }
        
    }
    return max;
}
*/
class Graph
{
private:
    string aminoacid1;
    string aminoacid2;
    int aminoacid1length;
    int aminoacid2length;
    int indel;
    int maximumscore;
    pair<int,int>pos;
    vector<vector<int>> PAM250;
    vector <vector<int>> GraphTable;
    vector <vector<char>> BackTrackTable;
    
public:
    Graph(string FirstAminoAcid , string SecondAminoAcid,vector<vector<int>> X)
    {
        aminoacid1 = FirstAminoAcid;
        aminoacid2 = SecondAminoAcid;
        PAM250 = X;
        aminoacid1length = FirstAminoAcid.length();
        aminoacid2length = SecondAminoAcid.length();
        indel = 5;
        maximumscore = -1;
        vector<int> init (aminoacid2length+1,0);
        vector<char> ini(aminoacid2length+1,'X');

        for (int i = 0; i <= aminoacid1length; i++)
        {
            GraphTable.push_back(init);
            BackTrackTable.push_back(ini);
        }
        
    }
    void FillGraph()
    {
        for (int i = 1; i <= aminoacid2length; i++)
        {
            GraphTable[0][i] = GraphTable[0][i-1] ;
            BackTrackTable[0][i] = 'S';  
        }
        for (int i = 1; i <= aminoacid1length; i++)
        {
            GraphTable[i][0] = GraphTable[i-1][0] ;
            BackTrackTable[i][0] = 'S';  
        }
        for (int i = 1; i <= aminoacid1length; i++)
        {
            for (int j = 1; j <= aminoacid2length; j++)
            {
                int UpNode = GraphTable[i-1][j];
                int LeftNode = GraphTable[i][j-1];
                int CornerNode = GraphTable[i-1][j-1];
                int freetaxiride = 0;
                int Vertical = UpNode - indel;
                int Horizontal = LeftNode - indel;
                char aminoacid1char = aminoacid1[i-1];
                char aminoacid2char = aminoacid2[j-1];
                int pos1 = GetPos(aminoacid1char);
                int pos2 = GetPos(aminoacid2char);
                int blval = PAM250[pos1][pos2];
                int Diagonal = CornerNode + blval;
                vector <int> Nodes = {Horizontal,Vertical ,Diagonal,freetaxiride};
                int maxx = *max_element(Nodes.begin(),Nodes.end());
                GraphTable[i][j] = maxx;
                if (GraphTable[i][j]> maximumscore)
                {
                    maximumscore = GraphTable[i][j];
                    pos.first = i;
                    pos.second = j;
                }
                
                if (maxx == Vertical)
                {
                    BackTrackTable[i][j] = 'D';
                    continue;
                }
                else if(maxx == Horizontal )
                {
                    BackTrackTable[i][j] = 'R';
                    continue;
                }
                else if(maxx == Diagonal)
                {
                    BackTrackTable[i][j] = 'C';
                    continue;
                }
                else if(maxx == freetaxiride)
                {
                    BackTrackTable[i][j] = 'S';
                    continue;
                }
            }
            
        }
        
           
    }
    vector <vector<int>> GetGraphTable()
    {
        return GraphTable ; 
    }
    vector <vector<char>> GetBackTrackTable()
    {
        return BackTrackTable ; 
    }
    pair<int,int> getpos()
    {
        return pos;
    }
    int getmaxscore()
    {
        return maximumscore;
    }

};

int main()
{
    string aminoacid1 = "MDQRCMLEMFSSLLQGCHNGTMPFLRDKNWNDGHLTMLPLNNKGCGRICVVTDVSYLCPVGCIDWICHYYVQREMMDDHATMFVYLRGMSCKSHRQTAVPGWGVTMDQSVKPANWDCCPCGIEAISVIKLNCERHGMWYQRTWTVAHENFYILFYAIGMSRMCQFIKNGDLQQVAMNMQAKYTSGMKEGPVDEYTNIYIFPQGFLDHDYDWMHFEQRDSPFGHGQVRDYVIKIMTIRNYNTDQLGCKARAHMFDFECLCMNWGHTHYCINWARVHSKKVLWLTKHRWFQTESNQRAKENGSNVDRNHLVYYAQRGSHWHCYGTVPVAIRMICNWEWPKTWFVIKFMQNVRTLWGRAMSGQNDDMRNQPESIEPLNSPFDFAPHFITCTSKQPITEWMWHENVVMIHHDFWMAAELNQFGGPHTVGGNCKAQHWEPQKSVCMIIAKQVDDHDIRRRTMKDQYGFCTFKPAVNDIVQLQCDKPCTDSAMYPAPPWAYSVLLFQVRRGNAYSTNGRETCFMVEDKTEWDTQCMNNCKAHHMCHPNMAHHCCIETRWSQHFTMGTAMWTTGQDAMIVWYSIWGDRWDRWRDDHSEVMRDCCIYCKEKVAGLICHYVTRNGAHDNCRPNRCNQCGMMQDAGDKSIYELFRKQGMGWRAQSNIDARVTKMRPYKREVFGVLWKKLTWPLNKEAWVDDARVYVGRPSCIWCPSYDRQTQLWFEVQDGDSRDAAKNFSQGTNRFWQLTTWHVPVPFLSMGFFTEYMCEWFGSYPNGYTNIFMVDANSPQVPTGRSIPFKLFCFCLECQQDQLYLPCIMCKMEMDKKAEPLTQVFFDHLGKEQYSHLWRRTASLNRLWEIYYCMWRAPLDMHEDMECPCPIAQFARIRTSHDMCRENEFLQPHCGCKHI";
    string aminoacid2 = "WQDACQKLAEVTEEITWWYLALGVRNFYERSQYQWTPWRHIFMYETEQWQMREESNWKEQYLSMWMNLFAPSQNTCGHNQLKVTYHSHQMFGCQSDIGSKYPDITHYRQPEGRKSPKALKKGVAPTGLADGHKMGFPIANVWRHGTQKIVQRETTAHPQAARYLKAAFTAHMILEFMKFMVMEALKACKSMMVIHFWLIFYGQGELHKSMYNMLQCTMAIPDNSFYKDGGMIHLPETQDRPLFWDSEDYRRPIISCQGLYHTQDFVLPKRIRALKTDFECLCMNWGHTHYCINWGATATPRVHSKKELFVAVTCGVLSMDSRFLSLESNQRAKENVSNVDRNHLVYYAQCIKNYFSSHWHCYGTDAMQAETPVAICNWEWPKEVINKEDGKAKSGQNDDPRQAKCNWAGLNSPFDFAFITCTSKQPITEWMWHENVVDFWNAAFFWYTNQFGGPHTVGGNCKAQHWEPQKSVCMIKAHLYNVDDHDIRRRYEACTFKETVLAPAVNTIVQLQCDKPCTDPAPPWAYSVLLFQVNAWSWNKRRPLVGFPDQCHNNCKHHHMCHHNMAHHCCIEPIANCRCLDRWSLGVHFTWMGCGTAMWTTGTDAMIVWYVIWGDKWIRWNDDHLYYEYRSFQRIGAAKVQVQKRWHYHWVTIRRITSWVLWTHAAAHMCRNWRSFERRQWTVHNPWGCPVTWEGLHIMCQLEDHVCDSTEDGYNVFKRVYKYTCYKHRSECIEWLWDSLPPYQHINMQECRFQHCILNYNCYSCHGMQVKYVCFEDTELYNSERGATWWMWAWAYGGNKTQNMAWYIGGQNGHSPGGQCKNNIINLWTWYSDQMYCLPAMCNVVKDPSQPLMLGKPPPGVLTYDEASQLAEFAEHIYNKPWNINAQEMNMNMVRYYTCMNDTPKGPNQWTIEKTSPQDNAMCQYIVPSHQSWKKVTSKEW";
    int a = aminoacid1.length();
    int b = aminoacid2.length();
    vector <vector<int>> PAM250 = inputMatrix(20,20);
    Graph Project(aminoacid1,aminoacid2,PAM250);
    Project.FillGraph();
    /* PAM250

  2 -2  0  0 -3  1 -1 -1 -1 -2 -1  0  1  0 -2  1  1  0 -6 -3
 -2 12 -5 -5 -4 -3 -3 -2 -5 -6 -5 -4 -3 -5 -4  0 -2 -2 -8  0
  0 -5  4  3 -6  1  1 -2  0 -4 -3  2 -1  2 -1  0  0 -2 -7 -4
  0 -5  3  4 -5  0  1 -2  0 -3 -2  1 -1  2 -1  0  0 -2 -7 -4
 -3 -4 -6 -5  9 -5 -2  1 -5  2  0 -3 -5 -5 -4 -3 -3 -1  0  7
  1 -3  1  0 -5  5 -2 -3 -2 -4 -3  0  0 -1 -3  1  0 -1 -7 -5
 -1 -3  1  1 -2 -2  6 -2  0 -2 -2  2  0  3  2 -1 -1 -2 -3  0
 -1 -2 -2 -2  1 -3 -2  5 -2  2  2 -2 -2 -2 -2 -1  0  4 -5 -1
 -1 -5  0  0 -5 -2  0 -2  5 -3  0  1 -1  1  3  0  0 -2 -3 -4
 -2 -6 -4 -3  2 -4 -2  2 -3  6  4 -3 -3 -2 -3 -3 -2  2 -2 -1
 -1 -5 -3 -2  0 -3 -2  2  0  4  6 -2 -2 -1  0 -2 -1  2 -4 -2
  0 -4  2  1 -3  0  2 -2  1 -3 -2  2  0  1  0  1  0 -2 -4 -2
  1 -3 -1 -1 -5  0  0 -2 -1 -3 -2  0  6  0  0  1  0 -1 -6 -5
  0 -5  2  2 -5 -1  3 -2  1 -2 -1  1  0  4  1 -1 -1 -2 -5 -4
 -2 -4 -1 -1 -4 -3  2 -2  3 -3  0  0  0  1  6  0 -1 -2  2 -4
  1  0  0  0 -3  1 -1 -1  0 -3 -2  1  1 -1  0  2  1 -1 -2 -3
  1 -2  0  0 -3  0 -1  0  0 -2 -1  0  0 -1 -1  1  3  0 -5 -3
  0 -2 -2 -2 -1 -1 -2  4 -2  2  2 -2 -1 -2 -2 -1  0  4 -6 -2
 -6 -8 -7 -7  0 -7 -3 -5 -3 -2 -4 -4 -6 -5  2 -2 -5 -6 17  0
 -3  0 -4 -4  7 -5  0 -1 -4 -1 -2 -2 -5 -4 -4 -3 -3 -2  0 10

    */
   cout << endl;
    
    vector <vector<int>> x = Project.GetGraphTable();
    vector <vector<char>> y = Project.GetBackTrackTable();
    //OutputMatrix(x);
    cout << endl ;
    //OutputMatrix(y);
    
    cout << endl;
    pair <int,int> cord = Project.getpos(); 
    cout << Project.getmaxscore() << endl;
    BackTracking(y,aminoacid1,aminoacid2,cord);
    cout << endl;

    cout << "Hello World";
}
