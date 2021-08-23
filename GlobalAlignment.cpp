#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>

using namespace std;

//https://bioinformaticsalgorithms.com/data/realdatasets/Alignment/BLOSUM62.txt

/*

   A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
A  4  0 -2 -1 -2  0 -2 -1 -1 -1 -1 -2 -1 -1 -1  1  0  0 -3 -2
C  0  9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2
D -2 -3  6  2 -3 -1 -1 -3 -1 -4 -3  1 -1  0 -2  0 -1 -3 -4 -3
E -1 -4  2  5 -3 -2  0 -3  1 -3 -2  0 -1  2  0  0 -1 -2 -3 -2
F -2 -2 -3 -3  6 -3 -1  0 -3  0  0 -3 -4 -3 -3 -2 -2 -1  1  3
G  0 -3 -1 -2 -3  6 -2 -4 -2 -4 -3  0 -2 -2 -2  0 -2 -3 -2 -3
H -2 -3 -1  0 -1 -2  8 -3 -1 -3 -2  1 -2  0  0 -1 -2 -3 -2  2
I -1 -1 -3 -3  0 -4 -3  4 -3  2  1 -3 -3 -3 -3 -2 -1  3 -3 -1
K -1 -3 -1  1 -3 -2 -1 -3  5 -2 -1  0 -1  1  2  0 -1 -2 -3 -2
L -1 -1 -4 -3  0 -4 -3  2 -2  4  2 -3 -3 -2 -2 -2 -1  1 -2 -1
M -1 -1 -3 -2  0 -3 -2  1 -1  2  5 -2 -2  0 -1 -1 -1  1 -1 -1
N -2 -3  1  0 -3  0  1 -3  0 -3 -2  6 -2  0  0  1  0 -3 -4 -2
P -1 -3 -1 -1 -4 -2 -2 -3 -1 -3 -2 -2  7 -1 -2 -1 -1 -2 -4 -3
Q -1 -3  0  2 -3 -2  0 -3  1 -2  0  0 -1  5  1  0 -1 -2 -2 -1
R -1 -3 -2  0 -3 -2  0 -3  2 -2 -1  0 -2  1  5 -1 -1 -3 -3 -2
S  1 -1  0  0 -2  0 -1 -2  0 -2 -1  1 -1  0 -1  4  1 -2 -3 -2
T  0 -1 -1 -1 -2 -2 -2 -1 -1 -1 -1  0 -1 -1 -1  1  5  0 -2 -2
V  0 -1 -3 -2 -1 -3 -3  3 -2  1  1 -3 -2 -2 -3 -2  0  4 -3 -1
W -3 -2 -4 -3  1 -2 -2 -3 -3 -2 -1 -4 -4 -2 -3 -3 -2 -3 11  2
Y -2 -2 -3 -2  3 -3  2 -1 -2 -1 -1 -2 -3 -1 -2 -2 -2 -1  2  7

*/
vector <vector<int>> BLOSUM62 = {{4  ,0 ,-2 ,-1 ,-2  ,0 ,-2 ,-1 ,-1 ,-1 ,-1 ,-2 ,-1 ,-1 ,-1  ,1  ,0  ,0 ,-3 ,-2},{0  ,9 ,-3 ,-4 ,-2 ,-3 ,-3 ,-1 ,-3 ,-1 ,-1 ,-3 ,-3 ,-3 ,-3 ,-1 ,-1 ,-1 ,-2 ,-2},{-2 ,-3  ,6  ,2 ,-3 ,-1 ,-1 ,-3 ,-1 ,-4 ,-3  ,1 ,-1  ,0 ,-2  ,0 ,-1 ,-3 ,-4 ,-3},{-1 ,-4 , 2 , 5 ,-3 ,-2 , 0 ,-3  ,1 ,-3 ,-2 , 0 ,-1  ,2  ,0  ,0 ,-1 ,-2 ,-3 ,-2},{-2 ,-2 ,-3 ,-3  ,6 ,-3 ,-1 , 0 ,-3 , 0 , 0 ,-3 ,-4 ,-3 ,-3 ,-2 ,-2 ,-1 , 1 , 3},{0 ,-3 ,-1 ,-2 ,-3 , 6 ,-2 ,-4 ,-2 ,-4, -3 , 0 ,-2 ,-2 ,-2  ,0 ,-2 ,-3 ,-2 ,-3},{-2 ,-3 ,-1 , 0 ,-1 ,-2 , 8 ,-3 ,-1 ,-3 ,-2 , 1, -2 , 0  ,0 ,-1 ,-2 ,-3 ,-2 , 2},{-1 ,-1 ,-3 ,-3  ,0 ,-4 ,-3  ,4 ,-3 , 2 , 1 ,-3 ,-3 ,-3 ,-3 ,-2 ,-1 , 3 ,-3 ,-1},{-1 ,-3 ,-1  ,1 ,-3 ,-2 ,-1 ,-3  ,5 ,-2 ,-1 , 0, -1 , 1 , 2 , 0 ,-1 ,-2 ,-3 ,-2},{-1, -1, -4, -3 , 0 ,-4 ,-3 , 2 ,-2 , 4 , 2 ,-3 ,-3 ,-2 ,-2 ,-2 ,-1 , 1 ,-2 ,-1},{-1 ,-1 ,-3 ,-2 , 0 ,-3 ,-2 , 1 ,-1 , 2 , 5 ,-2 ,-2 , 0 ,-1 ,-1 ,-1 , 1 ,-1 ,-1},{-2 ,-3 , 1 , 0 ,-3 , 0 , 1 ,-3  ,0 ,-3 ,-2 , 6 ,-2 , 0 , 0 , 1 , 0 ,-3 ,-4, -2},{-1 ,-3 ,-1 ,-1 ,-4 ,-2 ,-2 ,-3 ,-1 ,-3 ,-2 ,-2  ,7 ,-1 ,-2 ,-1 ,-1 ,-2 ,-4, -3},{-1 ,-3 , 0 , 2 ,-3 ,-2 , 0 ,-3 , 1, -2 , 0 , 0 ,-1 , 5 , 1 , 0 ,-1 ,-2, -2 ,-1},{-1 ,-3 ,-2 , 0 ,-3 ,-2 , 0 ,-3 , 2, -2, -1 , 0 ,-2 , 1 , 5 ,-1 ,-1 ,-3, -3 ,-2},{1 ,-1 , 0 , 0 ,-2 , 0 ,-1 ,-2 , 0 ,-2 ,-1 , 1, -1 , 0 ,-1 , 4 , 1 ,-2 ,-3 ,-2},{0 ,-1 ,-1 ,-1 ,-2 ,-2 ,-2 ,-1 ,-1, -1, -1 , 0 ,-1 ,-1 ,-1 , 1 , 5 , 0 ,-2 ,-2},{0 ,-1 ,-3 ,-2 ,-1 ,-3 ,-3 , 3 ,-2 , 1 , 1 ,-3 ,-2 ,-2 ,-3 ,-2 , 0 , 4 ,-3 ,-1},{-3 ,-2 ,-4 ,-3 , 1 ,-2 ,-2 ,-3 ,-3 ,-2 ,-1 ,-4 ,-4 ,-2 ,-3 ,-3 ,-2 ,-3 ,11 , 2},{-2 ,-2 ,-3 ,-2 , 3 ,-3 , 2 ,-1 ,-2 ,-1 ,-1 ,-2 ,-3 ,-1 ,-2 ,-2, -2, -1 , 2 , 7}};



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

void BackTracking(vector <vector<char>> Matrix,string firstAminoAcid,string secondAminoAcid)
{
   int rows = Matrix.size();
   int columns = Matrix[0].size();
   vector<char> c,d;
   int a = firstAminoAcid.length()-1;
   int b = secondAminoAcid.length()-1;
    // 0 is diagonal , 1 Horizontal , 2 Vertical;
   for (int i = rows-1; i >= 0 ; )
   {
      for (int j = columns-1; j >=0 ;)
      {
         char f = Matrix[i][j];
         if(Matrix[i][j] == 'D')
         {
            c.push_back(firstAminoAcid[a]);
            d.push_back(secondAminoAcid[b]);
            a--;
            b--;
            i--;
            j--;
            continue;
         }
         else if(Matrix[i][j] == 'V')
         {
            i--;
            c.push_back(firstAminoAcid[a]);
            d.push_back('-');
            a--;
            continue;
         }
         else if(Matrix[i][j] == 'H')
         {
            j--;
            c.push_back('-');
            d.push_back(secondAminoAcid[b]);
            b--;
            continue;
         }
         if (i==0 && j==0)
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

// Vertical (columns) is First Amino acid , Horizontal is the second amino acid
class Graph
{
private:
   string FirstAminoAcid;
   string SecondAminoAcid;
   int FirstAminoAcidLength;
   int SecondAminoAcidLLength;
   vector <vector <int>> GraphTable;
   vector<vector<char>> Track;
public:
   
   Graph(string x , string y)
   {
      FirstAminoAcid = x;
      SecondAminoAcid = y;
      FirstAminoAcidLength = x.length();
      SecondAminoAcidLLength = y.length(); 
   }

   void CreateGraph()
   {
      vector <int> init(SecondAminoAcidLLength+1,0);
      vector <char> ini(SecondAminoAcidLLength+1,'X');
      for (int i = 0; i <= FirstAminoAcidLength; i++)
      {
         GraphTable.push_back(init);
         Track.push_back(ini);
      }

      for (int i = 1; i < FirstAminoAcidLength+1 ; i++) // Gaps on Vertical Axis
      {
         GraphTable[i][0] = GraphTable[i-1][0] - 5;
         Track[i][0] = 'V';
      }

      for (int i = 1; i < SecondAminoAcidLLength+1; i++)
      {
         GraphTable[0][i] = GraphTable[0][i-1]-5;
         Track[0][i] = 'H';
      }
   }
   
   void FillGraph()
   {
      for (int i = 1; i <= FirstAminoAcidLength; i++)
      {
         for (int j = 1; j <= SecondAminoAcidLLength; j++)
         {
            int vals [3];
            vals[0] = GraphTable[i-1][j]-5;
            vals[1]= GraphTable[i][j-1] -5;
            int firstletterpos = GetPos(FirstAminoAcid[i-1]);
            int secondletterpos = GetPos(SecondAminoAcid[j-1]);
            int blval = BLOSUM62[firstletterpos][secondletterpos];
            vals[2]= GraphTable[i-1][j-1] + blval;

            int maxx = -5000;
            for (int x = 0; x < 3; x++)
            {
               if(maxx < vals[x])
               {
                  maxx = vals[x];
               }
            }
            GraphTable[i][j] = maxx;
            if(maxx == vals[0])
            {
               Track[i][j] = 'V';

            }
            else if(maxx == vals[1])
            {
               Track[i][j] = 'H';
            }
            else if(maxx == vals[2])
            {
               Track[i][j] = 'D';
            }
            
         }
         
      }
      

   }
   string GetFAA()
   {
      return FirstAminoAcid;
   }
   string GetSAA()
   {
      return SecondAminoAcid;
   }
   int GetFAAL()
   {
      return FirstAminoAcidLength;
   }
   int GetSAAL()
   {
      return SecondAminoAcidLLength;
   }
   vector <vector<int>> GetGraphTable()
   {
      return GraphTable;
   }
   vector<vector<char>> GetTrack()
   {
      return Track;
   }


};





int main()
{
   string c = "VYEHNTHSEEKRKGHQAYWTCWPYQESNHPMEREISRGTNTEAGWLNNYFCVVDFCVHVNKDRKGTFDMQMFRYIFYWEDGICEFVANLMGILCQKRCMTYSHSIPANFVFSRPEALEEYAHVNPPSKPGSRHKNKKNAPIQAVFNWWSREPNEVAWEQKGLKPTSVEKYRGPLYFACARHEYVMNNMMLQFMKIDCPAKAADYWTIDGYDLWNFEHYWCIMVAISPASAFKFRAARGAVELFCWVVWLYNKFRHPGDCSYIVTQVQQRLPVCRLIVWWFNHNQIRFVQWHSQKVGPVNSADYIQIHHDQFSSFGRNWPNGLWHMYDMRLFGELNIWHNQEVDADLQWVDQKNKMPSLCQVWHDAPLLQRDRNNTMWLRKMDCICKENMMNWLQKDPRNSGMHIKIFAVCEMQPWTWMTEWTYFGEMYNHAGMSLNDCPRTPVHMLHEDVRCPQKTMNTQHSKDKIIWAQFWCWWGWEEYFRMGSRAENISSNWDSMPVGEAHMDYCEHTVCVHHRWHIRDDMRLHQDRNIQYVGQWSDIMWFFPIWKEKLKSRYYAFANQSQSDKMGSIYMPVPQMREQKRKYRHAQWTDNHTFVNWHWTAPFLDHNSQGMLSELRPNLYVFMLMGIHLDAWPVSTYHDVILETVHRIQQWHKQTMGLKPKTSKHYTFIQPFGSTHCPIYDMCFRCISAWNRYVCQCYYIDELRETQEHANRQCLKRGNPHYVQPMSM";
   string y = "VYEHQTHSEEKRKHHQHGRNDGYMLRPPNGHGEATYLSCVSRGTMDAWLSGPTSKFHNYFCVVDTAPKFKVHFCVHFHIKDRKGTFDMQMWRYIFYWEMITHGIINSIICEVANLMGILCQNSAAWTCRCCLPMHDIHVHDEGHPYKLEEYAHVMSPSKPGHKFTAKKNAPVFNWWEENMVVREPCEWAQQWWYKCAWEQQGLKPTSVEKYRDWNMNGNTVNILYFACAEYMMLQCMKLDCPAKALWNDLDTEHYHCIMKAISFKFRAARGAVELFCWFVWLYNKNRHTQVQQWLPVCRLWWFNHNQIRFVQSHKVGPPGGSNSADYIQIHSHICWGDQFSSAGRCQLWHMYDLKTSVSRRRLFGELVDCWHNQADLQWVIAAPQKNFMCNVGIDCMPSLCQVWHDAPLLQRDRNNTMWLRKMDCICKWHAKYWFQTWWNSCICQQKDPRNSGHIKIFKPFPWGMREWLQKAEDYFGETCQGWNILGMSLNDMPRTFPYQIRPGSNHEDVREPQMNTQHSKDKIIWDADISGYATWLWFEYFRMGLRAENISSNEGSHGVYVYSSEDVSHHNPVGEAHPLPDYCEHTDACVACVYHKRDDMRYVGQWNDIQWFKLKSLPYAFANVSWYMQYYGSHYMRKYRHYHIRIGFEIQWTDNHTFVFLDHNSQYMLSMGPYDDKDSLRPNLYVGIHLDAWPVSTYHDPIQWHKQTMGLKPKTSKHYTFIQPFKYDMCFRCISDWNRYVCQCYYIDELQCAKRGNPHYVQPWSM";
   Graph x( c, y);
   x.CreateGraph();
   x.FillGraph();
   vector<vector<int>> Gr = x.GetGraphTable();//Matrix With Score
   vector<vector<char>> Tr = x.GetTrack();// Trackback Matrix
   int rows = Gr.size();
   int columns = Gr[0].size();


   /*for (int i = 0; i < rows; i++) // Function used in Debugging To Output Score Matrix
   {
      for (int j = 0; j < columns; j++)
      {
         cout << setw(5) << Gr[i][j] << " ";
      }
      cout << endl;
      
   }*/

   cout << endl;

   /*for (int i = 0; i < rows; i++) // Function used in Debugging To Output Trackback Matrix
   {
      for (int j = 0; j < columns; j++)
      {
         cout << setw(5) << Tr[i][j] << " ";
      }
      cout << endl;
      
   }*/

   cout << Gr[rows-1][columns-1] << endl;
   BackTracking(Tr,c,y);

 



   
   
   

   
}

/*

Input :

VYEHNTHSEEKRKGHQAYWTCWPYQESNHPMEREISRGTNTEAGWLNNYFCVVDFCVHVNKDRKGTFDMQMFRYIFYWEDGICEFVANLMGILCQKRCMTYSHSIPANFVFSRPEALEEYAHVNPPSKPGSRHKNKKNAPIQAVFNWWSREPNEVAWEQKGLKPTSVEKYRGPLYFACARHEYVMNNMMLQFMKIDCPAKAADYWTIDGYDLWNFEHYWCIMVAISPASAFKFRAARGAVELFCWVVWLYNKFRHPGDCSYIVTQVQQRLPVCRLIVWWFNHNQIRFVQWHSQKVGPVNSADYIQIHHDQFSSFGRNWPNGLWHMYDMRLFGELNIWHNQEVDADLQWVDQKNKMPSLCQVWHDAPLLQRDRNNTMWLRKMDCICKENMMNWLQKDPRNSGMHIKIFAVCEMQPWTWMTEWTYFGEMYNHAGMSLNDCPRTPVHMLHEDVRCPQKTMNTQHSKDKIIWAQFWCWWGWEEYFRMGSRAENISSNWDSMPVGEAHMDYCEHTVCVHHRWHIRDDMRLHQDRNIQYVGQWSDIMWFFPIWKEKLKSRYYAFANQSQSDKMGSIYMPVPQMREQKRKYRHAQWTDNHTFVNWHWTAPFLDHNSQGMLSELRPNLYVFMLMGIHLDAWPVSTYHDVILETVHRIQQWHKQTMGLKPKTSKHYTFIQPFGSTHCPIYDMCFRCISAWNRYVCQCYYIDELRETQEHANRQCLKRGNPHYVQPMSM
VYEHQTHSEEKRKHHQHGRNDGYMLRPPNGHGEATYLSCVSRGTMDAWLSGPTSKFHNYFCVVDTAPKFKVHFCVHFHIKDRKGTFDMQMWRYIFYWEMITHGIINSIICEVANLMGILCQNSAAWTCRCCLPMHDIHVHDEGHPYKLEEYAHVMSPSKPGHKFTAKKNAPVFNWWEENMVVREPCEWAQQWWYKCAWEQQGLKPTSVEKYRDWNMNGNTVNILYFACAEYMMLQCMKLDCPAKALWNDLDTEHYHCIMKAISFKFRAARGAVELFCWFVWLYNKNRHTQVQQWLPVCRLWWFNHNQIRFVQSHKVGPPGGSNSADYIQIHSHICWGDQFSSAGRCQLWHMYDLKTSVSRRRLFGELVDCWHNQADLQWVIAAPQKNFMCNVGIDCMPSLCQVWHDAPLLQRDRNNTMWLRKMDCICKWHAKYWFQTWWNSCICQQKDPRNSGHIKIFKPFPWGMREWLQKAEDYFGETCQGWNILGMSLNDMPRTFPYQIRPGSNHEDVREPQMNTQHSKDKIIWDADISGYATWLWFEYFRMGLRAENISSNEGSHGVYVYSSEDVSHHNPVGEAHPLPDYCEHTDACVACVYHKRDDMRYVGQWNDIQWFKLKSLPYAFANVSWYMQYYGSHYMRKYRHYHIRIGFEIQWTDNHTFVFLDHNSQYMLSMGPYDDKDSLRPNLYVGIHLDAWPVSTYHDPIQWHKQTMGLKPKTSKHYTFIQPFKYDMCFRCISDWNRYVCQCYYIDELQCAKRGNPHYVQPWSM


Output :

1567
VYEHNTHSEEKRKGHQ-AYWTCWPYQESN-H-PMER-E-ISRGT-NTE-AG-W--LNNYFCVVD----F----CVHVN-KDRKGTFDMQMFRYIFYWE--D-GI-----CEFVANLMGILCQKRCM-TYSHSIPAN-F-VFSR--PEALEEYAHVNPPSKPGSRHKNKKNAPIQAVFNWW--S---REPNE-V-------AWEQKGLKPTSVEKYR-----G-P---LYFACARHEYVMNNMMLQFMKIDCPAKAADYWTIDGYDLWNFEHYWCIMVAISPASAFKFRAARGAVELFCWVVWLYNKFRHPGDCSYIVTQVQQRLPVCRLIVWWFNHNQIRFVQWHSQKVGPV---NSADYIQIH-H----DQFSSFGRNWPNGLWHMYDM-----R--LFGEL-NIWHNQEVDADLQWV---DQKN---K-----MPSLCQVWHDAPLLQRDRNNTMWLRKMDCICK-E-N--MMNW----L-Q-KDPRNSGMHIKIFAVCE--MQPWTWMTEWTYFGE--M-YNHAGMSLNDCPRT-P--VHM-L-HEDVRCPQKTMNTQHSKDKIIW-AQF--WCWWGWEEYFRMGSRAENISSN------W-----D-SM--PVGEAH-M-DYCEHT-VCVHHRWHIRDDMRLHQDRNIQYVGQWSDIMWFFPIWKEKLKSRYYAFANQS-QSDKMGSIYMPVPQMREQK-R-KYRHAQWTDNHTFVNW-HWTAPF-LDHNSQGMLSELRPNLYVFMLMGIHLDAWPVSTYHDVILETVHRIQQWHKQTMGLKPKTSKHYTFIQPFGSTHCPIYDMCFRCISAWNRYVCQCYYIDELRETQEHANRQCLKRGNPHYVQPMSM
VYEHQTHSEEKRKHHQHGRNDGYMLRPPNGHGEATYLSCVSRGTMDAWLSGPTSKFHNYFCVVDTAPKFKVHFCVHFHIKDRKGTFDMQMWRYIFYWEMITHGIINSIICE-VANLMGILCQNSAAWTCRCCLPMHDIHVHDEGHPYKLEEYAHVMSPSKPGHKFTAKKNAP---VFNWWEENMVVREPCEWAQQWWYKCAWEQQGLKPTSVEKYRDWNMNGNTVNILYFACA--EY-M--M-LQCMKLDCPAKAL--WN-D-LD--T-EHYHCIMKAIS----FKFRAARGAVELFCWFVWLYNKNRH--------TQVQQWLPVCRL--WWFNHNQIRFVQSH--KVGPPGGSNSADYIQIHSHICWGDQFSSAGR-C-Q-LWHMYDLKTSVSRRRLFGELVDCWHNQ---ADLQWVIAAPQKNFMCNVGIDCMPSLCQVWHDAPLLQRDRNNTMWLRKMDCICKWHAKYWFQTWWNSCICQQKDPRNSG-HIKIFKPFPWGMREWLQKAE-DYFGETCQGWNILGMSLNDMPRTFPYQIRPGSNHEDVREPQ--MNTQHSKDKIIWDADISGYATWLWFEYFRMGLRAENISSNEGSHGVYVYSSEDVSHHNPVGEAHPLPDYCEHTDACVACVYHKRDDMR--------YVGQWNDIQWF----K--LKSLPYAFANVSWYMQYYGSHYMR--KYRHYHIRIGFE-IQWTDNHTFVFLDH-NSQYMLSMGPYDDKDSLRPNLYV----GIHLDAWPVSTYHD---P----IQ-WHKQTMGLKPKTSKHYTFIQPF-K-----YDMCFRCISDWNRYVCQCYYIDEL---Q------CAKRGNPHYVQPWSM

*/