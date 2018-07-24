/*
To perform paper "Broad cross-reactivity of the T-cell repertoire achieves specific and sufficiently rapid target searching." (https://arxiv.org/abs/1712.04633)
Author: Jin Xu and Junghyo Jo
*/



#include <map>
#include <vector>
#include <cstdlib>
#include <iostream> /*cin, cout, endl*/
#include <stdlib.h> /*srand(),rand()*/
#include <fstream> /*ofstream os*/
#include <stdio.h> /*printf, srand, rand*/
#include <string>
#include <iomanip> /*setw*/
#include <cmath> /*abs(), sqrt(), pow()*/
#include <time.h> /*time*/
using namespace std; /*cin, cout, endl*/

/*
Fisher-Yates shuffle   http://ideone.com/3A3cv
*/
#include <vector>
#include <algorithm>
#include <cstdlib>

template<class bidiiter>
bidiiter random_unique(bidiiter begin, bidiiter end, size_t num_random)
{
    size_t left = std::distance(begin, end);
    while (num_random--)
    {
        bidiiter r = begin;
        std::advance(r, rand()%left);
        std::swap(*begin, *r);
        ++begin;
        --left;
    }
    return begin;
}


/*global variables*/


int main()

{
     

    double aa[20][20];                  // a.a. contact energy in M-J matrix
   
    
    int h,p,k;
    int i,j,x,y,ii;

    int min,max;
    int direct;
   
    static int tcr[1000000][20]={0};    // record tcr CDR3 a.a. sequence with total max number of tcrs 10^6 and max length 20
    static int tcr_reverse[1000000][20]={0};
    static int pep[100000][200]={0};    // record peptide a.a. sequence with total max number of peptides 10^5 and max length 200

            
    string newSequence_file,newSequence_TCR,newSequence_peptide;
    char bb;
  
    static int tcr_length[1000000]={0};  // record the length of the certain tcr
    static int pep_length[100000]={0};   // record the length of the certain peptide

    int index_x, index_y;                // index for AA position
    int index_move_start_max;
    int tcr_move[300];
    int tcr_reverse_move[300];
    int pep_fix[300];

    double pair_E_new[230];
    int contact_len[230];
    double pair_reverse_E_new[230];
    int reverse_contact_len[230];
    int T;
    int P; // count the actual P(peptide) number recognized/responded tor
    int pair_len;
    double pair_E;
    int pre_T;
   
    double expscantimej;

    double Ec=33.;
    double Er;
    double Javg=4.;

    double scale=1.;

    //output: file and formats

    fstream file_TCR;
    file_TCR.open("TCR-12-natureCR.txt");	// read a.a. frequency from this file
    //output: format
    file_TCR.precision(6);
    file_TCR.setf(ios::scientific | ios::showpoint);
    cout.precision(6);
    cout.setf(ios::fixed | ios::showpoint);

    fstream file_peptide;
    file_peptide.open("infec-peptide-15-sample.txt");  // read a.a. frequency from this file
    //output: format
    file_peptide.precision(6);
    file_peptide.setf(ios::scientific | ios::showpoint);
    cout.precision(6);
    cout.setf(ios::fixed | ios::showpoint);


    fstream file_MJ;
    file_MJ.open("MJ-matrix.txt"); // read Miyazawa and Jernigan contact energies in RT units from this file
    //output: format
    file_MJ.precision(6);
    file_MJ.setf(ios::scientific | ios::showpoint);
    cout.precision(6);
    cout.setf(ios::fixed | ios::showpoint);


    ofstream file;
    file.open("scanningTime.txt");      
    //output: format
    file.precision(6);
    file.setf(ios::scientific | ios::showpoint);
    cout.precision(6);
    cout.setf(ios::fixed | ios::showpoint);


//Read MJ-matrix

    for (i=0; i<20; i++)
    { 
	file_MJ >> aa[i][0]  >> aa[i][1]  >> aa[i][2]  >> aa[i][3]  >> aa[i][4]  >> aa[i][5]  >> aa[i][6]  >> aa[i][7]  >> aa[i][8]  >> aa[i][9]
                >> aa[i][10] >> aa[i][11] >> aa[i][12] >> aa[i][13] >> aa[i][14] >> aa[i][15] >> aa[i][16] >> aa[i][17] >> aa[i][18] >> aa[i][19];
    }    




// read TCR file into max tcr[h][k], h is the total number of TCRs, tcr_length[h] record the length of the certain TCR with index h 

    h=0;  // count the total number of tcrs
    if(file_TCR.fail())
    {
        cout <<"error.\n";
    }

    else
    {
     	string temp;
        while (!file_TCR.eof())
        {
            getline(file_TCR,temp);
            newSequence_TCR.append(temp);
           
            while(!newSequence_TCR.empty())
            {
             	
		//cout << newSequence_TCR << endl;
		tcr_length[h]= newSequence_TCR.length();
		//printf("%d\n",tcr_length[h]);
		for(k = 0; k < tcr_length[h]; k++)
		{
		    bb=newSequence_TCR[k];       
		    
		    if ( bb=='C' ) tcr[h][k]=0;
		    else if ( bb=='M' ) tcr[h][k]=1;
		    else if ( bb=='F' ) tcr[h][k]=2;
		    else if ( bb=='I' ) tcr[h][k]=3;
		    else if ( bb=='L' ) tcr[h][k]=4;
		    else if ( bb=='V' ) tcr[h][k]=5;
		    else if ( bb=='W' ) tcr[h][k]=6;
		    else if ( bb=='Y' ) tcr[h][k]=7;
		    else if ( bb=='A' ) tcr[h][k]=8;
		    else if ( bb=='G' ) tcr[h][k]=9;
		    else if ( bb=='T' ) tcr[h][k]=10;
		    else if ( bb=='S' ) tcr[h][k]=11;
		    else if ( bb=='N' ) tcr[h][k]=12;
		    else if ( bb=='Q' ) tcr[h][k]=13;
		    else if ( bb=='D' ) tcr[h][k]=14;
		    else if ( bb=='E' ) tcr[h][k]=15;
		    else if ( bb=='H' ) tcr[h][k]=16;
		    else if ( bb=='R' ) tcr[h][k]=17;
		    else if ( bb=='K' ) tcr[h][k]=18;
		    else if ( bb=='P' ) tcr[h][k]=19;
		    else if ( bb=='-' ) tcr[h][k]=20;
		    //printf ("%d\t",tcr[h][k]);
		}
                //printf ("\n");

		h++;
		newSequence_TCR = "";
            }
	    
        }
    }

    //cout << h << endl;


// read peptide file into matrix pep[p][k], p is the total number of peptides, pep_length[p] record the length of the certain peptide with index p

    p=0;  // count the total number of peptides
    if(file_peptide.fail())
    {
     	cout <<"error.\n";
    }

    else
    {
     	string temp;
        while (!file_peptide.eof())
        {
            getline(file_peptide,temp);
            newSequence_peptide.append(temp);

            while(!newSequence_peptide.empty())
            {

                //cout << newSequence_peptide << endl;
                pep_length[p]= newSequence_peptide.length();
                //printf("%d\n",pep_length[p]);
                for(k = 0; k < pep_length[p]; k++)
                {
                    bb=newSequence_peptide[k];

                    if ( bb=='C' ) pep[p][k]=0;
                    else if ( bb=='M' ) pep[p][k]=1;
                    else if ( bb=='F' ) pep[p][k]=2;
                    else if ( bb=='I' ) pep[p][k]=3;
                    else if ( bb=='L' ) pep[p][k]=4;
                    else if ( bb=='V' ) pep[p][k]=5;
                    else if ( bb=='W' ) pep[p][k]=6;
                    else if ( bb=='Y' ) pep[p][k]=7;
                    else if ( bb=='A' ) pep[p][k]=8;
                    else if ( bb=='G' ) pep[p][k]=9;
                    else if ( bb=='T' ) pep[p][k]=10;
                    else if ( bb=='S' ) pep[p][k]=11;
                    else if ( bb=='N' ) pep[p][k]=12;
                    else if ( bb=='Q' ) pep[p][k]=13;
                    else if ( bb=='D' ) pep[p][k]=14;
                    else if ( bb=='E' ) pep[p][k]=15;
                    else if ( bb=='H' ) pep[p][k]=16;
                    else if ( bb=='R' ) pep[p][k]=17;
                    else if ( bb=='K' ) pep[p][k]=18;
                    else if ( bb=='P' ) pep[p][k]=19;
                    //printf ("%d\t",pep[p][k]);
                }
		//printf ("\n");

                p++;
                newSequence_peptide = "";
            }

	}
    }

    //cout << p << endl;
    srand((unsigned)time(NULL));


////calculate the binding energy for each pairwise pair[h][p]

//forward binding & reverse binding

    for (i=0; i<h; i++)
    {

        // do tcr reverse
        for(x=0;x<tcr_length[i];x++)
        {
            tcr_reverse[i][x]=tcr[i][tcr_length[i]-x-1];
        }
	
        
    }

   


    for (j=0; j<p; j++)
    {

	
	expscantimej=0.;
	
	T=0;
	pre_T=0;

	max=pep_length[j];

	
	std::vector<int> index_T(h);      
        for(int e=0; e<h; ++e)
        {
            index_T[e] = e;
        }
       	random_unique(index_T.begin(), index_T.end(), h);
        for(int f=0; f<h; ++f)
        {
            
        
	    i=index_T[f];
	    

            min=tcr_length[i];
	
            pair_E=0.;
            pair_len=0;
	    
            
            min=tcr_length[i];
            

            for(k=min;k<(max+1);k++)
            {
             	pair_E_new[k]=0.;
                
                contact_len[k]=0;
                
            }

            
            k=min + (rand()%(max-min+1));
           

            direct=rand()%2;
            


            for (y=0;y<pep_length[j];y++)
            {
             	index_y=y+tcr_length[i];
                pep_fix[index_y]=pep[j][y];

               

                if(direct==0)
                {

                    for(x=0;x<tcr_length[i];x++)
                    {

                     	index_x=x+k;
                        tcr_move[index_x]=tcr[i][x];

                        if (index_x==index_y)
                       	{
                            pair_E_new[k]+=aa[tcr_move[index_x]][pep_fix[index_y]];
                           
                            contact_len[k]++;
                        }
                    }
                }

                if(direct==1)
                {
                    for(x=0;x<tcr_length[i];x++)
                    {

                     	index_x=x+k;

                        //reverse
                        tcr_reverse_move[index_x]=tcr_reverse[i][x];
                        if (index_x==index_y)
                        {
                            pair_E_new[k]+=aa[tcr_reverse_move[index_x]][pep_fix[index_y]];
                            contact_len[k]++;
                        }

                    }
		}
                    
            }//y:pep_length
  

            pair_E=pair_E_new[k];
            pair_len=contact_len[k];
            

	    pair_E+=Ec;
	    Er=Ec+pair_len*Javg;

            expscantimej+=exp(scale*pair_E);

            pre_T++;

            if (pair_E>Er)
            {
		T+=1;
		f=h;
		file << setw(15) << j << setw(15) << pre_T << setw(15) << expscantimej;
		pre_T--;
		expscantimej-=exp(scale*pair_E);
                file << setw(15) << pre_T << setw(15) << expscantimej << endl;      
		
            }
            //printf("%d\n",h);
            

        }//f
        

    }//j

   
    
    return 0;


}
