#include<iostream>
#include<string>
#include<algorithm>
#include<vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <time.h>


using namespace std;

#define TABLE_SIZE 251

class Point {
private:
    
public:
    string name;
    string sequence;
    Point(string title = "Default", string seq = "")
    {
        name = title;
        sequence = seq;
    }
    
};
class Pointr {
private:

public:
	string name;
	int sequence;
	Pointr(string title = "Default", int rand=0)
	{
		name = title;
		sequence = rand;
	}

};

class Point2 {
private:

public:
	string name1;
	string name2;
	Point2(string n1 = "Default", string n2 = "Default")
	{
		name2 = n2;
		name1 = n1;
	}

};

class Point3 {
private:

public:
	string firsta;
	string secondb;
	double score;
	Point3(string a = "Default", string b = "Default", double scorea = 0)
	{
		firsta = a;
		secondb = b;
		score = scorea;
	}

};


double calcLengthDiff(double x, double y) {
	double val;
	val = x - y;
	return val;
}

int factorialfinder(int x)
{
	if (x == 1)        // HERE 5 is not equal to 1 so goes to else
	{
		return 1;
	}
	else
	{
		return x*factorialfinder(x - 1); // returns 5*4*3*2*1  when x==1 it returns 1
	}
}



string PhyloTreeAlgorithm(vector<Point> Seqs)
{
	vector<double> scorefinal;
	vector<Point3> scoresP;
	int cases = Seqs.size();
	int switched = 0;
	for (int i = 0; i < cases - 1; i++)
	{
		for (int j = i + 1; j < cases; j++)
		{
			if (Seqs[j].sequence.size() < Seqs[i].sequence.size()) {
				Point temp = Seqs[i];
				Seqs[i] = Seqs[j];
				Seqs[j] = temp;
				switched = 1;
			}


			double table[TABLE_SIZE][TABLE_SIZE];
			for (int p = 0; p < TABLE_SIZE; p++)
			{
				for (int q = 0; q < TABLE_SIZE; q++)
				{
					table[p][q] = 0;
				}
			}

			// go through entire string

			for (int m = 0; m < Seqs[i].sequence.size(); m++)
			{
				for (int n = 0; n < Seqs[j].sequence.size(); n++) //compare with entire other string
				{
					if (Seqs[i].sequence[m] == Seqs[j].sequence[n]) //if equal, score is 1
						table[m][n] = 1;
					else if ((Seqs[i].sequence[m] == 'A' || Seqs[i].sequence[m] == 'G') && (Seqs[j].sequence[n] == 'G' || Seqs[j].sequence[n] == 'A')) //if one is A and other is G (or vice verca)
						table[m][n] = -0.66;
					else if ((Seqs[i].sequence[m] == 'T' || Seqs[i].sequence[m] == 'C') && (Seqs[j].sequence[n] == 'C' || Seqs[j].sequence[n] == 'T')) //if one is A and other is G (or vice verca)
						table[m][n] = -0.66;
					else
						table[m][n] = -1;
				}
			}


			for (int p = 0; p < TABLE_SIZE; p++)
			{
				for (int q = 0; q < TABLE_SIZE; q++)
				{
					cout << setw(4) << table[p][q] << " ";
				}
				cout << endl;
			}
			cout << endl;


			int m = 0, n = 0;
			double score = table[m][n];
			while (m < TABLE_SIZE - 1 && n < TABLE_SIZE - 1)
			{
				//3 cases, go right if values equal or if right is better, go diagnol if diagnol is better
				if (table[m][n + 1] <= table[m + 1][n + 1])
				{
			
					score += table[m + 1][n + 1];
					m++;
					n++;
				}
				else
				{
					score += table[m][n + 1];
					n++;
				}
			}
			if (switched == 1) {
				Point temp = Seqs[i];
				Seqs[i] = Seqs[j];
				Seqs[j] = temp;
				switched = 0;
			}


			scorefinal.push_back(score);
		
			scoresP.push_back(Point3(Seqs[i].name, Seqs[j].name, score));
		}
	}
	double maxscore = -999999;
	double minscore = 999999;

	/*int answers = factorialfinder(cases - 1); //// this is old score output format
	for (int i = 0; i < answers; i++) {
		cout << scorefinal[i];
		cout << endl;
	}*/
	cout << "new score format" << endl;
	for (int i = 0; i < scoresP.size(); i++) {
		if (scoresP[i].score < minscore) {//min scoresss
			minscore = scoresP[i].score;
		}
		if (scoresP[i].score > maxscore) {//max score setting
			maxscore = scoresP[i].score;
		}
		cout << scoresP[i].firsta << "  " << scoresP[i].secondb << " " << scoresP[i].score;
		cout << endl;
	}
	
	/*
	vector <Pointr> tier;
	int rand = 111;
	string str = scoresP[0].firsta;
	tier.push_back(Pointr(str, rand));
	tier.push_back(Pointr(scoresP[0].secondb, rand));

	for (int i = 1; i < cases - 1; i++) {
		for (int j = 0; j < tier.size(); j++) {
			for (int k = 0; k < tier.size(); k++) {
				if (scoresP[i].firsta == tier[j].name) {
					if (scoresP[i].secondb == tier[k].name && tier[j].name == tier[k].name) {
						scoresP.erase(scoresP.begin() + i);
					}
					else tier.push_back(Pointr(scoresP[0].secondb, tier[j].sequence));

					j = tier.size();
					k = tier.size();
				}
				else {
					if (scoresP[i].secondb == tier[j].name) {
						tier.push_back(Pointr(scoresP[0].firsta, tier[j].sequence));
					}
					else {
						int rand2 = j;
						tier.push_back(Pointr(scoresP[0].firsta, rand2));
						tier.push_back(Pointr(scoresP[0].secondb, rand2));
						j = tier.size();
						k = tier.size();
					}

				}

			}
		}

	}

	cout << "deleted val score format" << endl;
	for (int i = 0; i < scoresP.size(); i++) {
		if (scoresP[i].score < minscore) {//min scoresss
			minscore = scoresP[i].score;
		}
		if (scoresP[i].score > maxscore) {//max score setting
			maxscore = scoresP[i].score;
		}
		cout << scoresP[i].firsta << "  " << scoresP[i].secondb << " " << scoresP[i].score;
		cout << endl;
	}

	*/

	cout << endl;
	cout << "Decimal Values: " << endl;
	//////// tree formatting

	///changing to decimal
	double totaldist = maxscore - minscore;

	for (int i = 0; i < scoresP.size(); i++) {
		if (totaldist != 0) {
			scoresP[i].score = (totaldist - scoresP[i].score) / totaldist; //replaces scores with decimal value
		}
		else {
			scoresP[i].score = 1;
		}
	}

	///sorting  bubble sort
	for (int i = 0; i < scoresP.size(); i++) {
		for (int j = 0; j < scoresP.size(); j++) {
			if (scoresP[i].score < scoresP[j].score) {

				string temp = scoresP[i].firsta;
				string temp2 = scoresP[i].secondb;
				scoresP[i].firsta = scoresP[j].firsta;
				scoresP[i].secondb = scoresP[j].secondb;
				scoresP[j].firsta = temp;
				scoresP[j].secondb = temp2;

				double temp1 = scoresP[i].score;
				scoresP[i].score = scoresP[j].score;
				scoresP[j].score = temp1;
			}
		}
	}

	for (int j = 0; j < scoresP.size(); j++) {
		cout << scoresP[j].firsta << "  " << scoresP[j].secondb << " " << scoresP[j].score << endl;
	}	///// all scores are sorted


	vector<Point2> connected;
	vector<double> lengths;
	vector<string> vals;
	//Always 1st leaf
	vals.push_back(scoresP[0].firsta);
	vals.push_back(scoresP[0].secondb);

	for (int j = 0; j < cases - 1; j++) {
		for (int i = 0; i < vals.size(); i++) {
			if (scoresP[j].firsta == vals[i] || scoresP[j].secondb == vals[i]) {
				double length = calcLengthDiff(scoresP[cases - 2].score, scoresP[i].score);
				if (length <= 0) {
					i = vals.size();
					break;
				}
				else {
					lengths.push_back(length);
				}
			}
			else {
				vals.push_back(scoresP[j].firsta);
				vals.push_back(scoresP[j].secondb);
			}
		}
	}
	//if (vals[2] == vals[1]) {
	//	vals[2] = vals[3];
	//}
	for (int i = 0; i < vals.size(); i++) {
		cout << "VALS: " << vals[i] << " " << endl;
	}
	for (int i = 0; i < lengths.size(); i++) {
		cout << "LENGTHS: " << lengths[i] << " " << endl;
	}

	// Apply data into newick format
	vector<string> treenodes;
	vector<string> vals2;

	// connects the points
	int tree = 0;
	for (int i = 0; i < cases-1; i++) {
		
		string tempa = scoresP[i].firsta;
		string tempb = scoresP[i].secondb;
		int check = 0;
		for (int p = 0; p < 2 * i; p++) {//works 4 even seq at (2*1) other wise 2*i -2
			for (int q = 0; q < 2 * i ; q++) {
				if (vals.size() > 3) {
					if (scoresP[i].firsta == vals[p] && scoresP[i].secondb == vals[q]) { // cases - 2 for working 4-2 seqs
						string addnode = "(" + treenodes[tree] + ":" + to_string(lengths[tree]) + ", " + treenodes[tree + 1] + ":" + to_string(lengths[tree + 1]) + ")"; // catanating to the previous nodes??
						cout << addnode << endl;
						treenodes.push_back(addnode);
						p = 2 * i;
						q = 2 * i;
						check = 1;
						tree++;
					}
				}
			}
		}
		for (int p = 0; p < 2 * i; p++) {//// if check==0 outside the double for'ss
			for (int q = 0; q < 2 * i; q++) {
			
				if (scoresP[i].firsta == vals[q] && check == 0) {
					string addnode = "(" + treenodes[tree] + ":" + to_string(lengths[tree]) + ", " + tempb + ":" + to_string(scoresP[i].score) + ")";
					treenodes.push_back(addnode);
					cout << addnode << endl;
					tree++;
					vals2.push_back(tempb);
					p = 2 * i;
					q = 2 * i;
					check = 1;
				}
			}
		}
		for (int p = 0; p < 2 * i; p++) {
			for (int q = 0; q < 2 * i; q++) {
				if (scoresP[i].secondb == vals[q] && check == 0) {
					string addnode = "(" + tempa + ":" + to_string(scoresP[i].score) + ", " + treenodes[tree] + ":" + to_string(lengths[tree]) + ")";
					treenodes.push_back(addnode);
					cout << addnode << endl;
					tree++;

					vals2.push_back(tempa);
					p = 2 * i;
					q = 2 * i;
					check = 1;

				}
			}
		}
	
			 if (check == 0 ) {
				string addnode = "(" + tempa + ":" + to_string(scoresP[i].score) + ", " + tempb + ":" + to_string(scoresP[i].score) + ")";
				treenodes.push_back(addnode);
				cout << addnode << endl;
				vals2.push_back(tempb);
				vals2.push_back(tempa);
				connected.push_back(Point2(tempa, tempb));
			
		}
	}

	cout << endl;
	for (int i = 0; i < treenodes.size(); i++) {
		cout << "treenodes stuff: " << treenodes[i] << " " << endl;
	}
	cout << " Final answer:  " << treenodes.back() << ";" << endl;

    return "succeed";
}



// Makes phylogenetic in newick format  
int main()
{
	clock_t init, final;
	init = clock();
    vector<Point> Seqs;
	string line;
    string name;
    string seq;
    
	ifstream ifs("250seq.txt");
	while (getline(ifs, line))
	{
		stringstream ss(line);
		ss >> name >> seq;
		//convert to capital letters
		for (int i = 0; i<seq.size(); i++)
			seq.at(i) = toupper(seq.at(i));
		//add to vector for Points
		Seqs.push_back(Point(name, seq));
	}

    //call algorithm PhyloTreeAlgorithm(GENE_SEQ)
    cout << PhyloTreeAlgorithm(Seqs);

	final = clock() - init;
	cout << "Total Runtime: " << (double)final / ((double)CLOCKS_PER_SEC);
	//for calculations of efficiency
    

    
}
