//
// koala_vneighbor.cpp
//
// Usage: vneighbor FileName1 FileName2 min_distance max_distance
//
//        FileName1:   bulk supercell
//        FileName2: defect supercell
//
// e.g.: vneighbor POSCAR1 POSCAR2 1.2 3.0 (for 1st NNs of fcc, a=3.54)
//       vneighbor POSCAR1 POSCAR2 3.0 3.9 (for 2nd NNs of fcc, a=3.54)
//
// The input files must be in POSCAR format. The code then outputs
// the neighbors of the vacancy within the given shell. The output
// file is "vneighbor_min-max".
//
//
// Written by Yifeng Wu (ywu36@ncsu.edu) on 2/14/2017.
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <iomanip>
using namespace std;

// read POSCARs
// POSCAR format must be rigidly correct!!
// need two POSCARs: POSCAR1 for bulk, POSCAR2 for defect
void read_poscar(string bStrArrElem[10], int bNumArrElem[10], int &bNumElem,
		int &bNumAtom, double &bUnivScal, double bArrLatVec[3][3],
		double bArrAtomFrac[500][3], string vStrArrElem[10],
		int vNumArrElem[10], int &vNumElem, int &vNumAtom, double &vUnivScal,
		double vArrLatVec[3][3], double vArrAtomFrac[500][3], char *FileName1,
		char *FileName2) {
	ifstream file1, file2;
	istringstream str;
	string line, ele;
	double a, b, c;
	int i, temp1, temp2;
	i = temp1 = temp2 = 0;

	file1.open(FileName1);
	if (!file1)
		cout << "File " << FileName1 << " NOT found!" << endl;
	else {
		getline(file1, line);
		getline(file1, line);
		str.str(line);
		if (!(str >> bUnivScal)) {
			cout << FileName1 << ": Scaling number NOT found!" << endl;
			exit(1);
		}
		str.clear();
		for (; i < 3; i++) {
			getline(file1, line);
			str.str(line);
			if (str >> a >> b >> c) {
				str.clear();
				bArrLatVec[i][0] = a;
				bArrLatVec[i][1] = b;
				bArrLatVec[i][2] = c;
			} else {
				cout << FileName1 << ": Lattice parameters NOT found!" << endl;
				exit(1);
			}
			str.clear();
		}
		i = 0;
		getline(file1, line);
		str.str(line);
		while (str >> bStrArrElem[i]) {
			bNumElem++;
			i++;
		}
		i = 0;
		str.clear();
		getline(file1, line);
		str.str(line);
		while (str >> bNumArrElem[i]) {
			temp1++;
			i++;
		}
		if (temp1 != bNumElem) {
			cout << FileName1 << ": # elements Not correct!" << endl;
			exit(1);
		}
		i = 0;
		str.clear();

		while (bNumArrElem[i] && i < bNumElem) {
			bNumAtom += bNumArrElem[i];
			i++;
		}
		getline(file1, line);
		for (i = 0; i < bNumAtom; i++) {
			getline(file1, line);
			str.str(line);
			if (str >> a >> b >> c) {
				str.clear();
				bArrAtomFrac[i][0] = a;
				bArrAtomFrac[i][1] = b;
				bArrAtomFrac[i][2] = c;
				temp2++;
			} else {
				cout << FileName1 << ": Atomic coordinates not found!" << endl;
				exit(1);
			}
			str.clear();
		}
		if (temp2 != bNumAtom) {
			cout << FileName1 << ": # atoms do NOT match!" << endl;
			exit(1);
		}
		i = temp1 = temp2 = 0;
	}

	file2.open(FileName2);
	if (!file2)
		cout << "File " << FileName2 << " NOT found!" << endl;
	else {
		getline(file2, line);
		getline(file2, line);
		str.str(line);
		if (!(str >> vUnivScal)) {
			cout << FileName2 << ": Scaling number NOT found!" << endl;
			exit(1);
		}
		str.clear();
		for (; i < 3; i++) {
			getline(file2, line);
			str.str(line);
			if (str >> a >> b >> c) {
				str.clear();
				vArrLatVec[i][0] = a;
				vArrLatVec[i][1] = b;
				vArrLatVec[i][2] = c;
			} else {
				cout << FileName2 << ": Lattice parameters NOT found!" << endl;
				exit(1);
			}
			str.clear();
		}
		i = 0;
		getline(file2, line);
		str.str(line);
		while (str >> vStrArrElem[i]) {
			vNumElem++;
			i++;
		}
		i = 0;
		str.clear();
		getline(file2, line);
		str.str(line);
		while (str >> vNumArrElem[i]) {
			temp1++;
			i++;
		}
		if (temp1 != vNumElem) {
			cout << FileName2 << ": # elements Not correct!" << endl;
			exit(1);
		}
		i = 0;
		str.clear();

		while (vNumArrElem[i] && i < vNumElem) {
			vNumAtom += vNumArrElem[i];
			i++;
		}
		getline(file2, line);
		for (i = 0; i < vNumAtom; i++) {
			getline(file2, line);
			str.str(line);
			if (str >> a >> b >> c) {
				str.clear();
				vArrAtomFrac[i][0] = a;
				vArrAtomFrac[i][1] = b;
				vArrAtomFrac[i][2] = c;
				temp2++;
			} else {
				cout << FileName2 << ": Atomic coordinates NOT found!" << endl;
				exit(1);
			}
			str.clear();
		}
		if (temp2 != vNumAtom) {
			cout << FileName2 << ": # atoms NOT correct!" << endl;
			exit(1);
		}
	}

	if (bNumElem != vNumElem) {
		cout << "# elements in " << FileName1 << " and " << FileName2
				<< " do NOT match!" << endl;
		exit(1);
	}
}

// change fractional coordinates to Cartesian coordinates
void frac_to_cart(int bNumAtom, double bUnivScal, double bArrLatVec[3][3],
		double bArrAtomFrac[500][3], double bArrAtomCart[500][3], int vNumAtom,
		double vUnivScal, double vArrLatVec[3][3], double vArrAtomFrac[500][3],
		double vArrAtomCart[500][3]) {
	for (int i = 0; i < 3; i++) {
		bArrLatVec[i][0] = bArrLatVec[i][0] * bUnivScal;
		bArrLatVec[i][1] = bArrLatVec[i][1] * bUnivScal;
		bArrLatVec[i][2] = bArrLatVec[i][2] * bUnivScal;
	}

	for (int i = 0; i < bNumAtom; i++) {
		double tempx = bArrAtomFrac[i][0];
		double tempy = bArrAtomFrac[i][1];
		double tempz = bArrAtomFrac[i][2];
		bArrAtomCart[i][0] = bArrLatVec[0][0] * tempx + bArrLatVec[1][0] * tempy
				+ bArrLatVec[2][0] * tempz;
		bArrAtomCart[i][1] = bArrLatVec[0][1] * tempx + bArrLatVec[1][1] * tempy
				+ bArrLatVec[2][1] * tempz;
		bArrAtomCart[i][2] = bArrLatVec[0][2] * tempx + bArrLatVec[1][2] * tempy
				+ bArrLatVec[2][2] * tempz;
	}

	for (int i = 0; i < 3; i++) {
		vArrLatVec[i][0] = vArrLatVec[i][0] * vUnivScal;
		vArrLatVec[i][1] = vArrLatVec[i][1] * vUnivScal;
		vArrLatVec[i][2] = vArrLatVec[i][2] * vUnivScal;
	}

	for (int i = 0; i < vNumAtom; i++) {
		double tempx = vArrAtomFrac[i][0];
		double tempy = vArrAtomFrac[i][1];
		double tempz = vArrAtomFrac[i][2];
		vArrAtomCart[i][0] = vArrLatVec[0][0] * tempx + vArrLatVec[1][0] * tempy
				+ vArrLatVec[2][0] * tempz;
		vArrAtomCart[i][1] = vArrLatVec[0][1] * tempx + vArrLatVec[1][1] * tempy
				+ vArrLatVec[2][1] * tempz;
		vArrAtomCart[i][2] = vArrLatVec[0][2] * tempx + vArrLatVec[1][2] * tempy
				+ vArrLatVec[2][2] * tempz;
	}
}

// calculate bonding length by periodic image correction
double atom_length(double *a, double *b, double LatPara[3][3]) {
	double MinVal[3] = { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
	double temp1[3] = { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
	double temp2[3] = { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
	double Val1, Val2, Val3, Val4;
	for (int i = 0; i < 3; i++) {
		temp1[0] = temp1[0] - LatPara[i][0];
		temp1[1] = temp1[1] - LatPara[i][1];
		temp1[2] = temp1[2] - LatPara[i][2];
		temp2[0] = temp2[0] + LatPara[i][0];
		temp2[1] = temp2[1] + LatPara[i][1];
		temp2[2] = temp2[2] + LatPara[i][2];
		Val1 = pow(MinVal[0], 2) + pow(MinVal[1], 2) + pow(MinVal[2], 2);
		Val2 = pow(temp1[0], 2) + pow(temp1[1], 2) + pow(temp1[2], 2);
		Val3 = pow(temp2[0], 2) + pow(temp2[1], 2) + pow(temp2[2], 2);
		Val4 = min(Val1, min(Val2, Val3));
		if (Val4 == Val1) {
			temp1[0] = MinVal[0];
			temp1[1] = MinVal[1];
			temp1[2] = MinVal[2];
			temp2[0] = MinVal[0];
			temp2[1] = MinVal[1];
			temp2[2] = MinVal[2];
		} else if (Val4 == Val2) {
			temp2[0] = temp1[0];
			temp2[1] = temp1[1];
			temp2[2] = temp1[2];
			MinVal[0] = temp1[0];
			MinVal[1] = temp1[1];
			MinVal[2] = temp1[2];
		} else {
			temp1[0] = temp2[0];
			temp1[1] = temp2[1];
			temp1[2] = temp2[2];
			MinVal[0] = temp2[0];
			MinVal[1] = temp2[1];
			MinVal[2] = temp2[2];
		}
	}
	return sqrt(Val4);
}

// sum the array of integers
int sum_array(int *arr, int n = 1) {
	int i, sum;
	i = sum = 0;
	while (arr[i] && i < n) {
		sum += arr[i];
		i++;
	}
	return sum;
}

// calculate atomic index and species of the vacancy in bulk
void find_vacancy(int bNumAtom, double bArrAtomFrac[500][3],
		double bArrAtomCart[500][3], double bArrLatVec[3][3], int vNumAtom,
		double vUnivScal, double vArrLatVec[3][3], double vArrAtomCart[500][3],
		int &vAtomIndex, double bArrAtomLen[500], double vArrAtomLen[500],
		int bNumElem, int bNumArrElem[10], string bStrArrElem[10],
		string &vStrElem) {
	int i, j;
	double vAtomLen, temp[3];
	for (i = 0; i < bNumAtom; i++) {
		temp[0] = bArrAtomFrac[i][0] * vArrLatVec[0][0]
				+ bArrAtomFrac[i][1] * vArrLatVec[1][0]
				+ bArrAtomFrac[i][2] * vArrLatVec[2][0];
		temp[1] = bArrAtomFrac[i][0] * vArrLatVec[0][1]
				+ bArrAtomFrac[i][1] * vArrLatVec[1][1]
				+ bArrAtomFrac[i][2] * vArrLatVec[2][1];
		temp[2] = bArrAtomFrac[i][0] * vArrLatVec[0][2]
				+ bArrAtomFrac[i][1] * vArrLatVec[1][2]
				+ bArrAtomFrac[i][2] * vArrLatVec[2][2];
		for (j = 0; j < vNumAtom; j++) {
			vAtomLen = atom_length(vArrAtomCart[j], temp, vArrLatVec);
			if (vAtomLen < 0.1 * vUnivScal)
				break;
		}
		if (j != vNumAtom)
			continue;
		vAtomIndex = i;
		for (int k = 0; k < bNumElem; k++) {
			if (i >= sum_array(bNumArrElem, k)
					&& i < sum_array(bNumArrElem, k + 1)) {
				vStrElem = bStrArrElem[k];
				break;
			}
		}
		for (int k = 0; k < bNumAtom; k++)
			bArrAtomLen[k] = atom_length(bArrAtomCart[k],
					bArrAtomCart[vAtomIndex], bArrLatVec);
		for (int k = 0; k < vNumAtom; k++)
			vArrAtomLen[k] = atom_length(vArrAtomCart[k], temp, vArrLatVec);
	}
}

// calculate atomic indices, species, lengths of the shell
void nn_calculator(int bNumAtom, string bStrArrElem[10], int bNumArrElem[10],
		int bNumElem, double bArrAtomLen[500], string nnStrArrElem[50],
		int nnNumArrElem[10], int bnnArrAtomIndex[50], double bnnArrAtomLen[50],
		int vNumAtom, double vArrAtomLen[500], int vnnArrAtomIndex[50],
		double vnnArrAtomLen[50], int &nnNumAtom, double dMin, double dMax) {
	int Num = 0;
	for (int i = 0; i < bNumAtom; i++) {
		if (bArrAtomLen[i] >= dMin && bArrAtomLen[i] <= dMax) {
			string temp;
			for (int j = 0; j < bNumElem; j++) {
				if (i >= sum_array(bNumArrElem, j)
						&& i < sum_array(bNumArrElem, j + 1)) {
					temp = bStrArrElem[j];
					break;
				}
			}
			nnStrArrElem[Num] = temp;
			bnnArrAtomLen[Num] = bArrAtomLen[i];
			bnnArrAtomIndex[Num] = i;
			Num++;
		}
	}
	Num = 0;

	for (int i = 0; i < vNumAtom; i++) {
		if (vArrAtomLen[i] >= dMin && vArrAtomLen[i] <= dMax) {
			vnnArrAtomLen[Num] = vArrAtomLen[i];
			vnnArrAtomIndex[Num] = i;
			Num++;
		}
	}
	nnNumAtom = Num;

	for (int i = 0; i < bNumElem; i++) {
		for (int j = 0; j < nnNumAtom; j++) {
			if (nnStrArrElem[j] == bStrArrElem[i])
				nnNumArrElem[i]++;
		}
	}
}

int main(int argc, char* argv[]) {
	if (argc != 5) {
		cout << endl << "Usage: vneighbor FileName1 FileName2 r_min r_max"
				<< endl << "       FileName1:   bulk supercell" << endl
				<< "       FileName2: defect supercell" << endl << endl;
		cout
				<< " e.g., vneighbor POSCAR1 POSCAR2 1.0 2.7 (1st NNs of fcc, a=3.5)"
				<< endl;
		cout
				<< "       vneighbor POSCAR1 POSCAR2 2.7 3.6 (2nd NNs of fcc, a=3.5)"
				<< endl;
		exit(1);
	}

	/// data of bulk
	string bStrArrElem[10];      // atomic species
	int bNumArrElem[10];         // # of atoms of each element
	int bNumElem = 0;            // # of elements
	int bNumAtom = 0;            // # of atoms
	double bUnivScal = 0;        // universal scaling parameter
	double bArrLatVec[3][3];     // lattice vectors
	double bArrAtomFrac[500][3]; // fractional atomic coordinates
	double bArrAtomCart[500][3]; // Cartesian atomic coordinates

	/// data of defect
	string vStrArrElem[10];      // atomic species
	int vNumArrElem[10];         // # of atoms of each element
	int vNumElem = 0;            // # of elements
	int vNumAtom = 0;            // # of atoms
	double vUnivScal = 0;        // universal scaling parameter
	double vArrLatVec[3][3];     // lattice vectors
	double vArrAtomFrac[500][3]; // fractional atomic coordinates
	double vArrAtomCart[500][3]; // Cartesian atomic coordinates

	/// data of vacancy
	string vStrElem;             // atomic species of vacancy in bulk
	int vAtomIndex;              // atomic index of vacancy in defect
	int nnNumAtom;               // # of atoms in the shell
	string nnStrArrElem[50];     // atomic species in the shell
	int nnNumArrElem[10];        // # of atoms of each element in the shell
	for (int i = 0; i < 10; i++)
		nnNumArrElem[i] = 0;
	double bArrAtomLen[500];     // bonding lengths of vacancy in bulk
	int bnnArrAtomIndex[50];     // atomic indices in the shell of bulk
	double bnnArrAtomLen[50];    // bonding lengths in the shell of bulk
	double vArrAtomLen[500];   // bonding lengths of vacancy in defect
	int vnnArrAtomIndex[50];     // atomic indices in the shell of defect
	double vnnArrAtomLen[50];    // bonding lengths in the shell of defect

	char *FileName1, *FileName2; // names of input files
	FileName1 = argv[1];
	FileName2 = argv[2];
	read_poscar(bStrArrElem, bNumArrElem, bNumElem, bNumAtom, bUnivScal,
			bArrLatVec, bArrAtomFrac, vStrArrElem, vNumArrElem, vNumElem,
			vNumAtom, vUnivScal, vArrLatVec, vArrAtomFrac, FileName1,
			FileName2);

	frac_to_cart(bNumAtom, bUnivScal, bArrLatVec, bArrAtomFrac, bArrAtomCart,
			vNumAtom, vUnivScal, vArrLatVec, vArrAtomFrac, vArrAtomCart);

	find_vacancy(bNumAtom, bArrAtomFrac, bArrAtomCart, bArrLatVec, vNumAtom,
			vUnivScal, vArrLatVec, vArrAtomCart, vAtomIndex, bArrAtomLen,
			vArrAtomLen, bNumElem, bNumArrElem, bStrArrElem, vStrElem);

	double dMin, dMax; // min and max distances of the shell
	dMin = atof(argv[3]);
	dMax = atof(argv[4]);
	nn_calculator(bNumAtom, bStrArrElem, bNumArrElem, bNumElem, bArrAtomLen,
			nnStrArrElem, nnNumArrElem, bnnArrAtomIndex, bnnArrAtomLen,
			vNumAtom, vArrAtomLen, vnnArrAtomIndex, vnnArrAtomLen, nnNumAtom,
			dMin, dMax);

	/// output
	ofstream os;
	char output[] = "vneighbor_";
	strcat(output, argv[3]);
	strcat(output, "-");
	strcat(output, argv[4]);
	os.open(output);

	os.setf(ios_base::fixed, ios_base::floatfield);
	os.precision(5);
	os << "Atomic index and species of the vacancy in bulk:" << endl;
	os << setw(4) << vAtomIndex + 1 << setw(4) << vStrElem << endl;
	os << "Number of bondings:" << endl;
	for (int i = 0; i < bNumElem; i++) {
		os << setw(4) << bStrArrElem[i];
	}
	os << endl;
	for (int i = 0; i < bNumElem; i++) {
		os << setw(4) << nnNumArrElem[i];
	}
	os << endl;

	os << "Atomic index of neighbors in the shell:" << endl;
	os << FileName1 << ":";
	for (int i = 0; i < nnNumAtom; i++) {
		os << setw(4) << bnnArrAtomIndex[i] + 1;
	}
	os << endl;
	os << FileName2 << ":";
	for (int i = 0; i < nnNumAtom; i++) {
		os << setw(4) << vnnArrAtomIndex[i] + 1;
	}
	os << endl;

	os << "and corresponding lengths:" << endl;
	os << FileName1 << ":";
	for (int i = 0; i < nnNumAtom; i++) {
		os << setw(8) << bnnArrAtomLen[i];
	}
	os << endl;
	os << FileName2 << ":";
	for (int i = 0; i < nnNumAtom; i++) {
		os << setw(8) << vnnArrAtomLen[i];
	}
	os << endl;
	os.close();

	exit(0);
}
