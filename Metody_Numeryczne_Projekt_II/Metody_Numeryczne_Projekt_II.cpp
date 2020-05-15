// Metody_Numeryczne_Projekt_II.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <ctime>
#include <string>
#include <fstream>
using namespace std;

double* stworzWektorB(int N, int f)
{
	double* wektorB = new double[N];
	for(int i =0; i < N; i++)
	{
		wektorB[i] = sin(i*(f + 1));
	}	

	return wektorB;
}

double* stworzWektorX(int N)
{
	double* wektorX = new double[N];
	for (int i = 0; i < N; i++)
	{
		wektorX[i] = 0;
	}
	return wektorX;
}

double** stworzIWypelnijUkladRownan(int N, int a1, int a2, int a3)
{
	double** ukladRownan = new double*[N];

	for (int i = 0; i < N; i++)
	{
		ukladRownan[i] = new double[N];
	}

	int cols = N;
	int rows = N;

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (j == i) {
				ukladRownan[j][i] = a1;
			}else if((j - 1 == i) || (j + 1 == i))
			{
				ukladRownan[j][i] = a2;				
			}else if((j - 2 == i) || (j + 2 == i))
			{
				ukladRownan[j][i] = a3;				
			}else
			{
				ukladRownan[j][i] = 0;				
			}
		}
	}
	return ukladRownan;
}

void wypiszUkladRownan(double** ukladRownan, int N)
{
	int cols = N;
	int rows = N;

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
				cout << ukladRownan[j][i] << " ";
		}
		cout << endl;
	}
}

void zwolnijUklad(double** ukladRownan, int N)
{
	for(int i = 0; i < N; i++)
	{
		delete[] ukladRownan[i];
	}
	delete[] ukladRownan;
}

void metodaJacobiego(double** ukladRownan, double* wektorB, double* wektorX, int N)
{
	double* wektorXNext = new double[N];

	for(int kolumna = 0; kolumna < N; kolumna++)
	{
		double firstSum = 0;
		double secondSum = 0;
		wektorXNext[kolumna] = wektorB[kolumna];
		for(int wiersz = 0; wiersz < kolumna-1 ; wiersz++)
		{
			firstSum += ukladRownan[kolumna][wiersz] * wektorX[wiersz];
		}
		for(int wiersz = kolumna+1; wiersz < N; wiersz++)
		{
			secondSum += ukladRownan[kolumna][wiersz] * wektorX[wiersz];
		}
		wektorXNext[kolumna] = (wektorXNext[kolumna] - firstSum - secondSum)/(ukladRownan[kolumna][kolumna]);		
	}

	for (int i = 0; i < N; i++)
	{
		wektorX[i] = wektorXNext[i];
	}

}

void metodaGaussaSeidela(double** ukladRownan, double* wektorB, double* wektorX, int N)
{
	for (int kolumna = 0; kolumna < N; kolumna++)
	{
		double firstSum = 0;
		double secondSum = 0;
		//wektorX[kolumna] = wektorB[kolumna];
		for (int wiersz = 1; wiersz < (kolumna - 1); wiersz++)
		{
			firstSum += ukladRownan[kolumna][wiersz] * wektorX[wiersz];
		}
		for (int wiersz = kolumna+1; wiersz < N; wiersz++)
		{
			secondSum += ukladRownan[kolumna][wiersz] * wektorX[wiersz];
			//cout << "ukladRownan[" << i << "][" << j <<"] * wektorX[" << j << "]" << endl;
		}
		wektorX[kolumna] = ((wektorB[kolumna] - firstSum - secondSum) / ukladRownan[kolumna][kolumna]);
	}
}

double* wektorResiduum(double** ukladRownan, double* wektorB, double*wektorX, int N)
{
	double* res = new double[N];
	
	//Mnozenie Macierzy Ax^k
	for(int kolumna = 0; kolumna < N; kolumna++)
	{
		res[kolumna] = 0;
		for(int wiersz = 0; wiersz < N; wiersz++)
		{
			res[kolumna] += ukladRownan[wiersz][kolumna] * wektorX[wiersz];
		}
	}
	//Odejmowanie wektoraB
	for(int i = 0; i < N; i++)
	{
		res[i] -= wektorB[i];
	}

	return res;
}

double norm(double* wektorResiduum, int N)
{
	double wynik = 0;

	for(int i = 0; i < N; i++)
	{
		wynik += wektorResiduum[i] * wektorResiduum[i];
	}
	wynik = sqrt(wynik);
	//cout << "wynik: " << wynik << endl;
	return wynik;
}

double module(double x)
{
	if(x < 0)
	{
		return -x;
	}
	return x;
}

int main()
{
	//indeks
	int indeks[6] = { 1,6,5,6,0,9 };
	//zmienne potrzebne do inicjalizacji
	int c = indeks[4], d = indeks[5], e = indeks[3], f = indeks[2];
	int a1 = 5 + e;
	int a2 = -1, a3 = -1;
//	int N = 9 * 100 + c * 10 + d;
	int N = 10;
	//kryteria i liczniki
	double kryterium = 0.000000001;
	int counterJacobi = 0;
	int counterGauss = 0;
	int counterCJacobi = 0;
	int counterCGauss = 0;
	double timer = clock();
	double timer2 = clock();
	//Stworznie i wypelnienie potrzebnych macierzy do metody Jacobiego
	double** ukladRownan = stworzIWypelnijUkladRownan(N, a1, a2, a3);
	double** ukladRownanC = stworzIWypelnijUkladRownan(N, 3, -1, -1);
	double* wektorB = stworzWektorB(N, f);
	double* wektorXJacobi = stworzWektorX(N);
	double* wektorXGauss = stworzWektorX(N);
	double* wektorXZadCJacobi = stworzWektorX(N);
	double* wektorXZadCGauss = stworzWektorX(N);
//	wypiszUkladRownan(ukladRownan,N);
//
//	cout << "Jacobi\n";
//	cout << "Norma przed\n";
//	cout << norm(wektorResiduum(ukladRownan, wektorB, wektorXJacobi, N), N) << endl;
//	metodaJacobiego(ukladRownan, wektorB, wektorXJacobi, N);
//	cout << "Norma po\n";
//	cout << norm(wektorResiduum(ukladRownan, wektorB, wektorXJacobi, N), N) << endl;
//
//	cout << "Gauss\n";
//	cout << "Norma przed\n";
//	cout << norm(wektorResiduum(ukladRownan, wektorB, wektorXGauss, N), N) << endl;
//	metodaGaussaSeidela(ukladRownan, wektorB, wektorXGauss, N);
//	cout << "Norma po\n";
//	cout << norm(wektorResiduum(ukladRownan, wektorB, wektorXGauss, N), N) << endl;
//

//	//Metoda Jacobiego	
//	timer = clock();
//	double tmp = norm(wektorResiduum(ukladRownan, wektorB, wektorXJacobi, N), N);
//	while (true)
//	{
//		metodaJacobiego(ukladRownan, wektorB, wektorXJacobi, N);
//		double tmp2 = norm(wektorResiduum(ukladRownan, wektorB, wektorXJacobi, N), N);
//		//Porownujemy wyliczone normy i sprawdzamy z kryterium
//		if (module(tmp2 - tmp) < kryterium) {
//			break;
//		}
//		counterJacobi++;
//		tmp = tmp2;
//	}
//	timer2 = clock();
//	cout << "ZAD A Jacobi:\t iterations = " << counterJacobi << " time = " << (timer2 - timer) / CLOCKS_PER_SEC << "s" << endl;
//
//	//Metoda Gaussa-Seidela
//	timer = clock();
//	tmp = norm(wektorResiduum(ukladRownan, wektorB, wektorXGauss, N), N);
//	while (true)
//	{
//		metodaGaussaSeidela(ukladRownan, wektorB, wektorXGauss, N);
//		double tmp2 = norm(wektorResiduum(ukladRownan, wektorB, wektorXGauss, N), N);
//		//Porownujemy wyliczone normy i sprawdzamy z kryterium
//		if (module(tmp2 - tmp) < kryterium) {
//			break;
//		}
//		counterGauss++;
//		tmp = tmp2;
//	}
//	timer2 = clock();
//	cout << "ZAD A Gauss:\t iterations = " << counterGauss << " time = " << (timer2 - timer) / CLOCKS_PER_SEC << "s" << endl;
//
//	//Zadanie C Jacobi
//	timer = clock();
//	tmp = norm(wektorResiduum(ukladRownanC, wektorB, wektorXZadCJacobi, N), N);
//	while (true)
//	{
//		metodaJacobiego(ukladRownanC, wektorB, wektorXZadCJacobi, N);
//		double tmp2 = norm(wektorResiduum(ukladRownanC, wektorB, wektorXZadCJacobi, N), N);
//		//Porownujemy wyliczone normy i sprawdzamy z kryterium
//		if (tmp2 - tmp < kryterium) {
//			break;
//		}
//		//cout << counterCJacobi << ": " << tmp << ", " << tmp2 << endl << endl << endl;
//		counterCJacobi++;
//		tmp = tmp2;
//	}
//	timer2 = clock();
//	cout << "ZAD C Jacobi:\t iterations = " << counterCJacobi << " time = " << (timer2 - timer) / CLOCKS_PER_SEC << "s" << endl;
//
//	//Zadanie C Jacobi
//	timer = clock();
//	tmp = norm(wektorResiduum(ukladRownanC, wektorB, wektorXZadCGauss, N), N);
//	while (true)
//	{
//		metodaGaussaSeidela(ukladRownanC, wektorB, wektorXZadCGauss, N);
//		double tmp2 = norm(wektorResiduum(ukladRownanC, wektorB, wektorXZadCGauss, N), N);
//		//Porownujemy wyliczone normy i sprawdzamy z kryterium
//		if (tmp2 - tmp < kryterium) {
//			break;
//		}
//		//cout << counterCGauss << ": " << tmp << ", " << tmp2 << endl << endl << endl;
//		counterCGauss++;
//		tmp = tmp2;
//	}
//	timer2 = clock();
//	cout << "ZAD C Gauss:\t iterations = " << counterCGauss << " time = " << (timer2 - timer) / CLOCKS_PER_SEC << "s" << endl;
//
//	//Zadanie D
	double** macierzA = stworzIWypelnijUkladRownan(N, 3, -1, -1);
	double** macierzL = stworzIWypelnijUkladRownan(N, 1, 0, 0);
	double** macierzU = stworzIWypelnijUkladRownan(N, 3, -1, -1);
	double* wektorY = stworzWektorX(N);
	double* wektorX = stworzWektorX(N);
	
	for (int k = 0; k < N - 1; k++)
	{
		for (int j = k + 1; j < N; j++)
		{
			macierzL[k][j] = macierzU[k][j] / macierzU[k][k];
			for (int i = k; i < N; i++) {
				macierzU[i][j] = macierzU[i][j] - macierzL[k][j] * macierzU[i][k];
			}
		}
	}

	//cout << "Macierz B" << endl;
	//wypiszUkladRownan(ukladRownan, N);
	cout << "Macierz A" << endl;
	wypiszUkladRownan(macierzA, N);
	cout << "Macierz L" << endl;
	wypiszUkladRownan(macierzL, N);
	cout << "Macierz U" << endl;
	wypiszUkladRownan(macierzU, N);

	//Obliczenie wektora y;
	for(int i = 0; i < N; i++)
	{
		double suma = 0;
		for(int j = 0; j < i; j++)
		{
			suma += macierzL[j][i] * wektorY[j];
		}
		wektorY[i] = (wektorB[i] - suma) / macierzL[i][i];
	}

	cout << "Wektor B" << endl;
	for (int i = 0; i < N; i++)
	{
		printf("%f\n", wektorB[i]);
	}

	cout << "Wektor Y" << endl;
	for (int i = 0; i < N; i++)
	{
		printf("%f\n", wektorY[i]);
	}

	//Obliczenie wektora x;
	for (int i = N - 1 ; i >= 0 ; i--)
	{
		double suma = 0;
		for (int j = N-1 ; j > i ; j--)
		{
			suma += macierzU[j][i] * wektorX[j];
		}
		wektorX[i] = (wektorY[i] - suma) / macierzU[i][i];
	}

	cout << "Wektor X" << endl;
	for(int i = 0; i < N; i++)
	{
		printf("%f\n", wektorX[i]);
	}


//		//Tworzenie wykresow dla roznych danych
////	fstream file;
////
////	//Zapis do pliku.csv
////	file.open("WynikiKoncowe.csv", ios::out | ios::app);
////	file << "Type , Size , Iterations , Time \n";
////
////	const int ROZMIAR = 12;
////	int wielkoscUkladow[ROZMIAR] = { 100,500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000 };
////
////	for (int i = 0; i < ROZMIAR; i++)
////	{
////		double** ukladRownanE = stworzIWypelnijUkladRownan(wielkoscUkladow[i], a1, a2, a3);
////		double* wektorXJacobiE = stworzWektorX(wielkoscUkladow[i]);
////		double* wektorBE = stworzWektorB(wielkoscUkladow[i], f);
////		int counterJacobiE = 0;
////
////		//Jacobi
////		timer = clock();
////		double tmpE = norm(wektorResiduum(ukladRownanE, wektorBE, wektorXJacobiE, wielkoscUkladow[i]), wielkoscUkladow[i]);
////		while (true)
////		{
////			metodaJacobiego(ukladRownanE, wektorBE, wektorXJacobiE, wielkoscUkladow[i]);
////			double tmp2E = norm(wektorResiduum(ukladRownanE, wektorBE, wektorXJacobiE, wielkoscUkladow[i]), wielkoscUkladow[i]);
////			//Porownujemy wyliczone normy i sprawdzamy z kryterium
////			if (module(tmp2E - tmpE) < kryterium) {
////				break;
////			}
////			counterJacobiE++;
////			tmpE = tmp2E;
////		}
////		timer2 = clock();
////		file << "Jacobi" << "," << wielkoscUkladow[i] << "," << counterJacobiE << "," << ((timer2 - timer) / CLOCKS_PER_SEC) << "\n";
////
////		//Zwalnianie Pamieci
////		zwolnijUklad(ukladRownanE, wielkoscUkladow[i]);
////		delete wektorBE;
////		delete wektorXJacobiE;
////	}
////
////	for (int i = 0; i < ROZMIAR; i++)
////	{
////		double** ukladRownanE = stworzIWypelnijUkladRownan(wielkoscUkladow[i], a1, a2, a3);
////		double* wektorXGaussE = stworzWektorX(wielkoscUkladow[i]);
////		double* wektorBE = stworzWektorB(wielkoscUkladow[i], f);
////		int counterGaussE = 0;
////
////		//Gauss-Siedel
////		timer = clock();
////		double tmpE = norm(wektorResiduum(ukladRownanE, wektorBE, wektorXGaussE, wielkoscUkladow[i]), wielkoscUkladow[i]);
////		while (true)
////		{
////			metodaGaussaSeidela(ukladRownanE, wektorBE, wektorXGaussE, wielkoscUkladow[i]);
////			double tmp2E = norm(wektorResiduum(ukladRownanE, wektorBE, wektorXGaussE, wielkoscUkladow[i]), wielkoscUkladow[i]);
////			//Porownujemy wyliczone normy i sprawdzamy z kryterium
////			if (module(tmp2E - tmpE) < kryterium) {
////				break;
////			}
////			counterGaussE++;
////			tmpE = tmp2E;
////		}
////		timer2 = clock();
////		file << "Gauss" << "," << wielkoscUkladow[i] << "," << counterGaussE << "," << ((timer2 - timer) / CLOCKS_PER_SEC) << "\n";
////
////
////		//Zwalnianie Pamieci
////		zwolnijUklad(ukladRownanE, wielkoscUkladow[i]);
////		delete wektorBE;
////		delete wektorXGaussE;
////	}
////	for (int g = 0; g < 6; g++)
////	{
////		//faktoryzacja LU
////		double** macierzAE = stworzIWypelnijUkladRownan(wielkoscUkladow[g], 3, -1, -1);
////		double** macierzLE = stworzIWypelnijUkladRownan(wielkoscUkladow[g], 1, 0, 0);
////		double** macierzUE = stworzIWypelnijUkladRownan(wielkoscUkladow[g], 3, -1, -1);
////
////		timer = clock();
////		for (int i = 1; i < wielkoscUkladow[g] - 1; i++)
////		{
////			//			cout << "test" << endl;
////			for (int j = i + 1; j < wielkoscUkladow[g]; j++)
////			{
////				macierzLE[i][j] = macierzUE[i][j] / macierzUE[j][j];
////				for (int k = j; k < wielkoscUkladow[g]; k++) {
////					macierzUE[i][k] = macierzUE[i][k] - macierzLE[i][k] * macierzUE[j][k];
////				}
////			}
////		}
////		timer2 = clock();
////		file << "LU" << "," << wielkoscUkladow[g] << "," << "x" << "," << ((timer2 - timer) / CLOCKS_PER_SEC) << "\n";
////
////		//Zwalnianie Pamieci
////		zwolnijUklad(macierzAE, wielkoscUkladow[g]);
////		zwolnijUklad(macierzLE, wielkoscUkladow[g]);
////		zwolnijUklad(macierzUE, wielkoscUkladow[g]);
////	}
////	file.close();
//
//	cout << "Koniec" << endl;
//	return 0;
}

