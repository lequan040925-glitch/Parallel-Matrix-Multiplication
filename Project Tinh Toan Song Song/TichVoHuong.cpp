#include "TichVoHuong.h" // Nhúng header của chính nó
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <omp.h>
#include <iomanip>

using namespace std;

void chayPhanATichVoHuong() {
    long long N = 100000000;
    int num_threads = 4;

    cout << "====================================================\n";
    cout << "   BAI TAP CHUNG: TICH VO HUONG 2 VECTOR (OPENMP)\n";
    cout << "====================================================\n";
    cout << "Kich thuoc Vector (N): " << N << " phan tu.\n";
    cout << "So luong luong (Threads) su dung: " << num_threads << "\n";
    cout << "----------------------------------------------------\n";

    cout << "1. Dang cap phat va khoi tao du lieu ngau nhien...\n";
    vector<int> A(N);
    vector<int> B(N);
    srand(time(NULL));
    for (long long i = 0; i < N; i++) {
        A[i] = rand() % 10 + 1;
        B[i] = rand() % 10 + 1;
    }
    cout << "   -> Da khoi tao xong!\n\n";

    cout << "2. Bat dau tinh toan TUAN TU (1 Core)...\n";
    long long dot_product_seq = 0;
    double start_seq = omp_get_wtime();
    for (long long i = 0; i < N; i++) {
        dot_product_seq += (long long)A[i] * B[i];
    }
    double end_seq = omp_get_wtime();
    double time_seq = end_seq - start_seq;
    cout << "   -> Ket qua: " << dot_product_seq << "\n";
    cout << fixed << setprecision(5);
    cout << "   -> Thoi gian (Ts): " << time_seq << " giay.\n\n";

    cout << "3. Bat dau tinh toan SONG SONG (" << num_threads << " Threads)...\n";
    long long dot_product_par = 0;
    omp_set_num_threads(num_threads);
    double start_par = omp_get_wtime();

#pragma omp parallel for reduction(+:dot_product_par)
    for (long long i = 0; i < N; i++) {
        dot_product_par += (long long)A[i] * B[i];
    }

    double end_par = omp_get_wtime();
    double time_par = end_par - start_par;
    cout << "   -> Ket qua: " << dot_product_par << "\n";
    cout << "   -> Thoi gian (Tp): " << time_par << " giay.\n\n";

    if (dot_product_seq == dot_product_par) {
        cout << "=> [THANH CONG] Ket qua Tuan tu va Song song khop nhau!\n\n";
    }
    else {
        cout << "=> [LOI] Ket qua khong khop! Hay kiem tra lai Race Condition!\n\n";
    }

    double speedup = time_seq / time_par;
    double efficiency = speedup / num_threads;

    cout << "====================================================\n";
    cout << "                 BANG DANH GIA                      \n";
    cout << "====================================================\n";
    cout << " - Speedup (He so tang toc) : " << speedup << " lan\n";
    cout << " - Efficiency (Hieu suat)   : " << efficiency * 100 << " %\n";
    cout << "====================================================\n";
}