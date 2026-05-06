#include "TichVoHuong.h"
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <omp.h>
#include <iomanip>

using namespace std;

void chayPhanATichVoHuong() {
    cout << "====================================================\n";
    cout << "    CHAY THUC NGHIEM PHAN A - NHAP TU BAN PHIM\n";
    cout << "====================================================\n\n";

    // Khởi tạo mầm số ngẫu nhiên một lần duy nhất
    srand(time(NULL));

    char continue_test = 'y';

    // Vòng lặp while để cho phép người dùng test nhiều lần liên tiếp
    while (continue_test == 'y' || continue_test == 'Y') {

        long long N;
        int num_threads;

        // 1. NGƯỜI DÙNG NHẬP KÍCH THƯỚC N
        cout << "-> Nhap kich thuoc Vector N (Goi y: 10000000, 50000000...): ";
        cin >> N;

        if (N <= 0) {
            cout << "[LOI] Kich thuoc phai lon hon 0. Vui long thu lai!\n\n";
            continue;
        }

        // 2. NGƯỜI DÙNG NHẬP SỐ LUỒNG
        cout << "-> Nhap so luong Luong (Threads) muon su dung (Goi y: 2, 4, 8...): ";
        cin >> num_threads;

        if (num_threads <= 0) {
            cout << "[LOI] So luong phai lon hon 0. Vui long thu lai!\n\n";
            continue;
        }

        // ==========================================
        cout << "\nDang cap phat va khoi tao du lieu (" << N << " phan tu)...\n";
        vector<int> A(N);
        vector<int> B(N);
        for (long long j = 0; j < N; j++) {
            A[j] = rand() % 10 + 1;
            B[j] = rand() % 10 + 1;
        }

        // ==========================================
        cout << "1. Dang chay Tuan Tu (1 Core)...\n";
        long long dot_product_seq = 0;
        double start_seq = omp_get_wtime();

        for (long long j = 0; j < N; j++) {
            dot_product_seq += (long long)A[j] * B[j];
        }
        double time_seq = omp_get_wtime() - start_seq;

        // ==========================================
        cout << "2. Dang chay Song Song (" << num_threads << " Threads)...\n";
        omp_set_num_threads(num_threads);
        long long dot_product_par = 0;
        double start_par = omp_get_wtime();

#pragma omp parallel for reduction(+:dot_product_par)
        for (long long j = 0; j < N; j++) {
            dot_product_par += (long long)A[j] * B[j];
        }
        double time_par = omp_get_wtime() - start_par;

        // ==========================================
        // Kiểm tra đúng đắn
        if (dot_product_seq != dot_product_par) {
            cout << "\n[LOI TINH TOAN] Ket qua Tuan tu va Song song KHONG KHOP!\n";
            cout << "Tuan tu ra: " << dot_product_seq << ", Song song ra: " << dot_product_par << "\n";
        }
        else {
            // TÍNH TOÁN HIỆU NĂNG VÀ IN KẾT QUẢ
            double speedup = time_seq / time_par;
            double efficiency = (speedup / num_threads) * 100.0;

            cout << "\n---------------- KET QUA ----------------\n";
            cout << fixed << setprecision(5);
            cout << "Thoi gian Tuan tu (Ts) : " << time_seq << " s\n";
            cout << "Thoi gian Song song(Tp): " << time_par << " s\n";
            cout << "Speedup (S)            : " << speedup << " lan\n";
            cout << "Efficiency (E)         : " << efficiency << " %\n";
            cout << "-----------------------------------------\n";
        }

        // Giải phóng bộ nhớ trước khi test vòng lặp mới
        A.clear(); A.shrink_to_fit();
        B.clear(); B.shrink_to_fit();

        // Hỏi người dùng có muốn tiếp tục test không
        cout << "\nBan co muon tiep tuc test voi so khac khong? (y/n): ";
        cin >> continue_test;
        cout << "\n";
    }

    cout << "=> KET THUC PHAN A.\n";
}