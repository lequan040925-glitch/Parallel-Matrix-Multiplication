#include "Strassen.h"
#include <iostream>
#include <vector>
#include <omp.h>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <cmath>

using namespace std;

// Hàm bổ trợ: Cộng hai ma trận song song
void congMatran(const vector<vector<int>>& A, const vector<vector<int>>& B, vector<vector<int>>& C, int n) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

// Hàm bổ trợ: Trừ hai ma trận song song
void truMatran(const vector<vector<int>>& A, const vector<vector<int>>& B, vector<vector<int>>& C, int n) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}

// Hàm nhân ma trận chuẩn (Thuật toán cơ sở phục vụ cho Base Case và chạy tuần tự kiểm tra)
void nhanMatranChuan(const vector<vector<int>>& A, const vector<vector<int>>& B, vector<vector<int>>& C, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int sum = 0;
            for (int k = 0; k < n; k++) {
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }
}

// Thuật toán đệ quy Strassen kết hợp OpenMP nâng cao
void thuatToanStrassen(const vector<vector<int>>& A, const vector<vector<int>>& B, vector<vector<int>>& C, int n, int depth) {
    // Base case: Nếu kích thước ma trận nhỏ (<= 64), dùng nhân ma trận thường song song hóa 
    // để tránh hao phí tạo luồng đệ quy quá nhỏ (Thread Overhead)
    if (n <= 64) {
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int sum = 0;
                for (int k = 0; k < n; k++) {
                    sum += A[i][k] * B[k][j];
                }
                C[i][j] = sum;
            }
        }
        return;
    }

    int k = n / 2;
    vector<vector<int>> A11(k, vector<int>(k)), A12(k, vector<int>(k)), A21(k, vector<int>(k)), A22(k, vector<int>(k));
    vector<vector<int>> B11(k, vector<int>(k)), B12(k, vector<int>(k)), B21(k, vector<int>(k)), B22(k, vector<int>(k));

    // Bước 1: Chia ma trận gốc thành 4 ma trận con song song
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j + k];
            A21[i][j] = A[i + k][j];
            A22[i][j] = A[i + k][j + k];

            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j + k];
            B21[i][j] = B[i + k][j];
            B22[i][j] = B[i + k][j + k];
        }
    }

    vector<vector<int>> P1(k, vector<int>(k)), P2(k, vector<int>(k)), P3(k, vector<int>(k)), P4(k, vector<int>(k)),
                        P5(k, vector<int>(k)), P6(k, vector<int>(k)), P7(k, vector<int>(k));

    // Bước 2: Tính toán song song 7 tích trung gian P1 -> P7 sử dụng OpenMP Sections.
    // Chỉ kích hoạt song song hóa tầng đệ quy đầu (depth < 2) để tránh bùng nổ số luồng ảo trên Core i5.
    #pragma omp parallel sections if(depth < 2)
    {
        #pragma omp section
        {
            vector<vector<int>> AR(k, vector<int>(k)), BR(k, vector<int>(k));
            congMatran(A11, A22, AR, k);
            congMatran(B11, B22, BR, k);
            thuatToanStrassen(AR, BR, P1, k, depth + 1);
        }
        #pragma omp section
        {
            vector<vector<int>> AR(k, vector<int>(k));
            congMatran(A21, A22, AR, k);
            thuatToanStrassen(AR, B11, P2, k, depth + 1);
        }
        #pragma omp section
        {
            vector<vector<int>> BR(k, vector<int>(k));
            truMatran(B12, B22, BR, k);
            thuatToanStrassen(A11, BR, P3, k, depth + 1);
        }
        #pragma omp section
        {
            vector<vector<int>> BR(k, vector<int>(k));
            truMatran(B21, B11, BR, k);
            thuatToanStrassen(A22, BR, P4, k, depth + 1);
        }
        #pragma omp section
        {
            vector<vector<int>> AR(k, vector<int>(k));
            congMatran(A11, A12, AR, k);
            thuatToanStrassen(AR, B22, P5, k, depth + 1);
        }
        #pragma omp section
        {
            vector<vector<int>> AR(k, vector<int>(k)), BR(k, vector<int>(k));
            truMatran(A21, A11, AR, k);
            congMatran(B11, B12, BR, k);
            thuatToanStrassen(AR, BR, P6, k, depth + 1);
        }
        #pragma omp section
        {
            vector<vector<int>> AR(k, vector<int>(k)), BR(k, vector<int>(k));
            truMatran(A12, A22, AR, k);
            congMatran(B21, B22, BR, k);
            thuatToanStrassen(AR, BR, P7, k, depth + 1);
        }
    }

    // Các ma trận con kết quả
    vector<vector<int>> C11(k, vector<int>(k)), C12(k, vector<int>(k)), C21(k, vector<int>(k)), C22(k, vector<int>(k));

    // Bước 3: Tổng hợp kết quả từ P1->P7 vào ma trận đích song song
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            C11[i][j] = P1[i][j] + P4[i][j] - P5[i][j] + P7[i][j];
            C12[i][j] = P3[i][j] + P5[i][j];
            C21[i][j] = P2[i][j] + P4[i][j];
            C22[i][j] = P1[i][j] - P2[i][j] + P3[i][j] + P6[i][j];
        }
    }

    // Bước 4: Ráp 4 phần C11, C12, C21, C22 vào ma trận lớn C ban đầu
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            C[i][j] = C11[i][j];
            C[i][j + k] = C12[i][j];
            C[i + k][j] = C21[i][j];
            C[i + k][j + k] = C22[i][j];
        }
    }
}

// Hàm giao diện điều khiển thực nghiệm phần B
void chayPhanBStrassen() {
    cout << "====================================================\n";
    cout << "   CHAY THUC NGHIEM PHAN B - MA TRAN STRASSEN (OpenMP)\n";
    cout << "====================================================\n\n";

    srand(time(NULL));
    char tiepTuc = 'y';

    while (tiepTuc == 'y' || tiepTuc == 'Y') {
        int N, num_threads;
        
        cout << "-> Nhap kich thuoc ma tran vuong N (Phai la luy thua cua 2 như 128, 256, 512, 1024): ";
        cin >> N;

        // Kiểm tra điều kiện lũy thừa của 2
        if (N <= 0 || (N & (N - 1)) != 0) {
            cout << "[LOI] Kich thuoc ma tran phai la luy thua cua 2. Vui long thu lai!\n\n";
            continue;
        }

        cout << "-> Nhap so luong Luong (Threads) muon su dung (Goi y: 2, 4, 8): ";
        cin >> num_threads;

        if (num_threads <= 0) {
            cout << "[LOI] So luong luong phai lon hon 0. Vui long thu lai!\n\n";
            continue;
        }

        cout << "\nDang cap phat va khoi tao ma tran ngau nhien cỡ " << N << "x" << N << "...\n";
        vector<vector<int>> A(N, vector<int>(N));
        vector<vector<int>> B(N, vector<int>(N));
        vector<vector<int>> C_seq(N, vector<int>(N, 0)); // Kết quả chạy tuần tự
        vector<vector<int>> C_par(N, vector<int>(N, 0)); // Kết quả chạy song song Strassen

        // Khởi tạo dữ liệu ngẫu nhiên song song
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] = rand() % 5 + 1;
                B[i][j] = rand() % 5 + 1;
            }
        }
// 1. CHẠY TUẦN TỰ (Sử dụng chính thuật toán Strassen nhưng chạy 1 luồng để làm mốc nền)
        cout << "1. Dang chay Strassen Tuan tu (1 Core)... \n";
        omp_set_num_threads(1); // Ép hệ thống chạy 1 luồng duy nhất
        double start_seq = omp_get_wtime();
        thuatToanStrassen(A, B, C_seq, N, 0); // Gọi Strassen tuần tự
        double time_seq = omp_get_wtime() - start_seq;

        // 2. CHẠY SONG SONG STRASSEN (Đa luồng OpenMP)
        cout << "2. Dang chay Strassen Song song (" << num_threads << " Threads)... \n";
        omp_set_num_threads(num_threads); // Bật lại đa luồng
        double start_par = omp_get_wtime();
        thuatToanStrassen(A, B, C_par, N, 0);
        double time_par = omp_get_wtime() - start_par;

        // 3. KIỂM TRA ĐỘ CHÍNH XÁC KẾT QUẢ
        bool chinhXac = true;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (C_seq[i][j] != C_par[i][j]) {
                    chinhXac = false;
                    break;
                }
            }
        }

        if (!chinhXac) {
            cout << "\n[LOI TINH TOAN] Ket qua Tuan tu va Strassen Song song bi LECH NHAU!\n";
        } else {
            double speedup = time_seq / time_par;
            double efficiency = (speedup / num_threads) * 100.0;

            cout << "\n---------------- KET QUA KHAO SAT ----------------\n";
            cout << fixed << setprecision(5);
            cout << "Thoi gian Tuan tu (Ts)      : " << time_seq << " s\n";
            cout << "Thoi gian Song song (Tp)    : " << time_par << " s\n";
            cout << "Speedup (S)                 : " << speedup << " lan\n";
            cout << "Efficiency (E)              : " << efficiency << " %\n";
            cout << "--------------------------------------------------\n";
        }

        cout << "\nBan co muon tiep tuc test voi kích thuoc khac khong? (y/n): ";
        cin >> tiepTuc;
        cout << "\n";
    }
    cout << "=> KET THUC PHAN B.\n";
}