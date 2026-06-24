#include "Strassen.h"
#include <iostream>
#include <vector>
#include <omp.h>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <cmath>

using namespace std;

// ==========================================
// CÁC HÀM TUẦN TỰ
// ==========================================
void congMatranTuanTu(const vector<vector<int>>& A, const vector<vector<int>>& B, vector<vector<int>>& C, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = A[i][j] + B[i][j];
}

void truMatranTuanTu(const vector<vector<int>>& A, const vector<vector<int>>& B, vector<vector<int>>& C, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = A[i][j] - B[i][j];
}

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

void thuatToanStrassenTuanTu(const vector<vector<int>>& A, const vector<vector<int>>& B, vector<vector<int>>& C, int n)
{
    if (n <= 64)
    {
        nhanMatranChuan(A, B, C, n);
        return;
    }
    int k = n / 2;
    vector<vector<int>> A11(k, vector<int>(k)), A12(k, vector<int>(k)), A21(k, vector<int>(k)), A22(k, vector<int>(k));
    vector<vector<int>> B11(k, vector<int>(k)), B12(k, vector<int>(k)), B21(k, vector<int>(k)), B22(k, vector<int>(k));
    
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < k; j++)
        {
            A11[i][j] = A[i][j];         A12[i][j] = A[i][j + k];
            A21[i][j] = A[i + k][j];     A22[i][j] = A[i + k][j + k];
            B11[i][j] = B[i][j];         B12[i][j] = B[i][j + k];
            B21[i][j] = B[i + k][j];     B22[i][j] = B[i + k][j + k];
        }
    }
    vector<vector<int>> P1(k, vector<int>(k)), P2(k, vector<int>(k)), P3(k, vector<int>(k)), P4(k, vector<int>(k)),
                        P5(k, vector<int>(k)), P6(k, vector<int>(k)), P7(k, vector<int>(k));
    vector<vector<int>> AR(k, vector<int>(k)), BR(k, vector<int>(k));

    congMatranTuanTu(A11, A22, AR, k); congMatranTuanTu(B11, B22, BR, k); thuatToanStrassenTuanTu(AR, BR, P1, k);
    congMatranTuanTu(A21, A22, AR, k); thuatToanStrassenTuanTu(AR, B11, P2, k);
    truMatranTuanTu(B12, B22, BR, k);  thuatToanStrassenTuanTu(A11, BR, P3, k);
    truMatranTuanTu(B21, B11, BR, k);  thuatToanStrassenTuanTu(A22, BR, P4, k);
    congMatranTuanTu(A11, A12, AR, k); thuatToanStrassenTuanTu(AR, B22, P5, k);
    truMatranTuanTu(A21, A11, AR, k);  congMatranTuanTu(B11, B12, BR, k); thuatToanStrassenTuanTu(AR, BR, P6, k);
    truMatranTuanTu(A12, A22, AR, k);  congMatranTuanTu(B21, B22, BR, k); thuatToanStrassenTuanTu(AR, BR, P7, k);

    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < k; j++)
        {
            C[i][j]       = P1[i][j] + P4[i][j] - P5[i][j] + P7[i][j];
            C[i][j + k]   = P3[i][j] + P5[i][j];
            C[i + k][j]   = P2[i][j] + P4[i][j];
            C[i + k][j+k] = P1[i][j] - P2[i][j] + P3[i][j] + P6[i][j];
        }
    }
}

// ==========================================
// CÁC HÀM BỔ TRỢ SONG SONG
// ==========================================
void congMatran(const vector<vector<int>>& A, const vector<vector<int>>& B, vector<vector<int>>& C, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

void truMatran(const vector<vector<int>>& A, const vector<vector<int>>& B, vector<vector<int>>& C, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}

// ==========================================
// THUẬT TOÁN STRASSEN SONG SONG
// ==========================================
void thuatToanStrassen(const vector<vector<int>>& A, const vector<vector<int>>& B, vector<vector<int>>& C, int n, int depth) {
    if (n <= 64) {
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

    // Bước 1: Chia nhỏ ma trận tuần tự
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            A11[i][j] = A[i][j];         A12[i][j] = A[i][j + k];
            A21[i][j] = A[i + k][j];     A22[i][j] = A[i + k][j + k];
            B11[i][j] = B[i][j];         B12[i][j] = B[i][j + k];
            B21[i][j] = B[i + k][j];     B22[i][j] = B[i + k][j + k];
        }
    }

    vector<vector<int>> P1(k, vector<int>(k)), P2(k, vector<int>(k)), P3(k, vector<int>(k)), P4(k, vector<int>(k)),
                        P5(k, vector<int>(k)), P6(k, vector<int>(k)), P7(k, vector<int>(k));

    // Kiểm soát độ sâu đệ quy: Chỉ tạo luồng mới ở hai tầng đầu tiên
    bool spawn_task = (depth < 2);

    #pragma omp task shared(P1) final(!spawn_task) firstprivate(k, depth)
    {
        vector<vector<int>> AR(k, vector<int>(k)), BR(k, vector<int>(k));
        congMatran(A11, A22, AR, k);  congMatran(B11, B22, BR, k);
        thuatToanStrassen(AR, BR, P1, k, depth + 1);
    }
    
    #pragma omp task shared(P2) final(!spawn_task) firstprivate(k, depth)
    {
        vector<vector<int>> AR(k, vector<int>(k));
        congMatran(A21, A22, AR, k);
        thuatToanStrassen(AR, B11, P2, k, depth + 1);
    }
    
    #pragma omp task shared(P3) final(!spawn_task) firstprivate(k, depth)
    {
        vector<vector<int>> BR(k, vector<int>(k));
        truMatran(B12, B22, BR, k);
        thuatToanStrassen(A11, BR, P3, k, depth + 1);
    }
    
    #pragma omp task shared(P4) final(!spawn_task) firstprivate(k, depth)
    {
        vector<vector<int>> BR(k, vector<int>(k));
        truMatran(B21, B11, BR, k);
        thuatToanStrassen(A22, BR, P4, k, depth + 1);
    }
    
    #pragma omp task shared(P5) final(!spawn_task) firstprivate(k, depth)
    {
        vector<vector<int>> AR(k, vector<int>(k));
        congMatran(A11, A12, AR, k);
        thuatToanStrassen(AR, B22, P5, k, depth + 1);
    }
    
    #pragma omp task shared(P6) final(!spawn_task) firstprivate(k, depth)
    {
        vector<vector<int>> AR(k, vector<int>(k)), BR(k, vector<int>(k));
        truMatran(A21, A11, AR, k);   congMatran(B11, B12, BR, k);
        thuatToanStrassen(AR, BR, P6, k, depth + 1);
    }
    
    #pragma omp task shared(P7) final(!spawn_task) firstprivate(k, depth)
    {
        vector<vector<int>> AR(k, vector<int>(k)), BR(k, vector<int>(k));
        truMatran(A12, A22, AR, k);   congMatran(B21, B22, BR, k);
        thuatToanStrassen(AR, BR, P7, k, depth + 1);
    }

    // Đợi rào chắn động xử lý xong toàn bộ các Task nhánh con
    #pragma omp taskwait

    // Bước 3: Tổng hợp kết quả về ma trận C
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            C[i][j]     = P1[i][j] + P4[i][j] - P5[i][j] + P7[i][j];
            C[i][j+k]   = P3[i][j] + P5[i][j];
            C[i+k][j]   = P2[i][j] + P4[i][j];
            C[i+k][j+k] = P1[i][j] - P2[i][j] + P3[i][j] + P6[i][j];
        }
    }
}

// ==========================================
// HÀM ĐIỀU KHIỂN THỰC NGHIỆM PHẦN B
// ==========================================
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
        vector<vector<int>> C_seq(N, vector<int>(N, 0)); 
        vector<vector<int>> C_par(N, vector<int>(N, 0)); 

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] = rand() % 5 + 1;
                B[i][j] = rand() % 5 + 1;
            }
        }

        // 1. CHẠY STRASSEN TUẦN TỰ
        cout << "1. Dang chay Strassen Tuan tu...\n";
        double start_seq = omp_get_wtime();
        thuatToanStrassenTuanTu(A, B, C_seq, N);
        double time_seq = omp_get_wtime() - start_seq;

        // 2. CHẠY SONG SONG STRASSEN
        cout << "2. Dang chay Strassen Song song (" << num_threads << " Threads)... \n";
        omp_set_num_threads(num_threads); 
        double start_par = omp_get_wtime();
        
        #pragma omp parallel
        {
            
            #pragma omp single
            {
                thuatToanStrassen(A, B, C_par, N, 0);
            }
        }
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
