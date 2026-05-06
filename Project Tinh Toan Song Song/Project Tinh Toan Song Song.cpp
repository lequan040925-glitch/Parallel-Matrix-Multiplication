#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <omp.h>

using namespace std;

int main() {
    // 1. KÍCH THƯỚC VECTOR (Thay đổi số này khi làm báo cáo)
    int N = 50000000; // 50 triệu phần tử

    cout << "Kich thuoc Vector: " << N << " phan tu.\n" << endl;

    // 2. KHỞI TẠO VÀ SINH DỮ LIỆU
    vector<long long> A(N);
    vector<long long> B(N);
    srand(time(NULL));

    cout << "Dang khoi tao du lieu..." << endl;
    for (int i = 0; i < N; i++) {
        A[i] = rand() % 10; // Sinh số từ 0 đến 9
        B[i] = rand() % 10;
    }

    // ==========================================
    // 3. TÍNH TUẦN TỰ (Sequential)
    // ==========================================
    long long dot_product_seq = 0;
    double start_seq = omp_get_wtime();

    for (int i = 0; i < N; i++) {
        dot_product_seq += A[i] * B[i];
    }

    double end_seq = omp_get_wtime();
    double time_seq = end_seq - start_seq;
    cout << "[Tuan tu] Ket qua: " << dot_product_seq << " | Thoi gian: " << time_seq << " s" << endl;

    // ==========================================
    // 4. TÍNH SONG SONG (Parallel bằng OpenMP)
    // ==========================================
    long long dot_product_par = 0;
    int num_threads = 4; // Thử thay đổi số này (2, 4, 8) khi làm báo cáo

    omp_set_num_threads(num_threads); // Thiết lập số luồng

    double start_par = omp_get_wtime();

    // OpenMP: Chia vòng lặp for và cộng dồn an toàn vào biến dot_product_par
#pragma omp parallel for reduction(+:dot_product_par)
    for (int i = 0; i < N; i++) {
        dot_product_par += A[i] * B[i];
    }

    double end_par = omp_get_wtime();
    double time_par = end_par - start_par;
    cout << "[Song song (" << num_threads << " luong)] Ket qua: " << dot_product_par << " | Thoi gian: " << time_par << " s" << endl;

    // ==========================================
    // 5. SO SÁNH HIỆU NĂNG
    // ==========================================
    double speedup = time_seq / time_par;
    double efficiency = speedup / num_threads;

    cout << "\n=== KET QUA DANH GIA ===" << endl;
    cout << "Speedup (S): " << speedup << " lan" << endl;
    cout << "Efficiency (E): " << efficiency * 100 << " %" << endl;

    return 0;
}