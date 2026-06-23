

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
