#include <iostream>
#include "TichVoHuong.h" // Nhúng file bạn vừa tạo vào đây

using namespace std;

int main() {
    int luaChon;
    do {
        cout << "====================================================\n";
        cout << "               MENU DO AN SONG SONG                 \n";
        cout << "====================================================\n";
        cout << "1. Phan A: Tinh tich vo huong vector (OpenMP)\n";
        cout << "2. Phan B: Nhan ma tran bang thuat toan Strassen (OpenMP)\n";
        cout << "0. Thoat chuong trinh\n";
        cout << "----------------------------------------------------\n";
        cout << "Nhap lua chon cua ban: ";
        cin >> luaChon;

        // Bẫy lỗi nếu người dùng gõ chữ thay vì gõ số
        if (cin.fail()) {
            cout << "[LOI] Vui long chi nhap ky tu SO (1, 2, 0) tren Menu!\n\n";
            cin.clear();             // Xóa trạng thái lỗi của cin
            cin.ignore(10000, '\n');   // Xóa ký tự bẩn trong bộ nhớ đệm
            luaChon = -1;            // Gán giá trị tạm để switch-case nhảy vào default
            continue;
        }

        switch (luaChon) {
            case 1:
                chayPhanATichVoHuong(); 
                break;
            case 2:
                chayPhanBStrassen(); 
                break;
            case 0:
                cout << "Tam biet!\n";
                break;
            default:
                cout << "Lua chon khong hop le! Vui long nhap lai.\n\n";
        }
    } while (luaChon != 0);

    return 0;
}