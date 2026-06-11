#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <cstdlib>

using namespace std;
using namespace std::chrono;

typedef complex<double> cd;
const double PI = acos(-1.0);

/**
 * 1. Thuật toán Biến đổi Fourier rời rạc thông thường (Brute-force DFT)
 * Độ phức tạp: O(N^2)
 */
vector<cd> bruteForceDFT(const vector<cd>& a) {
    int n = a.size();
    vector<cd> res(n, 0);
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            // Hệ số xoay W_N^(j*k) = e^(2*pi*i*j*k/N)
            double angle = 2.0 * PI * k * j / n;
            cd w(cos(angle), sin(angle));
            res[k] += a[j] * w;
        }
    }
    return res;
}

int reverseBits(int num, int log2n) {
    int res = 0;
    for (int i = 0; i < log2n; ++i) {
        if ((num & (1 << i)) != 0) {
            res |= (1 << (log2n - 1 - i));
        }
    }
    return res;
}

/**
 * 2. Thuật toán Cooley--Tukey FFT (Radix-2 Iterative)
 * Ghi đè trực tiếp tại chỗ (In-place) sau khi đảo bit
 * Độ phức tạp: O(N log_2 N)
 */
vector<cd> cooleyTukeyFFT(vector<cd> a) {
    int n = a.size();
    int log2n = 0;
    while ((1 << log2n) < n) log2n++;
    
    // Bước 1: Digit-Reversal
    for (int i = 0; i < n; ++i) {
        int rev = reverseBits(i, log2n);
        if (i < rev) {
            swap(a[i], a[rev]);
        }
    }
    
    // Bước 2: Vòng lặp tính toán các khối từ dưới lên (Iterative)
    for (int len = 2; len <= n; len <<= 1) {
        double angle = 2.0 * PI / len;
        cd wlen(cos(angle), sin(angle)); // Hệ số xoay cơ sở của khối
        
        for (int i = 0; i < n; i += len) {
            cd w(1.0, 0.0);
            for (int j = 0; j < len / 2; ++j) {
                cd u = a[i + j];
                cd v = a[i + j + len / 2] * w;
                
                // Cập nhật giá trị tại chỗ bằng phép bướm (Butterfly operation)
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    return a;
}

int main() {
    cout << fixed << setprecision(2);
    // Các kích thước N cần chạy thực nghiệm (phải là lũy thừa của 2)
    vector<int> test_sizes = {256, 1024, 4096, 16384};
    
    cout << "=========================================================\n";
    cout << "   KHAO SAT HIEU NANG: BRUTE-FORCE DFT VS COOLEY-TUKEY  \n";
    cout << "=========================================================\n";
    cout << setw(12) << "Kich thuoc (N)" 
         << setw(22) << "Brute-force DFT (us)" 
         << setw(22) << "Cooley-Tukey (us)" << endl;
    cout << "---------------------------------------------------------\n";

    for (int n : test_sizes) {
        // Sinh dữ liệu số phức ngẫu nhiên cho mảng đầu vào
        vector<cd> a(n);
        for (int i = 0; i < n; ++i) {
            double real_part = rand() % 100;
            double imag_part = rand() % 100;
            a[i] = cd(real_part, imag_part);
        }
        // 1. Đo thời gian chạy của Brute-force DFT
        auto start_dft = high_resolution_clock::now();
        vector<cd> res_dft = bruteForceDFT(a);
        auto stop_dft = high_resolution_clock::now();
        auto duration_dft = duration_cast<microseconds>(stop_dft - start_dft);
        // 2. Đo thời gian chạy của Cooley--Tukey FFT
        auto start_fft = high_resolution_clock::now();
        vector<cd> res_fft = cooleyTukeyFFT(a);
        auto stop_fft = high_resolution_clock::now();
        auto duration_fft = duration_cast<microseconds>(stop_fft - start_fft);
        // In kết quả
        cout << setw(12) << n 
             << setw(22) << duration_dft.count() 
             << setw(22) << duration_fft.count() << endl;
    }
    cout << "=========================================================\n";
    return 0;
}