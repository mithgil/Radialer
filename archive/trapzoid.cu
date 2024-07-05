#include <iostream>

// CUDA kernel for Trapezoidal rule integration
__global__ void trapezoidalIntegration(float *x, float *y, float *z, int n, float *result) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < n) {
        result[tid] = 0.5 * (z[tid] + z[tid + 1]) * (x[tid + 1] - x[tid]) * (y[tid + 1] - y[tid]);
    }
}

int main() {
    const int n = 1000; // Number of intervals
    const int size = n + 1; // Number of points

    // Host arrays
    float *x = new float[size];
    float *y = new float[size];
    float *z = new float[size];

    // Initialize x, y, and z arrays on the host

    // Device arrays
    float *d_x, *d_y, *d_z, *d_result;
    cudaMalloc(&d_x, size * sizeof(float));
    cudaMalloc(&d_y, size * sizeof(float));
    cudaMalloc(&d_z, size * sizeof(float));
    cudaMalloc(&d_result, n * sizeof(float));

    // Copy input data from host to device
    cudaMemcpy(d_x, x, size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, y, size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_z, z, size * sizeof(float), cudaMemcpyHostToDevice);

    // Configure grid and block dimensions
    int block_size = 256;
    int num_blocks = (n + block_size - 1) / block_size;

    // Launch CUDA kernel
    trapezoidalIntegration<<<num_blocks, block_size>>>(d_x, d_y, d_z, n, d_result);

    // Copy result from device to host
    float *result = new float[n];
    cudaMemcpy(result, d_result, n * sizeof(float), cudaMemcpyDeviceToHost);

    // Calculate the total integral result
    float total_integral = 0.0;
    for (int i = 0; i < n; ++i) {
        total_integral += result[i];
    }

    // Clean up
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] result;
    cudaFree(d_x);
    cudaFree(d_y);
    cudaFree(d_z);
    cudaFree(d_result);

    std::cout << "Total Integral: " << total_integral << std::endl;

    return 0;
}
