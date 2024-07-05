#include <iostream>
#include <complex>
#include <vector>
#include <thread>
#include <chrono>

// Define your 4D array type for convenience
using Complex4DArray = std::complex<double>****;

// Function to process a segment of the 4D array
void processSegment(Complex4DArray array, int dim1Size, int dim2Size, int dim3Size, int dim4Size, int start, int end) {
    for (int i = start; i < end; ++i) {
        // Perform your computations on array[i]
        // Example: Summing the real parts of all elements
        double sum = 0.0;
        for (int j = 0; j < dim2Size; ++j) {
            for (int k = 0; k < dim3Size; ++k) {
                for (int l = 0; l < dim4Size; ++l) {
                    sum += array[i][j][k][l].real();
                }
            }
        }
        // Example: Assigning a value to array[i][j][k][l]
        array[i][0][0][0] = std::complex<double>(1.0, 2.0); // Example assignment
    }
}

void process4DArraySingleThreaded(Complex4DArray array, int dim1Size, int dim2Size, int dim3Size, int dim4Size) {
    // Single-threaded processing
    processSegment(array, dim1Size, dim2Size, dim3Size, dim4Size, 0, dim1Size);
}

void process4DArrayMultithreaded(Complex4DArray array, int dim1Size, int dim2Size, int dim3Size, int dim4Size, int numThreads) {
    std::vector<std::thread> threads;
    int chunkSize = dim1Size / numThreads;
    
    // Create threads and divide work
    for (int i = 0; i < numThreads; ++i) {
        int start = i * chunkSize;
        int end = (i == numThreads - 1) ? dim1Size : start + chunkSize;
        threads.emplace_back(processSegment, array, dim1Size, dim2Size, dim3Size, dim4Size, start, end);
    }

    // Join threads
    for (auto& thread : threads) {
        thread.join();
    }
}

int main() {
    // Example dimensions
    int dim1Size = 1000, dim2Size = 100, dim3Size = 100, dim4Size = 100;

    // Example 4D array initialization
    Complex4DArray array = new std::complex<double>***[dim1Size];
    for (int i = 0; i < dim1Size; ++i) {
        array[i] = new std::complex<double>**[dim2Size];
        for (int j = 0; j < dim2Size; ++j) {
            array[i][j] = new std::complex<double>*[dim3Size];
            for (int k = 0; k < dim3Size; ++k) {
                array[i][j][k] = new std::complex<double>[dim4Size];
                for (int l = 0; l < dim4Size; ++l) {
                    array[i][j][k][l] = std::complex<double>(i + j + k + l, i - j - k + l); // Example values
                }
            }
        }
    }

    // Benchmark single-threaded processing
    auto startSingle = std::chrono::high_resolution_clock::now();
    process4DArraySingleThreaded(array, dim1Size, dim2Size, dim3Size, dim4Size);
    auto endSingle = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedSingle = endSingle - startSingle;
    std::cout << "Single-threaded processing time: " << elapsedSingle.count() << " seconds\n";

    // Benchmark multithreaded processing
    int numThreads = 4; // Number of threads to use
    auto startMulti = std::chrono::high_resolution_clock::now();
    process4DArrayMultithreaded(array, dim1Size, dim2Size, dim3Size, dim4Size, numThreads);
    auto endMulti = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedMulti = endMulti - startMulti;
    std::cout << "Multithreaded processing time (" << numThreads << " threads): " << elapsedMulti.count() << " seconds\n";

    // Clean up allocated memory
    for (int i = 0; i < dim1Size; ++i) {
        for (int j = 0; j < dim2Size; ++j) {
            for (int k = 0; k < dim3Size; ++k) {
                delete[] array[i][j][k];
            }
            delete[] array[i][j];
        }
        delete[] array[i];
    }
    delete[] array;

    return 0;
}

