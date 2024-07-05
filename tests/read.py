import numpy as np

def parse_complex(s):
    # Assuming the format "real+imagi" or "real-imagi"
    if '+' in s:
        real, imag = s.split('+')
        imag = imag.rstrip('i')
    elif '-' in s:
        real, imag = s.split('-')
        imag = '-' + imag.rstrip('i')
    else:
        raise ValueError("Invalid complex number format")
    return complex(float(real), float(imag))

def read_4d_complex_array_from_file(filename, dim1, dim2, dim3, dim4):
    with open(filename, 'r') as file:
        lines = file.readlines()

    array = np.zeros((dim1, dim2, dim3, dim4), dtype=complex)
    idx = 0 # initialization of linear index
    
    for line in lines:
        stripped_line = line.strip()
        if stripped_line:
            # Split the line into individual complex numbers
            complex_numbers = stripped_line.split()
            for complex_str in complex_numbers:
                # Calculate 4D indices from the linear index
                i = idx // (dim2 * dim3 * dim4)
                j = (idx // (dim3 * dim4)) % dim2
                k = (idx // dim4) % dim3
                l = idx % dim4
                # Assign the complex number to the appropriate position in the array
                array[i, j, k, l] = parse_complex(complex_str)
                idx += 1  # Increment the linear index

    return array

# Example usage:
filename = '../build/radial_array.txt'
dim1, dim2, dim3, dim4 = 3, 81, 81, 81
complex_array = read_4d_complex_array_from_file(filename, dim1, dim2, dim3, dim4)
print("Shape of the loaded array:", complex_array.shape)
