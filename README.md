# Towards Practical and Privacy-Preserving CNN Inference Service for Cloud-Based Medical Imaging Analysis: A Homomorphic Encryption-based Approach

## Project Description
This project uses a convolutional neural network (CNN) model hosted on a cloud server to classify private body-related radiographic images. The focus of the project is on cloud-based medical imaging analysis and privacy protection, aiming to ensure the privacy of medical imaging data while efficiently performing medical image classification tasks.

## Environment Dependencies

### Operating System
- Ubuntu 18.04 or later

### Programming Language and Compiler
- C++ 17 standard

### Required Libraries
1. **GMP (GNU Multiple Precision Arithmetic Library)**  
   Used for arbitrary-precision integer and floating-point arithmetic.  
   Required version: 6.1.2  
   Official website: [https://gmplib.org/](https://gmplib.org/)

2. **NTL (Number Theory Library)**  
   Used for efficient polynomial and lattice-based computations.  
   Required version: 11.3.2  
   Official website: [https://libntl.org/](https://libntl.org/)

3. **Other dependencies**
   - Additional dependencies can be installed via `make` or `CMake`, as specified in the project's `Makefile` or `CMakeLists.txt`.

## Installation Guide

1. **Clone the repository**
   ```bash
   git clone https://github.com/username/project-name.git
   ```

2. **Install dependencies**
   - Install GMP and NTL libraries using the package manager:
     ```bash
     sudo apt-get install libgmp-dev libntl-dev
     ```

   - Ensure you have a compiler that supports C++ 17 standard. If not, you can install it with:
     ```bash
     sudo apt-get install g++-17
     ```

3. **Build the project**
   Navigate to the project directory and run the following command to compile:
   ```bash
   make
   ```

## Usage

1. Run the project:
   ```bash
   ./your_project_executable
   ```

2. Provide your private medical imaging for classification:
   The project supports uploading imaging files to the cloud for classification processing. Please refer to the project documentation for API usage details.

## Directory Structure

```
LICENSE                - GNU General Public License v3.0
code.rar               - Contains all the project code
README.md              - Project documentation
```

## Contributing

We welcome community contributions! To contribute to the development of this project, please follow these steps:

1. Fork this repository
2. Create a new branch
3. Commit your changes and create a Pull Request

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).

## Contact

If you have any questions, feel free to reach out to us at [253863586@qq.com].
